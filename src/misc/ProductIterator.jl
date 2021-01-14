"""
    ProductIterator(iters::AbstractArray; inplace=false)

Return an iterator over the product of iterators in `iters`. Each generated element is
an array of the same shape as `iters` whose `i`-th element comes from the `i`-th iterator
in `iters`. The first iterator changes the fastest.
If `inplace` is `true`, each generated element re-uses the same array.

# Examples
```julia
julia> collect(ProductIterator([1:2, 1:3]))
6-element Vector{Vector{Int64}}:
 [1, 1]
 [2, 1]
 [1, 2]
 [2, 2]
 [1, 3]
 [2, 3]

julia> collect(ProductIterator([(1, 2) (3,)]))
2-element Vector{Matrix{Int64}}:
 [1 3]
 [2 3]

julia> p = ProductIterator([(1, 2) (3,)], inplace=true); collect(p)
2-element Vector{Matrix{Int64}}:
 [2 3]
 [2 3]

julia> for x in p; println(x); end
[1 3]
[2 3]
```
"""
struct ProductIterator{T, N}
   iters::Array{T, N}
   inplace::Bool

   ProductIterator(iters::AbstractArray{T, N}; inplace::Bool=false) where {T, N} =
      new{T, N}(iters, inplace)
end

ProductIterator(iter, n::Integer; inplace::Bool=false) =
   ProductIterator(fill(iter, n); inplace=inplace)

Base.eltype(::Type{ProductIterator{T, N}}) where {T, N} = Array{eltype(T), N}

function Base.IteratorSize(::Type{ProductIterator{T, N}}) where {T, N}
   if Base.IteratorSize(T) isa Base.HasShape
      Base.HasLength()
   else
      Base.IteratorSize(T)
   end
end

Base.length(p::ProductIterator) = prod(length, p.iters)

function Base.iterate(p::ProductIterator{T, N}) where {T, N}
   iters = p.iters
   xs = eltype(p)(undef, size(iters))

   x_st = iterate(iters[1])
   x_st === nothing && return nothing
   x, st = x_st
   states = Array{typeof(st), N}(undef, size(iters))
   states[1] = st
   xs[1] = x

   for i = 2:length(iters)
      x_st = iterate(iters[i])
      x_st === nothing && return nothing
      xs[i] = x_st[1]
      states[i] = x_st[2]
   end

   value = p.inplace ? xs : copy(xs)
   value, (xs, states)
end

function Base.iterate(p::ProductIterator, (xs, states))
   iters = p.iters
   n = length(iters)
   for i = 1:n
      x_st = iterate(iters[i], states[i])
      if x_st === nothing
         i == n ? (return nothing) : continue
      end
      for j = 1:i-1
         # all previous entries had reached end-of-iteration, which therefore need
         # to be restarted; we don't do it in the `x_st === nothing` branch to not do
         # useless work when p's iteration is over
         xs[j], states[j] = iterate(iters[j])
                            # !== nothing, otherwise first iteration of p would have
                            # yielded nothing (if p.iters are "well" behaved...)
      end
      xs[i], states[i] = x_st
      break
   end

   value = p.inplace ? xs : copy(xs)
   value, (xs, states)
end
