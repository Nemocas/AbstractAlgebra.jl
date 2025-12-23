#################################################################################
#
#   Factor.jl : Factorization
#
#################################################################################

################################################################################
#
#   Type
#
################################################################################

@doc raw"""
    Fac{T <: RingElement}

Type for factored ring elements. The structure holds a unit of type `T` and is
an iterable collection of `T => Int` pairs for the factors and exponents.

See [`unit(a::Fac)`](@ref), [`evaluate(a::Fac)`](@ref).
"""
mutable struct Fac{T <: RingElement}
   fac::Dict{T, Int}
   unit::T
   arr::Vector{Pair{T, Int}}

   function Fac{T}() where {T <: RingElement}
     return new{T}(Dict{T, Int}())
   end

   function Fac{T}(u::T, d::Dict{T, Int}) where {T <: RingElement}
     return new{T}(d, u)
   end

   function Fac{T}(u::T, a::Vector{Pair{T, Int}}) where {T <: RingElement}
     z = new{T}()
     z.unit = u
     z.arr = a
     return z
   end
end

function Fac(u::T, d::Dict{T, Int}) where {T}
   f = Fac{T}(u, d)
   return f
end

function Fac(u::T, d::Vector{Pair{T, Int}}) where {T}
  return Fac{T}(u, d)
end

function Fac(u::T, d::Vector{Tuple{T, Int}}) where {T}
  return Fac{T}(u, [x[1] => x[2] for x in d])
end

_is_legal(a::Fac) = xor(isdefined(a, :fac), isdefined(a, :arr))

_is_dic(a::Fac) = _is_legal(a) && isdefined(a, :fac)

function Base.:(==)(F1::Fac, F2::Fac)
  error("Equality testing of factorizations is not supported")
end

@doc raw"""
    unit(a::Fac{T}) -> T

Return the unit of the factorization.
"""
unit(a::Fac) = a.unit

#primes(a::Fac) = collect(keys(a.fac))

@doc raw"""
    evaluate(a::Fac{T}) -> T

Multiply out the factorization into a single element.
"""
function evaluate(a::Fac)
   r = a.unit
   if _is_dic(a)
     for (p, e) in a.fac
        r *= p^e
     end
   else
     for (p, e) in a.arr
        r *= p^e
     end
   end
   return r
end

################################################################################
#
#   Syntax sugar
#
################################################################################

@doc raw"""
    in(a, b::Fac)

Test whether $a$ is a factor of $b$.
"""
function Base.in(a, b::Fac{T}) where {T}
   # convert is necessary when T == ZZRingElem, because hash on ZZRingElem
   # doesn't coincide with hash on Integer
   if _is_dic(b)
     convert(T, a) in keys(b.fac)
   else
     return in(a, first.(b.arr))
   end
end

@doc raw"""
    getindex(a::Fac, b) -> Int

If $b$ is a factor of $a$, the corresponding exponent is returned. Otherwise
an error is thrown.
"""
function getindex(a::Fac{T}, b) where {T}
  if _is_dic(a)
    b = convert(T, b)
    if haskey(a.fac, b)
      return a.fac[b]
    else
      error("$b is not a factor of $a")
    end
  else
    i = findfirst(==(b), first.(a.arr))
    i === nothing && error("$b is not a factor of $a")
    return a.arr[i::Int][2]
  end
end

@doc raw"""
    setindex!(a::Fac{T}, c::Int, b::T)

If $b$ is a factor of $a$, the corresponding entry is set to $c$.
"""
function setindex!(a::Fac{T}, c::Int, b::T) where {T}
  if _is_dic(a)
    if haskey(a.fac, b)
      error("$b is already set (to $(a[b]))")
    else
      setindex!(a.fac, c, b)
    end
  else
    error("not supported")
  end
end

################################################################################
#
#   String I/O
#
################################################################################

function expressify(@nospecialize(a::Fac); context = nothing)
   prod = Expr(:call, :cdot)
   if isdefined(a, :unit)
      if !is_one(a.unit)
         push!(prod.args, expressify(a.unit, context = context))
      end
   else
      push!(prod.args, Expr(:call, :*, "[unit not set]"))
   end
   if _is_dic(a)
     it = a.fac
   else
     it = a.arr
   end
   for (p, i) in it
      ep = expressify(p, context = context)
      if isone(i)
         push!(prod.args, ep)
      else
         push!(prod.args, Expr(:call, :^, ep, i))
      end
   end
   return prod
end

@enable_all_show_via_expressify Fac

################################################################################
#
#   Make Factor objects iterable
#
################################################################################

Base.iterate(a::Fac) = _is_dic(a) ? Base.iterate(a.fac) : Base.iterate(a.arr)

Base.iterate(a::Fac, b) = _is_dic(a) ? Base.iterate(a.fac, b) : Base.iterate(a.arr, b)

Base.eltype(::Type{Fac{T}}) where {T} = Pair{T, Int}

@doc raw"""
    length(a::Fac) -> Int

Return the number of factors of $a$, not including the unit.
"""
Base.length(a::Fac) = _is_dic(a) ? Base.length(a.fac) : Base.length(a.arr)
