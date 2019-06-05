###############################################################################
#
#   Type and parent object methods
#
###############################################################################

@doc Markdown.doc"""
    parent_type(::Type{perm{T}}) where T = PermGroup{T}
> Return the type of the parent of a permutation.
"""
parent_type(::Type{perm{T}}) where T = PermGroup{T}

@doc Markdown.doc"""
    elem_type(::Type{PermGroup{T}}) where T = perm{T}
> Return the type of elements of a permutation group.
"""
elem_type(::Type{PermGroup{T}}) where T = perm{T}

@doc Markdown.doc"""
    parent(g::perm{T}) where T = PermGroup
> Return the parent of the permutation `g`.

```jldoctest
julia> G = PermutationGroup(5); g = perm([3,4,5,2,1])
(1,3,5)(2,4)

julia> parent(g) == G
true
```
"""
parent(g::perm{T}) where T = PermGroup(T(length(g.d)))

###############################################################################
#
#   Low-level manipulation
#
###############################################################################

# hash(perm) = 0x0d9939c64ab650ca
Base.hash(g::perm, h::UInt) = xor(hash(g.d, h), 0x0d9939c64ab650ca)

function getindex(g::perm, n::Integer)
   return g.d[n]
end

function setindex!(g::perm, v::Integer, n::Integer)
   g.modified = true
   g.d[n] = v
   return g
end

Base.promote_rule(::Type{perm{I}}, ::Type{perm{J}}) where {I,J} =
   perm{promote_type(I,J)}

convert(::Type{perm{T}}, p::perm) where T = perm(convert(Vector{T}, p.d), false)

convert(::Type{Vector{T}}, p::perm{T}) where {T} = p.d

###############################################################################
#
#   Basic functions
#
###############################################################################

@doc Markdown.doc"""
    parity(g::perm{T}) where T
> Return the parity of the given permutation, i.e. the parity of the number of
> transpositions in any decomposition of `g` into transpositions.
>
> `parity` returns $1$ if the number is odd and $0$ otherwise. `parity` uses
> cycle decomposition of `g` if already available, but will not compute
> it on demand. Since cycle structure is cached in `g` you may call
> `cycles(g)` before calling `parity`.

# Examples:
```jldoctest
julia> g = perm([3,4,1,2,5])
(1,3)(2,4)

julia> parity(g)
0

julia> g = perm([3,4,5,2,1,6])
(1,3,5)(2,4)

julia> parity(g)
1
```
"""
function parity(g::perm{T}) where T
   if isdefined(g, :cycles) && !g.modified
      return T(sum([(length(c)+1)%2 for c in cycles(g)])%2)
   end
   to_visit = trues(size(g.d))
   parity = false
   k = 1
   @inbounds while any(to_visit)
      k = findnext(to_visit, k)
      to_visit[k] = false
      next = g[k]
      while next != k
         parity = !parity
         to_visit[next] = false
         next = g[next]
      end
   end
   return T(parity)
end

@doc Markdown.doc"""
    sign(g::perm{T}) where T
> Return the sign of permutation.
>
> `sign` returns $1$ if `g` is even and $-1$ if `g` is odd. `sign` represents
> the homomorphism from the permutation group to the unit group of $\mathbb{Z}$
> whose kernel is the alternating group.

# Examples:
```jldoctest
julia> g = perm([3,4,1,2,5])
(1,3)(2,4)

julia> sign(g)
1

julia> g = perm([3,4,5,2,1,6])
(1,3,5)(2,4)

julia> sign(g)
-1
```
"""
sign(g::perm{T}) where T = (-one(T))^parity(g)

###############################################################################
#
#   Iterator protocol for CycleDec
#
###############################################################################

function Base.iterate(cd::CycleDec)
   this = cd.cptrs[1]
   next = cd.cptrs[2]
   return (view(cd.ccycles, this:next-1), 2)
end

function Base.iterate(cd::CycleDec, state)
   if state > cd.n
      return nothing
   end

   this = cd.cptrs[state]
   next = cd.cptrs[state + 1]

   return (view(cd.ccycles, this:next-1), state + 1)
end

function Base.getindex(cd::CycleDec, n::Int)
   1 <= n <= length(cd) || throw(BoundsError([cd.cptrs], n+1))
   return cd.ccycles[cd.cptrs[n]:cd.cptrs[n+1]-1]
end

Base.getindex(cd::CycleDec, i::Number) = cd[convert(Int, i)]
Base.getindex(cd::CycleDec, I) = [cd[i] for i in I]

Base.length(cd::CycleDec) = cd.n
Base.lastindex(cd::CycleDec) = cd.n
Base.eltype(::Type{CycleDec{T}}) where T = Vector{T}

function Base.show(io::IO, cd::CycleDec)
   a = [join(c, ",") for c in cd]
   print(io, "Cycle Decomposition: ("*join(a, ")(")*")")
end

@doc Markdown.doc"""
    cycles(g::perm{T}) where T<:Integer
> Decompose permutation `g` into disjoint cycles.
>
> Returns a `CycleDec` object which iterates over disjoint cycles of `g`. The
> ordering of cycles is not guaranteed, and the order within each cycle is
> computed up to a cyclic permutation.
> The cycle decomposition is cached in `g` and used in future computation of
> `permtype`, `parity`, `sign`, `order` and `^` (powering).

# Examples:
```jldoctest
julia> g = perm([3,4,5,2,1,6])
(1,3,5)(2,4)

julia> collect(cycles(g))
3-element Array{Array{Int64,1},1}:
 [1, 3, 5]
 [2, 4]
 [6]
```
"""
function cycles(g::perm{T}) where T<:Integer
   if !isdefined(g, :cycles) || g.modified
      ccycles, cptrs = cycledec(g.d)
      g.cycles = CycleDec{T}(ccycles, cptrs, length(cptrs)-1)
      g.modified = false
   end
   return g.cycles
end

function cycledec(v::Vector{T}) where T<:Integer
  to_visit = trues(size(v))
   ccycles = similar(v) # consecutive cycles entries
   cptrs = [1] # pointers to where cycles start
   # ccycles[cptrs[i], cptrs[i+1]-1] contains i-th cycle

   # expected number of cycles - (overestimation of) the harmonic
   sizehint!(cptrs, 5 + ceil(Int, log(length(v))))
   # cptrs[1] = one(T)

   k = 1
   i = 1

   while any(to_visit)
      k = findnext(to_visit, k)
      to_visit[k] = false
      next = v[k]

      ccycles[i] = T(k)
      i += 1
      while next != k
         ccycles[i] = next
         to_visit[next] = false
         next = v[next]
         i += 1
      end
      push!(cptrs, i)
   end
   return ccycles, cptrs
end

@doc Markdown.doc"""
    permtype(g::perm)
> Return the type of permutation `g`, i.e. lengths of disjoint cycles in cycle
> decomposition of `g`.
>
> The lengths are sorted in decreasing order by default. `permtype(g)` fully
> determines the conjugacy class of `g`.

# Examples:
```jldoctest
julia> g = perm([3,4,5,2,1,6])
(1,3,5)(2,4)

julia> permtype(g)
3-element Array{Int64,1}:
 3
 2
 1

julia> G = PermGroup(5); e = parent(g)()
()

julia> permtype(e)
6-element Array{Int64,1}:
 1
 1
 1
 1
 1
 1
```
"""
permtype(g::perm) = sort(diff(cycles(g).cptrs), rev=true)

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::PermGroup)
   print(io, "Permutation group over $(G.n) elements")
end

mutable struct PermDisplayStyle
   format::Symbol
end

const _permdisplaystyle = PermDisplayStyle(:cycles)

@doc Markdown.doc"""
    setpermstyle(format::Symbol)
> Select the style in which permutations are displayed (in REPL or in general
> as string). This can be either
> * `:array` - as vectors of integers whose $n$-th position represents the
> value at $n$), or
> * `:cycles` - as, more familiar for mathematicians, decomposition into
> disjoint cycles, where the value at $n$ is represented by the entry
> immediately following $n$ in a cycle (the default).

> The difference is purely esthetical.

# Examples:
```jldoctest
julia> Generic.setpermstyle(:array)
:array

julia> perm([2,3,1,5,4])
[2, 3, 1, 5, 4]

julia> Generic.setpermstyle(:cycles)
:cycles

julia> perm([2,3,1,5,4])
(1,2,3)(4,5)
```
"""
function setpermstyle(format::Symbol)
   if format in (:array, :cycles)
      _permdisplaystyle.format = format
   else
      throw("Permutations can be displayed only as :array or :cycles.")
   end
   return format
end

function show(io::IO, g::perm)
   if _permdisplaystyle.format == :array
      print(io, "[" * join(g.d, ", ") * "]")
   elseif _permdisplaystyle.format == :cycles
      cd = cycles(g)
      if g == parent(g)()
         print(io, "()")
      else
         print(io, join(["("*join(c, ",")*")" for c in cd if length(c)>1],""))
      end
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

@doc Markdown.doc"""
    ==(g::perm, h::perm)
> Return `true` if permutations are equal, otherwise return `false`.
>
> Permutations parametrized by different integer types are considered equal if
> they define the same permutation in the abstract permutation group.

# Examples:
```
julia> g = perm(Int8[2,3,1])
(1,2,3)

julia> h = perm"(3,1,2)"
(1,2,3)

julia> g == h
true
```
"""
==(g::perm, h::perm) = g.d == h.d

@doc Markdown.doc"""
    ==(G::PermGroup, H::PermGroup)
> Return `true` if permutation groups are equal, otherwise return `false`.
>
> Permutation groups on the same number of letters, but parametrized
> by different integer types are considered different.

# Examples:
```
julia> G = PermGroup(UInt(5))
Permutation group over 5 elements

julia> H = PermGroup(5)
Permutation group over 5 elements

julia> G == H
false
```
"""
==(G::PermGroup, H::PermGroup) = typeof(G) == typeof(H) && G.n == H.n

###############################################################################
#
#   Binary operators
#
###############################################################################
function mul!(out::perm, g::perm, h::perm)
   @inbounds for i in eachindex(out.d)
      out[i] = h[g[i]]
   end
   return out
end

@doc Markdown.doc"""
    *(g::perm{T}, h::perm{T}) where T
> Return the composition ``h ∘ g`` of two permutations.
>
> This corresponds to the action of permutation group on the set `[1..n]`
> **on the right** and follows the convention of GAP.
>
> If `g` and `h` are parametrized by different types, the result is promoted
> accordingly.

# Examples:
```jldoctest
julia> perm([2,3,1,4])*perm([1,3,4,2]) # (1,2,3)*(2,3,4)
(1,3)(2,4)
```
"""
function *(g::perm{T}, h::perm{T}) where T
   res = perm(similar(g.d), false)
   return mul!(res, g, h)
end

*(g::perm{S}, h::perm{T}) where {S,T} = *(promote(g,h)...)

@doc Markdown.doc"""
    ^(g::perm{T}, n::Integer) where T
> Return the $n$-th power of a permutation `g`.
>
> By default `g^n` is computed by cycle decomposition of `g` if `n > 3`.
> `Generic.power_by_squaring` provides a different method for powering which
> may or may not be faster, depending on the particuar case. Due to caching of
> the cycle structure, repeated powering of `g` will be faster with the default
> method.

# Examples:
```jldoctest
julia> g = perm([2,3,4,5,1])
(1,2,3,4,5)

julia> g^3
(1,4,2,5,3)

julia> g^5
()
```
"""
function ^(g::perm{T}, n::Integer) where T
   if n < 0
      return inv(g)^-n
   elseif n == 0
      return perm(T(length(g.d)))
   elseif n == 1
      return deepcopy(g)
   elseif n == 2
      return perm(g.d[g.d], false)
   elseif n == 3
      return perm(g.d[g.d[g.d]], false)
   else
      new_perm = similar(g.d)

      @inbounds for cycle in cycles(g)
         l = length(cycle)
         k = n % l
         for (idx,j) in enumerate(cycle)
            idx += k
            idx = (idx > l ? idx-l : idx)
            new_perm[j] = cycle[idx]
         end
      end
      p = perm(new_perm, false)
      return p
   end
end

function power_by_squaring(g::perm{I}, n::Integer) where {I}
   if n < 0
      return inv(g)^-n
   elseif n == 0
      return perm(T(length(g.d)))
   elseif n == 1
      return deepcopy(g)
   elseif n == 2
      return perm(g.d[g.d], false)
   elseif n == 3
      return perm(g.d[g.d[g.d]], false)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & n) == 0
         bit >>= 1
      end
      cache1 = deepcopy(g.d)
      cache2 = deepcopy(g.d)
      bit >>= 1
      while bit != 0
         cache2 = cache1[cache1]
         cache1 = cache2
         if (UInt(bit) & n) != 0
            cache1 = cache1[g.d]
         end
         bit >>= 1
      end
      return perm(cache1, false)
   end
end

###############################################################################
#
#   Inversion
#
###############################################################################

@doc Markdown.doc"""
    inv(g::perm)
> Return the inverse of the given permutation, i.e. the permuation $g^{-1}$
> such that $g ∘ g^{-1} = g^{-1} ∘ g$ is the identity permutation.
"""
function inv(g::perm)
   d = similar(g.d)
   @inbounds for i in 1:length(d)
      d[g[i]] = i
   end
   return perm(d, false)
end

# TODO: See M. Robertson, Inverting Permutations In Place
# n+O(log^2 n) space, O(n*log n) time
function inv!(a::perm)
   d = similar(a.d)
   @inbounds for i in 1:length(d)
      d[a[i]] = i
   end
   a.d = d
   a.modified = true
   return a
end

###############################################################################
#
#   Iterating over all permutations
#
###############################################################################

@inline function Base.iterate(A::AllPerms{<: Integer})
   return A.elts, 1
end

@inline function Base.iterate(A::AllPerms{<: Integer}, count)
   if count >= A.all
     # Reset the iterator to make it iterable again
     for i in 1:length(A.c)
       A.c[i] = 1
     end
     return nothing
   end

   k = 0
   n = 1

   @inbounds while true
      if A.c[n] < n
         k = ifelse(isodd(n), 1, A.c[n])
         A.elts[k], A.elts[n] = A.elts[n], A.elts[k]
         A.c[n] += 1
         return A.elts, count + 1
      else
         A.c[n] = 1
         n += 1
      end
   end
end

#Base.start(A::AllPerms{<:Integer}) = 0
#
#function Base.next(A::AllPerms, count)
#    count = nextperm(A, count)
#    return A.elts, count
#end
#
#Base.done(A::AllPerms, count) = count >= A.all

Base.eltype(::Type{AllPerms{T}}) where T<:Integer = perm{T}

Base.length(A::AllPerms) = A.all

#function nextperm(A::AllPerms{<:Integer}, count)
#   if count == 0
#        return count+1
#   end
#
#   k = 0
#   n = 1
#
#   while true
#      if A.c[n] < n
#         k = ifelse(isodd(n), 1, A.c[n])
#         A.elts[k], A.elts[n] = A.elts[n], A.elts[k]
#         A.c[n] += 1
#         return count+1
#      else
#         A.c[n] = 1
#         n += 1
#      end
#   end
#end

@doc Markdown.doc"""
    Generic.elements!(G::PermGroup)
> Return an unsafe iterator over all permutations in `G`. Only one permutation
> is allocated and then modified in-place using the non-recursive
> [Heaps algorithm](https://en.wikipedia.org/wiki/Heap's_algorithm).
>
> Note: you need to explicitely copy permutations intended to be stored or
> modified.

# Examples:
```jldoctest
julia> elts = Generic.elements!(PermGroup(5));

julia> length(elts)
120

julia> for p in Generic.elements!(PermGroup(3))
         println(p)
       end
()
(1,2)
(1,3,2)
(2,3)
(1,2,3)
(1,3)

julia> A = collect(Generic.elements!(PermGroup(3))); A
6-element Array{AbstractAlgebra.Generic.perm{Int64},1}:
 (1,3)
 (1,3)
 (1,3)
 (1,3)
 (1,3)
 (1,3)

julia> unique(A)
1-element Array{AbstractAlgebra.Generic.perm{Int64},1}:
 (1,3)
```
"""
elements!(G::PermGroup)= (p for p in AllPerms(G.n))

function Base.iterate(G::PermGroup)
  A = AllPerms(G.n)
  a, b = iterate(A)
  return deepcopy(A.elts), (A, b)
end

function Base.iterate(G::PermGroup, S)
  A = S[1]
  c = S[2]

  s = iterate(A, c)

  if s === nothing
    return nothing
  end

  return deepcopy(s[1]), (A, s[2])
end

Base.eltype(::Type{PermGroup{T}}) where T = perm{T}

Base.length(G::PermGroup) = order(G)

###############################################################################
#
#   Misc
#
###############################################################################

@doc Markdown.doc"""
    order(G::PermGroup) -> BigInt
> Return the order of the full permutation group as `BigInt`.
"""
order(G::PermGroup) = order(BigInt, G)
order(::Type{T}, G::PermGroup) where T = factorial(T(G.n))

@doc Markdown.doc"""
    order(a::perm) -> BigInt
> Return the order of permutation `a` as `BigInt`.
>
> If you are sure that computation over `T` (or its `Int` promotion) will not
> overflow you may use the method `order(T::Type, a::perm)` which bypasses
> computation with BigInts and returns `promote(T, Int)`.
"""
order(a::perm) = order(BigInt, a)
function order(::Type{T}, a::perm) where T
   TT = promote_type(T, Int)
   return reduce(lcm, diff(cycles(a).cptrs), init = one(TT))
end

@doc Markdown.doc"""
    matrix_repr(a::perm{T}) where T
> Return the permutation matrix as sparse matrix representing `a` via natural
> embedding of the permutation group into general linear group over $\mathbb{Z}$.

# Examples:
```jldoctest
julia> p = perm([2,3,1])
(1,2,3)

julia> matrix_repr(p)
3×3 SparseMatrixCSC{Int64,Int64} with 3 stored entries:
  [3, 1]  =  1
  [1, 2]  =  1
  [2, 3]  =  1

julia> Array(ans)
3×3 Array{Int64,2}:
 0  1  0
 0  0  1
 1  0  0
```
"""
matrix_repr(a::perm{T}) where T = sparse(collect(T, 1:length(a.d)), a.d, ones(T,length(a.d)))

@doc Markdown.doc"""
    emb!(result::perm, p::perm, V)
> Embed permutation `p` into permutation `result` on the indices given by `V`.
>
> This corresponds to the natural embedding of $S_k$ into $S_n$ as the
> subgroup permuting points indexed by `V`.

# Examples:
```jldoctest
julia> p = perm([2,1,4,3])
(1,2)(3,4)

julia> Generic.emb!(perm(collect(1:5)), p, [3,1,4,5])
(1,3)(4,5)
```
"""
function emb!(result::perm, p::perm, V)
   result.d[V] = (result.d[V])[p.d]
   return result
end

@doc Markdown.doc"""
    emb(G::PermGroup, V::Vector{Int}, check::Bool=true)
> Return the natural embedding of a permutation group into `G` as the
> subgroup permuting points indexed by `V`.

# Examples:
```jldoctest
julia> p = perm([2,3,1])
(1,2,3)

julia> f = Generic.emb(PermGroup(5), [3,2,5]);

julia> f(p)
(2,5,3)
```
"""
function emb(G::PermGroup, V::Vector{Int}, check::Bool=true)
   if check
      @assert length(Base.Set(V)) == length(V)
      @assert all(V .<= G.n)
   end
   return p -> Generic.emb!(G(), p, V)
end

@doc Markdown.doc"""
    rand(G::PermGroup{T}) where {T}
> Return a random permutation from `G`.
"""
rand(G::PermGroup{T}) where {T} = perm(randperm!(Vector{T}(undef, G.n)), false)

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

(G::PermGroup)() = perm(G.n)

function (G::PermGroup{T})(a::Vector{T}, check::Bool=true) where T<:Integer
   if check
      G.n == length(a) || throw("Cannot coerce to $G: lengths differ")
   end
   return perm(a, check)
end

function (G::PermGroup{T})(a::Vector{S}, check::Bool=true) where {S<:Integer,T}
   return G(convert(Vector{T}, a), check)
end

function (G::PermGroup{T})(p::perm{S}, check::Bool=true) where {S<:Integer, T}
   if parent(p) == G
      return p
   end
   return G(convert(Vector{T}, p.d), check)
end

function (G::PermGroup)(str::String, check::Bool=true)
   return G(cycledec(parse_cycles(str)..., G.n), check)
end

function (G::PermGroup{T})(cdec::CycleDec{T}, check::Bool=true) where T
   if check
      length(cdec.ccycles) == G.n || throw("Can not coerce to $G: lengths differ")
   end

   elt = perm(G.n)
   for c in cdec
      for i in 1:length(c)-1
         elt[c[i]] = c[i+1]
      end
      elt[c[end]] = c[1]
   end

   elt.cycles = cdec
   return elt
end

###############################################################################
#
#   Parsing strings/GAP output
#
###############################################################################

function parse_cycles(str::AbstractString)
   ccycles = Int[]
   cptrs = Int[1]
   if occursin(r"\n|\s| ", str)
      str = replace(str, r"\n|\s| " => "")
   end
   cycle_regex = r"\(\d+(,\d+)*\)?"
   for cycle_str in (m.match for m = eachmatch(cycle_regex, str))
      cycle = [parse(Int, a) for a in split(cycle_str[2:end-1], ",")]
      append!(ccycles, cycle)
      push!(cptrs, cptrs[end]+length(cycle))
   end
   return ccycles, cptrs
end

function cycledec(ccycles::Vector{Int}, cptrs::Vector{Int}, n::T,
   check::Bool=true) where T
   if check
      if length(ccycles) != 0
         maximum(ccycles) <= n || throw("elts in $ccycles larger than $n")
      end
      length(Set(ccycles)) == length(ccycles) || throw("Non-unique entries in $ccycles")
   end

   if length(ccycles) != n
      sizehint!(ccycles, n)
      to_append = filter(x -> !(x in ccycles), 1:n)
      l = length(ccycles)
      append!(cptrs, l+2:l+length(to_append)+1)
      append!(ccycles, to_append)
   end

   return CycleDec{T}(ccycles, cptrs, length(cptrs)-1)
end

@doc Markdown.doc"""
    perm"..."
> String macro to parse disjoint cycles into `perm{Int}`.
>
> Strings for the output of GAP could be copied directly into `perm"..."`.
> Cycles of length $1$ are not necessary, but could be included. A permutation
> of the minimal support is constructed, i.e. the maximal $n$ in the
> decomposition determines the parent group $S_n$.

# Examples:
```jldoctest
julia> p = perm"(1,3)(2,4)"
(1,3)(2,4)

julia> typeof(p)
AbstractAlgebra.Generic.perm{Int64}

julia> parent(p) == PermutationGroup(4)
true

julia> p = perm"(1,3)(2,4)(10)"
(1,3)(2,4)

julia> parent(p) == PermutationGroup(10)
true
```
"""
macro perm_str(s)
   c, p = parse_cycles(s)
   if length(c) == 0
      n = 1
   else
      n = maximum(c)
   end
   cdec = cycledec(c, p, n)
   return PermGroup(cdec.cptrs[end]-1)(cdec)
end

###############################################################################
#
#   PermGroup constructor
#
###############################################################################

# handled by inner constructors

##############################################################################
#
#   Irreducible Characters
#
##############################################################################

const _charvalsTable = Dict{Tuple{BitVector,Vector{Int}}, Int}()
const _charvalsTableBig = Dict{Tuple{BitVector,Vector{Int}}, BigInt}()

@doc Markdown.doc"""
    character(lambda::Partition)
> Return the $\lambda$-th irreducible character of permutation group on
> `sum(lambda)` symbols. The returned character function is of the following signature:
> > `chi(p::perm[, check::Bool=true]) -> BigInt`
> The function checks (if `p` belongs to the appropriate group) can be switched
> off by calling `chi(p, false)`. The values computed by $\chi$ are cached in
> look-up table.
>
> The computation follows the Murnaghan-Nakayama formula:
> $$\chi_\lambda(\sigma) = \sum_{\text{rimhook }\xi\subset \lambda}(-1)^{ll(\lambda\backslash\xi)} \chi_{\lambda \backslash\xi}(\tilde\sigma)$$
> where $\lambda\backslash\xi$ denotes the skew diagram of $\lambda$ with $\xi$
> removed, $ll$ denotes the leg-length (i.e. number of rows - 1) and
> $\tilde\sigma$ is permutation obtained from $\sigma$ by the removal of the
> longest cycle.
>
> For more details see e.g. Chapter 2.8 of *Group Theory and Physics* by
> S.Sternberg.

# Examples
```jldoctest
julia> G = PermutationGroup(4)
Permutation group over 4 elements

julia> chi = character(Partition([3,1])) # character of the regular representation
(::char) (generic function with 2 methods)

julia> chi(G())
3

julia> chi(perm"(1,3)(2,4)")
-1
```
"""
function character(lambda::Partition)
   R = partitionseq(lambda)

   function char(p::perm, check::Bool=true)
      if check
         sum(lambda) == length(p.d) || throw(ArgumentError("Can't evaluate character on $p : lengths differ."))
      end
      return MN1inner(R, Partition(permtype(p)), 1, _charvalsTableBig)
   end

   return char
end

@doc Markdown.doc"""
    character(lambda::Partition, p::perm, check::Bool=true) -> BigInt
> Returns the value of `lambda`-th irreducible character of the permutation
> group on permutation `p`.
"""
function character(lambda::Partition, p::perm, check::Bool=true)
   if check
      sum(lambda) == length(p.d) || throw("lambda-th irreducible character can be evaluated only on permutations of length $(sum(lambda)).")
   end
   return character(BigInt, lambda, Partition(permtype(p)))
end

function character(::Type{T}, lambda::Partition, p::perm) where T <: Integer
   return character(T, lambda, Partition(permtype(p)))
end

@doc Markdown.doc"""
    character(lambda::Partition, mu::Partition, check::Bool=true) -> BigInt
> Returns the value of `lambda-th` irreducible character on the conjugacy class
> represented by partition `mu`.
"""
function character(lambda::Partition, mu::Partition, check::Bool=true)
   if check
      sum(lambda) == sum(mu) || throw("Cannot evaluate $lambda on the conjugacy class of $mu: lengths differ.")
   end
   return character(BigInt, lambda, mu)
end

function character(::Type{BigInt}, lambda::Partition, mu::Partition)
   return MN1inner(partitionseq(lambda), mu, 1, _charvalsTableBig)
end

function character(::Type{T}, lambda::Partition, mu::Partition) where T<:Union{Signed, Unsigned}
   return MN1inner(partitionseq(lambda), mu, 1, _charvalsTable)
end
