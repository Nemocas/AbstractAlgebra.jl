export PermGroup, perm, parity, elements, cycles, character, setpermstyle,
       order

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent_type(::Type{perm{T}}) where T = PermGroup{T}

elem_type(::Type{PermGroup{T}}) where T = perm{T}

doc"""
    parent(a::perm)
> Return the parent of the given permutation group element.
"""
parent(a::perm) = a.parent

###############################################################################
#
#   Low-level manipulation
#
###############################################################################

function deepcopy_internal(a::perm, dict::ObjectIdDict)
   G = parent(a)
   return G(deepcopy(a.d), false)
end

function Base.hash(a::perm, h::UInt)
   b = 0x595dee0e71d271d0%UInt
   for i in 1:length(a.d)
      b = xor(b, xor(hash(a[i], h), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

function getindex(a::perm{T}, n::S) where {T<:Integer, S<:Integer}
   return a.d[n]
end

function setindex!(a::perm{T}, v::T, n::S) where {T<:Integer, S<:Integer}
   a.d[n] = v
end

Base.promote_rule(::Type{perm{I}}, ::Type{perm{J}}) where {I,J} =
   perm{promote_type(I,J)}

convert(::Type{perm{I}}, p::perm{J}) where {I,J} = (PermGroup(I(parent(p).n)))(p.d, false)

###############################################################################
#
#   Basic functions
#
###############################################################################

doc"""
    parity(a::perm)
> Return the parity of the given permutation, i.e. the parity of the number of
> transpositions that compose it. The function returns $1$ if the parity is odd
> and $0$ otherwise. By default `parity` will uses the cycle decomposition if
> it is already available, but will not compute it on demand. If You intend to
> use parity, or the cycle decomposition of a permutation later You may force
> `parity` to compute the cycle structure by calling
> `parity(a, Val{:cycles})``.
"""
# TODO: 2x slower than Flint
function parity(a::perm{T}) where T
   if isdefined(a, :cycles)
      return T(sum([(length(c)+1)%2 for c in cycles(a)])%2)
   end
   to_visit = trues(a.d)
   parity = false
   k = 1
   @inbounds while any(to_visit)
      k = findnext(to_visit, k)
      to_visit[k] = false
      next = a[k]
      while next != k
         parity = !parity
         to_visit[next] = false
         next = a[next]
      end
   end
   return T(parity)
end

function parity(a::perm, ::Type{Val{:cycles}})
   cycles(a)
   return parity(a)
end

doc"""
    sign(a::perm)
> Returns the sign of the given permutations, i.e. `1` if `a` is even and `-1`
> if `a` is odd.
"""
sign(a::perm{T}) where T = (-one(T))^parity(a)

function sign(a::perm, ::Type{Val{:cycles}})
   cycles(a)
   return sign(a)
end

###############################################################################
#
#   Iterator protocol for CycleDec
#
###############################################################################

Base.start(cd::CycleDec)= 1
Base.done(cd::CycleDec, state::Int) = state > cd.n

function Base.next(cd::CycleDec, state::Int)
   this = cd.cptrs[state]
   next = cd.cptrs[state+1]
   return (view(cd.ccycles, this:next-1), state+1)
end

Base.IndexStyle(::Type{CycleDec}) = IndexLinear()

function Base.getindex(cd::CycleDec, n::Int)
   1 <= n <= length(cd) || throw(BoundsError([cd.cptrs], n+1))
   return cd.ccycles[cd.cptrs[n]:cd.cptrs[n+1]-1]
end

Base.getindex(cd::CycleDec, i::Number) = S[convert(Int, i)]
Base.getindex(cd::CycleDec, I) = [cd[i] for i in I]

Base.length(cd::CycleDec) = cd.n
Base.endof(cd::CycleDec) = cd.n
Base.eltype(::Type{CycleDec{T}}) where T = Vector{T}

function Base.show(io::IO, cd::CycleDec)
   a = [join(c, ",") for c in cd]
   print(io, "Cycle Decomposition: ("*join(a, ")(")*")")
end

doc"""
    cycles(a::perm)
> Decomposes permutation into disjoint cycles.
"""
function cycles(a::perm{T}) where T<:Integer
   if !isdefined(a, :cycles)
      a.cycles = CycleDec(cycledec(a.d)...)
   end
   return a.cycles
end

function cycledec(v::Vector{T}) where T<:Integer
   to_visit = trues(v)
   ccycles = similar(v) # consecutive cycles entries
   cptrs = Vector{Int}() # pointers to where cycles start
   # ccycles[cptrs[i], cptrs[i+1]-1] contains i-th cycle

   # expected number of cycles - (overestimation of) the harmonic
   sizehint!(cptrs, 5 + ceil(Int, log(length(v))))
   # cptrs[1] = one(T)

   k = 1
   i = 1

   while any(to_visit)
      push!(cptrs, i)
      k = findnext(to_visit, k)
      to_visit[k] = false
      next = v[k]

      ccycles[i] = T(k)
      i += 1
      while next != k
         ccycles[i] = next
         to_visit[next] = false
         next = v[next]
         i += one(T)
      end

   end
   return ccycles, cptrs
end

doc"""
    permtype(a::perm, rev=true)
> Returns the type of permutation `a`, i.e. lengths of disjoint cycles building
> `a`. This fully determines the conjugacy class of `a`. The lengths are sorted
> in reverse order by default.
"""
permtype(a::perm{T}) where T = sort(diff(cycles(a).cptrs), rev=true)

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

const _perm_display_style = PermDisplayStyle(:array)

doc"""
    setpermstyle(format::Symbol)
> Nemo can display (in REPL or in general as string) permutations by either
> array of integers whose `n`-th position represents value at `n`, or as
> disjoint cycles, where the value of permutation at `n` is represented as the
> entry immediately following `n` in the cycle.
> The style can be switched by calling `setpermstyle` with `:array` or
> `:cycles` argument.
"""
function setpermstyle(format::Symbol)
   format in (:array, :cycles) || throw("Permutations can be displayed
   only as :array or :cycles.")
   _perm_display_style.format = format
   return format
end

function show(io::IO, a::perm)
   if _perm_display_style.format == :array
      print(io, "[" * join(a.d, ", ") * "]")
   elseif _perm_display_style.format == :cycles
      if a == parent(a)()
         print(io, "()")
      else
         print(io, join(["("*join(c, ",")*")" for c in cycles(a) if
            length(c)>1],""))
      end
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

doc"""
    ==(a::perm, b::perm)
> Return `true` if the given permutations are equal, otherwise return `false`.
> Permutations parametrized by different integer types are considered equal if
> they define the same permutation in the abstract permutation group.
"""
==(a::perm, b::perm) = a.d == b.d

###############################################################################
#
#   Binary operators
#
###############################################################################

doc"""
    *(a::perm, b::perm)
> Return the composition of the two permutations, i.e. $a\circ b$. In other
> words, the permutation corresponding to applying $b$ first, then $a$, is
> returned. If `a` and `b` are parametrized by different types, the result is
> promoted accordingly.
"""
function *(a::perm{T}, b::perm{T}) where T
   d = similar(a.d)
   @inbounds for i in 1:length(d)
      d[i] = a[b[i]]
   end
   return parent(a)(d, false)
end

*(a::perm{S}, b::perm{T}) where {S,T} = *(promote(a,b)...)

doc"""
    ^(a::perm, n::Int)
> Return the `n`-th power of a permutation `a`. By default `a^n` is computed
> by cycle decomposition of `a`. `Generic.power_by_squaring` provides a
> different method for powering which may or may not be faster, depending on the
> particuar case. Due to caching of the cycle structure, repeated powering of
> `a` should be faster with the default method. NOTE: cycle structure is not
> computed for `n<4`.
"""
function ^(a::perm{T}, n::Integer) where T
   if n < 0
      return inv(a)^-n
   elseif n == 0
      return parent(a)()
   elseif n == 1
      return deepcopy(a)
   elseif n == 2
      return parent(a)(a.d[a.d], false)
   elseif n == 3
      return parent(a)(a.d[a.d[a.d]], false)
   else
      new_perm = similar(a.d)

      @inbounds for cycle in cycles(a)
         l = length(cycle)
         k = n % l
         for (idx,j) in enumerate(cycle)
            idx += k
            idx = (idx > l ? idx%l : idx)
            new_perm[j] = cycle[idx]
         end
      end
      p = parent(a)(new_perm, false)
      return p
   end
end

function power_by_squaring(a::perm{I}, n::Integer) where {I}
   if n < 0
      return inv(a)^-n
   elseif n == 0
      return parent(a)()
   elseif n == 1
      return deepcopy(a)
   elseif n == 2
      return parent(a)(a.d[a.d], false)
   elseif n == 3
      return parent(a)(a.d[a.d[a.d]], false)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & n) == 0
         bit >>= 1
      end
      cache1 = deepcopy(a.d)
      cache2 = deepcopy(a.d)
      bit >>= 1
      while bit != 0
         cache2 = cache1[cache1]
         cache1 = cache2
         if (UInt(bit) & n) != 0
            cache1 = cache1[a.d]
         end
         bit >>= 1
      end
      return parent(a)(cache1, false)
   end
end

###############################################################################
#
#   Inversion
#
###############################################################################

doc"""
    inv(a::perm)
> Return the inverse of the given permutation, i.e. the permuation $a^{-1}$
> such that $a\circ a^{-1} = a^{-1}\circ a$ is the identity permutation.
"""
function inv(a::perm)
   d = similar(a.d)
   @inbounds for i in 1:length(d)
      d[a[i]] = i
   end
   return parent(a)(d, false)
end

# TODO: can we do that in place??
function inv!(a::perm)
   d = similar(a.d)
   @inbounds for i in 1:length(d)
      d[a[i]] = i
   end
   a.d = d
   return a
end

###############################################################################
#
#   Iterating over all permutations
#
###############################################################################

Base.start(A::AllPerms{T}) where T<:Integer = (collect(T, 1:A.n), one(T), one(T), ones(T, A.n))
Base.next(A::AllPerms, state) = all_perms(state...)
Base.done(A::AllPerms, state) = state[2] > A.all
Base.eltype(::Type{AllPerms{T}}) where T<:Integer = Vector{T}
length(A::AllPerms) = A.all

function all_perms(elts, counter, i, c)
   if counter == 1
      return (copy(elts), (elts, counter+1, i, c))
   end
   n = length(elts)
   @inbounds while i <= n
      if c[i] < i
         if isodd(i)
            elts[1], elts[i] = elts[i], elts[1]
         else
            elts[c[i]], elts[i] = elts[i], elts[c[i]]
         end
         c[i] += 1
         i = 1
         return (copy(elts), (elts, counter+1, i, c))
      else
         c[i] = 1
         i += 1
      end
   end
end

doc"""
    elements(G::PermGroup)
> Returns an iterator over all elements in the group $G$. You may use
> `collect(elements(G))` to get an array of all elements.
"""
elements(G::PermGroup) = (G(p, false) for p in AllPerms(G.n))

###############################################################################
#
#   Misc
#
###############################################################################

doc"""
    order(G::PermGroup)
> Returns the order of the full permutation group.
"""
order(G::PermGroup) = order(BigInt, G)
order(::Type{T}, G::PermGroup) where T = factorial(T(G.n))

doc"""
    order(a::perm, outputType::Type)
> Returns the order of permutation `a` as `T`, or its widening.
"""
order(a::perm) = order(BigInt, a)
order(::Type{T}, a::perm) where T = reduce(lcm, one(T), diff(cycles(a).cptrs))

doc"""
    matrix_repr(a::perm)
> Return the permutation matrix representing `a` via natural embedding of the
> permutation group into general linear group over ZZ
"""
function matrix_repr(a::perm)
   A = eye(Int, parent(a).n)
   return A[a.d,:]
end

doc"""
    emb!(result::perm, p::perm, V::Vector{Int})
> Embedds permutation `p` into permutation `result` on the indices given by `V`. This
> corresponds to natural embedding of $S_k$ into $S_n$ as the subgroup
> permuting points indexed by `V`.
"""
function emb!(result::perm, p::perm, V::Vector{Int})
    result.d[V] = (result.d[V])[p.d]
    return result
end

doc"""
    emb(G::PermGroup, V::Vector{Int})
> Returns the natural embedding of a permutation group into `G` as the
> subgroup permuting points indexed by `V`.
"""
function emb(G::PermGroup, V::Vector{Int}, check::Bool=true)
   if check
      @assert length(Base.Set(V)) == length(V)
      @assert all(V .<= G.n)
   end
   return p -> Nemo.emb!(G(), p, V)
end

doc"""
    rand(G::PermGroup)
> Returns a random element from group `G`.
"""
rand(G::PermGroup) = G(randperm(G.n), false)

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (G::PermGroup)()
   z = perm(G.n)
   z.parent = G
   return z
end

function (G::PermGroup{T})(a::Vector{T}, check::Bool=true) where T<:Integer
   if check
      length(a) != G.n && error("Unable to coerce to permutation: lengths differ")
      Base.Set(a) != Base.Set(1:length(a)) && error("Unable to coerce to
         permutation: non-unique elements in array")
   end
   z = perm(a)
   z.parent = G
   return z
end

function (G::PermGroup{T})(a::Vector{S}, check::Bool=true) where {S<:Integer,T}
   return G(convert(Vector{T}, a), check)
end

function (G::PermGroup{T})(p::perm{S}, check::Bool=true) where {S<:Integer, T}
   return G(convert(Vector{T}, p.d), check)
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

doc"""
    character(lambda::Partition)
> Returns the $\lambda$-th irreducible character of permutation group on
> `lambda.n` generators. The returned character function
>     `chi(p::perm, check::Bool=true)`
> can be evaluated on a permutation `p::perm`. Checking in $\chi$ if `p`
> belongs to the appropriate group can be switched off by calling
> `chi(p, false)`. The values computed by $\chi$ are stored in look-up table.
>
> The computation follows the Murnaghan-Nakayama formula:
> $$\chi_\lambda(\sigma) = \sum_{\text{rimhook }\xi\subset \lambda}(-1)^{ll(\lambda\backslash\xi)} \chi_{\lambda \backslash\xi}(\tilde\sigma)$$
> where $\lambda\backslash\xi$ denotes partition associated with Young Diagram
> of $\lambda$ with $\xi$ removed, $ll$ denotes the leg-length (i.e. number of
> rows - 1) and $\tilde\sigma$ is permutation obtained from $\sigma$ by the
> removal of the longest cycle.
>
> For more details see e.g. Chapter 2.8 of _Group Theory and Physics_ by
> S.Sternberg.
"""
function character(lambda::Partition)
   R = partitionseq(lambda)

   function char(p::perm, check::Bool=true)
      if check
         lambda.n == length(p.d) || throw("Can't evaluate character on $p : lengths differ.")
      end
      return MN1inner(R, Partition(permtype(p)), 1, _charvalsTableBig)
   end

   return char
end

doc"""
    character(lambda::Partition, p::perm, check::Bool=true)
> Returns the value of `lambda`-th irreducible character on permutation `p`.
"""
function character(lambda::Partition, p::perm{T}, check::Bool=true) where {T}
   if check
      lambda.n == length(p.d) || throw("lambda-th irreducible character can be evaluated only on permutations of length $(lambda.n).")
   end

   return MN1inner(partitionseq(lambda), Partition(permtype(p)), 1, _charvalsTableBig)
end

function character(::Type{T}, lambda::Partition, p::perm) where T <: Int
   return MN1inner(partitionseq(lambda), Partition(permtype(p)), 1, _charvalsTable)
end


doc"""
    character(lambda::Partition, mu::Partition)
> Returns the value of `lambda-th` irreducible character on the conjugacy class
> represented by partition `mu`. Values of characters computed by this method
> are not cached _globally_.
"""
character(lambda::Partition, mu::Partition) = MN1inner(partitionseq(lambda), mu, 1, Dict{Tuple{BitVector,Vector{Int}}, BigInt}())
