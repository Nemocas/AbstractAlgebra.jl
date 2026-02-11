###############################################################################
#
#   Module.jl : Functionality for modules over Euclidean domains
#
###############################################################################

###############################################################################
#
#   Basic manipulation
#
###############################################################################

Base.eltype(T::Type{<:FPModule{<:FinFieldElem}}) = elem_type(T)

function zero(M::FPModule{T}) where T <: RingElement
   R = base_ring(M)
   return M(zero_matrix(R, 1, ngens(M)))
end

function iszero(v::FPModuleElem{T}) where T <: RingElement
   return iszero(Generic._matrix(v))
end

function check_parent(M::FPModule{T}, N::FPModule{T}) where T <: RingElement
   base_ring(M) !== base_ring(N) && error("Incompatible modules")
end

is_finite(M::FPModule{<:FinFieldElem}) = true
is_known(::typeof(is_finite), ::FPModule{<:FinFieldElem}) = true

is_finitely_generated(M::FPModule) = true
is_finitely_generated(M::Module) = isfinite(ngens(M)) || throw(NotImplementedError(:is_finitely_generated, M))

@doc raw"""
    is_noetherian(M::Module)

Check if the module $M$ is Noetherian.

# Examples
```jldoctest
julia> R, x = polynomial_ring(ZZ, [:x]);

julia> M = free_module(R, 2)
Free module of rank 2 over R

julia> is_noetherian(M)
true
```
"""
function is_noetherian(M::Module)
  is_finitely_generated(M) || return false
  is_noetherian(base_ring(M)) && return true
  throw(NotImplementedError(:is_noetherian, M))
end

function is_sub_with_data(M::FPModule{T}, N::FPModule{T}) where T <: RingElement
  fl = is_submodule(N, M)
  if fl
    return fl, hom(M, N, elem_type(N)[N(m) for m = gens(M)])
  else
    return fl, hom(M, N, elem_type(N)[zero(N) for m = gens(M)])
  end
end

Base.issubset(M::FPModule{T}, N::FPModule{T}) where T <: RingElement = is_submodule(M, N)

order(M::FPModule{<:FinFieldElem}) = order(base_ring(M))^dim(M)

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(v::FPModuleElem{T}) where T <: RingElement
   N = parent(v)
   return N(-Generic._matrix(v))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(v1::FPModuleElem{T}, v2::FPModuleElem{T}) where T <: RingElement
   check_parent(v1, v2)
   N = parent(v1)
   return N(Generic._matrix(v1) + Generic._matrix(v2))
end

function -(v1::FPModuleElem{T}, v2::FPModuleElem{T}) where T <: RingElement
   check_parent(v1, v2)
   N = parent(v1)
   return N(Generic._matrix(v1) - Generic._matrix(v2))
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(v::FPModuleElem{T}, c::T) where T <: RingElem
   base_ring(v) != parent(c) && error("Incompatible rings")
   N = parent(v)
   return N(Generic._matrix(v) * c)
end

function *(v::FPModuleElem{T}, c::U) where {T <: RingElement, U <: Union{Rational, Integer}}
   N = parent(v)
   return N(Generic._matrix(v) * c)
end

function *(c::T, v::FPModuleElem{T}) where T <: RingElem
   base_ring(v) != parent(c) && error("Incompatible rings")
   N = parent(v)
   return N(c * Generic._matrix(v))
end

function *(c::U, v::FPModuleElem{T}) where {T <: RingElement, U <: Union{Rational, Integer}}
   N = parent(v)
   return N(c * Generic._matrix(v))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(m::FPModuleElem{T}, n::FPModuleElem{T}) where T <: RingElement
   check_parent(m, n)
   return Generic._matrix(m) == Generic._matrix(n)
end

function hash(m::FPModuleElem{T}, h::UInt)  where T <: RingElement
   b = 0xe08f5b4ea1cd9a12%UInt
   return xor(hash(Generic._matrix(m), h), b)
end

###############################################################################
#
#   Intersection
#
###############################################################################

@doc raw"""
    intersect(M::FPModule{T}, N::FPModule{T}) where T <: RingElement

Return the intersection of the modules $M$ as a submodule of $M$. Note that
$M$ and $N$ must be (constructed as) submodules (transitively) of some common
module $P$.
"""
function intersect(M::FPModule{T}, N::FPModule{T}) where T <: RingElement
   check_parent(M, N)
   # Compute the common supermodule P of M and N
   flag, P = is_compatible(M, N)
   !flag && error("Modules not compatible")
   # Compute the generators of M as elements of P
   G1 = gens(M)
   M1 = M
   while M1 !== P
      _map1 = M1.map
      G1 = [_map1(v) for v in G1]
      M1 = supermodule(M1)
   end
   # Compute the generators of N as elements of P
   G2 = gens(N)
   M2 = N
   while M2 !== P
      _map2 = M2.map
      G2 = [_map2(v) for v in G2]
      M2 = supermodule(M2)
   end
   # Make matrix containing all generators and relations as rows
   r1 = ngens(M)
   r2 = ngens(N)
   prels = rels(P)
   r3 = length(prels)
   c = ngens(P)
   mat = zero_matrix(base_ring(M), r1 + r2 + r3, c)
   # We flip the rows of the matrix so the input to Submodule is in upper
   # triangular form
   rn = r1 + r2 + r3
   for i = 1:r1
      for j = 1:c
         mat[rn - i + 1, j] = Generic._matrix(G1[i])[1, j]
      end
   end
   for i = 1:r2
      for j = 1:c
         mat[rn - i - r1 + 1, j] = Generic._matrix(G2[i])[1, j]
      end
   end
   for i = 1:r3
      for j = 1:c
         mat[rn - i - r1 - r2 + 1, j] = prels[i][1, j]
      end
   end
   # Find the left kernel space of the matrix
   K = kernel(mat)
   nc = nrows(K)
   # Last r1 elements of a row correspond to a generators of intersection
   # We flip the rows of K so the input to Submodule is upper triangular
   # and the columns so that they correspond to the original order before
   # flipping above
   I = [M(T[K[nc - j + 1, rn - i + 1] for i in 1:r1]) for j in 1:nc]
   return sub(M, I)
end

###############################################################################
#
#   Comparison
#
###############################################################################

@doc raw"""
    ==(M::FPModule{T}, N::FPModule{T}) where T <: RingElement

Return `true` if the modules are (constructed to be) the same module
elementwise. This is not object equality and it is not isomorphism. In fact,
each method of constructing modules (submodules, quotient modules, products,
etc.) must extend this notion of equality to the modules they create.
"""
function ==(M::FPModule{T}, N::FPModule{T}) where T <: RingElement
   M === N && return true  #object equality is sufficient
   check_parent(M, N)
   # Compute the common supermodule P of M and N
   flag, P = is_compatible(M, N)
   !flag && error("Modules not compatible")
   # Compute the generators of M as elements of P
   G1 = gens(M)
   M1 = M
   while M1 !== P
      _map1 = M1.map
      G1 = [_map1(v) for v in G1]
      M1 = supermodule(M1)
   end
   # Compute the generators of N as elements of P
   G2 = gens(N)
   M2 = N
   while M2 !== P
      _map2 = M2.map
      G2 = [_map2(v) for v in G2]
      M2 = supermodule(M2)
   end
   # Put (rewritten) gens of M and N into matrices with relations of P
   prels = rels(P)
   c = ngens(P)
   r1 = ngens(M)
   r2 = ngens(N)
   mat1 = zero_matrix(base_ring(M), r1 + length(prels), c)
   for i = 1:r1
      for j = 1:c
         mat1[i, j] = Generic._matrix(G1[i])[1, j]
      end
   end
   mat2 = zero_matrix(base_ring(M), r2 + length(prels), c)
   for i = 1:r2
      for j = 1:c
         mat2[i, j] = Generic._matrix(G2[i])[1, j]
      end
   end
   for i = 1:length(prels)
      for j = 1:c
         mat1[i + r1, j] = prels[i][1, j]
         mat2[i + r2, j] = prels[i][1, j]
      end
   end
   sol_ctx1 = solve_init(mat1)
   sol_ctx2 = solve_init(mat2)
   # Check containment of rewritten gens of M in row space of mat2
   for v in G1
      if !can_solve(sol_ctx2, Generic._matrix(v))
         return false
      end
   end
   # Check containment of rewritten gens of N in row space of mat1
   for v in G2
      if !can_solve(sol_ctx1, Generic._matrix(v))
         return false
      end
   end
   return true
end

# Equal `FPModule`s must have the same `hash` value.
Base.hash(M::FPModule, h::UInt) = h # FIXME

###############################################################################
#
#   Isomorphism
#
###############################################################################

@doc raw"""
    is_isomorphic(M::FPModule{T}, N::FPModule{T}) where T <: RingElement

Return `true` if the modules $M$ and $N$ are isomorphic.
"""
function is_isomorphic(M::FPModule{T}, N::FPModule{T}) where T <: RingElement
   return invariant_factors(M) == invariant_factors(N)
end

###############################################################################
#
#   Module element access
#
###############################################################################

@doc raw"""
    getindex(v::FPModuleElem{T}, i::Int) where T <: RingElement

Return the $i$-th coefficient of the module element $v$.
"""
function getindex(v::FPModuleElem{T}, i::Int) where T <: RingElement
   return Generic._matrix(v)[1, i]
end

@doc raw"""
    coordinates(v::FPModuleElem{T}, i::Int) where T <: RingElement

Return the coordinates of the module element $v$ as a `Vector{T}`.
"""
function coordinates(v::FPModuleElem{T}) where T <: RingElement
   return collect(Generic._matrix(v)[1, :])
end

###############################################################################
#
#   Random generation
#
###############################################################################

RandomExtensions.maketype(M::FPModule, _) = elem_type(M)

function RandomExtensions.make(M::FPModule, vs...)
   R = base_ring(M)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(M, vs[1]) # forward to default Make constructor
   else
      Make(M, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{
                 <:FPModuleElem, <:FPModule}})
   M, vals = sp[][1:end]
   M(rand(rng, vals, ngens(M)))
end

function rand(rng::AbstractRNG, M::FPModule{T}, vals...) where T <: RingElement
   rand(rng, make(M, vals...))
end

rand(M::FPModule, vals...) = rand(Random.default_rng(), M, vals...)

###############################################################################
#
#   Iteration
#
###############################################################################

Base.length(M::FPModule{T}) where T <: FinFieldElem = Int(order(M))

function Base.iterate(M::FPModule{T}) where T <: FinFieldElem
  k = base_ring(M)
  if dim(M) == 0
    return zero(M), iterate([1])
  end
  p = Base.Iterators.ProductIterator(Tuple([k for i=1:dim(M)]))
  f = iterate(p)
  @assert f !== nothing
  return M(elem_type(k)[f[1][i] for i=1:dim(M)]), (f[2], p)
end

function Base.iterate(M::FPModule{T}, st::Tuple{<:Tuple, <:Base.Iterators.ProductIterator}) where T <: FinFieldElem
  n = iterate(st[2], st[1])
  if n === nothing
    return n
  end
  return M(elem_type(base_ring(M))[n[1][i] for i=1:dim(M)]), (n[2], st[2])
end

function Base.iterate(::FPModule{<:FinFieldElem}, ::Tuple{Int64, Int64})
  return nothing
end
