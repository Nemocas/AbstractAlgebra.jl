###############################################################################
#
#   SNFModule.jl : Invariant factor decomposition of modules
#
###############################################################################

export SNFModule, snf_module_elem,
       invariant_factors

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{snf_module_elem{T}}) where T <: RingElement = SNFModule{T}

elem_type(::Type{SNFModule{T}}) where T <: RingElement = snf_module_elem{T}

parent(v::snf_module_elem) = v.parent

base_ring(N::SNFModule{T}) where T <: RingElement = N.base_ring

base_ring(v::snf_module_elem{T}) where T <: RingElement = base_ring(v.parent)

ngens(N::SNFModule{T}) where T <: RingElement = length(N.invariant_factors)

gens(N::SNFModule{T}) where T <: RingElement = [gen(N, i) for i = 1:ngens(N)]

function gen(N::SNFModule{T}, i::Int) where T <: RingElement
   R = base_ring(N)
   return N([(j == i ? one(R) : zero(R)) for j = 1:ngens(N)])
end

invariant_factors(N::SNFModule{T}) where T <: RingElement = N.invariant_factors

function rels(N::SNFModule{T}) where T <: RingElement
   T1 = dense_matrix_type(T)
   R = base_ring(N)
   invs = invariant_factors(N)
   # count nonzero invariant factors
   num = length(invs)
   while num > 0
      if !iszero(invs[num])
         break
      end
      num -= 1
   end
   n = length(invs)
   r = T1[matrix(R, 1, n, T[i == j ? invs[i] : zero(R) for j in 1:n]) for i in 1:num]
   return r 
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, N::SNFModule{T}) where T <: RingElement
   print(io, "Invariant factor decomposed module over ")
   print(IOContext(io, :compact => true), base_ring(N))
   print(io, " with invariant factors ")
   print(IOContext(io, :compact => true), invariant_factors(N))
end

function show(io::IO, N::SNFModule{T}) where T <: FieldElement
   print(io, "Vector space over ")
   print(IOContext(io, :compact => true), base_ring(N))
   print(io, " with dimension ")
   print(io, ngens(N))
end

function show(io::IO, v::snf_module_elem)
   print(io, "(")
   len = ngens(parent(v))
   for i = 1:len - 1
      print(IOContext(io, :compact => true), v.v[1, i])
      print(io, ", ")
   end
   if len > 0
      print(IOContext(io, :compact => true), v.v[1, len])
   end
   print(io, ")")
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(v::snf_module_elem{T}) where T <: RingElement
   N = parent(v)
   return N(-v.v)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(v1::snf_module_elem{T}, v2::snf_module_elem{T}) where T <: RingElement
   check_parent(v1, v2)
   N = parent(v1)
   return N(v1.v + v2.v)
end

function -(v1::snf_module_elem{T}, v2::snf_module_elem{T}) where T <: RingElement
   check_parent(v1, v2)
   N = parent(v1)
   return N(v1.v - v2.v)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(v::snf_module_elem{T}, c::T) where T <: RingElem
   base_ring(v) != parent(c) && error("Incompatible rings")
   N = parent(v)
   return N(v.v*c)
end

function *(v::snf_module_elem{T}, c::U) where {T <: RingElement, U <: Union{Rational, Integer}}
   N = parent(v)
   return N(v.v*c)
end

function *(c::T, v::snf_module_elem{T}) where T <: RingElem
   base_ring(v) != parent(c) && error("Incompatible rings")
   N = parent(v)
   return N(c*v.v)
end

function *(c::U, v::snf_module_elem{T}) where {T <: RingElement, U <: Union{Rational, Integer}}
   N = parent(v)
   return N(c*v.v)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(m::snf_module_elem{T}, n::snf_module_elem{T}) where T <: RingElement
   check_parent(m, n)
   return m.v == n.v
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function reduce_mod_invariants(v::AbstractAlgebra.MatElem{T}, invars::Vector{T}) where T <: RingElement
   v = deepcopy(v) # don't modify input
   for i = 1:length(invars)
      if !iszero(invars[i])
         q, v[1, i] = AbstractAlgebra.divrem(v[1, i], invars[i])
      end
   end
   return v
end

function (N::SNFModule{T})(v::Vector{T}) where T <: RingElement
   length(v) != ngens(N) && error("Length of vector does not match number of generators")
   mat = matrix(base_ring(N), 1, length(v), v)
   mat = reduce_mod_invariants(mat, invariant_factors(N))
   return snf_module_elem{T}(N, mat)
end

function (M::SNFModule{T})(a::Vector{Any}) where T <: RingElement
   length(a) != 0 && error("Incompatible element")
   return M(T[])
end

function (N::SNFModule{T})(v::AbstractAlgebra.MatElem{T}) where T <: RingElement
   ncols(v) != ngens(N) && error("Length of vector does not match number of generators")
   nrows(v) != 1 && ("Not a vector in snf_module_elem constructor")
   v = reduce_mod_invariants(v, invariant_factors(N))
   return snf_module_elem{T}(N, v)
end

###############################################################################
#
#   SNFModule constructor
#
###############################################################################

@doc Markdown.doc"""
    snf(m::AbstractAlgebra.FPModule{T}) where T <: RingElement
> Return a pair `M, f` consisting of the invariant factor decomposition $M$ of
> the module `m` and a module homomorphism (isomorphisms) $f : M \to m$. The
> module `M` is itself a module which can be manipulated as any other module
> in the system.
"""
function snf(m::AbstractAlgebra.FPModule{T}) where T <: RingElement
   R = base_ring(m)
   old_rels = rels(m)
   # put the relations into a matrix
   r = length(old_rels)
   s = ngens(m)
   A = matrix(R, r, s, T[old_rels[i][1, j] for i in 1:r for j in 1:s])
   # compute the snf
   S, U, K = snf_with_transform(A)
   # compute K^-1
   K = inv(K)
   # count unit invariant factors
   nunits = 0
   while nunits < min(nrows(S), ncols(S))
      nunits += 1
      if !isunit(S[nunits, nunits])
         nunits -= 1
         break
      end
   end
   num_gens = nrows(S) - nunits
   # Make generators out of cols of matrix K
   # throwing away ones corresponding to unit invariant factors
   T2 = elem_type(m)
   gens = Vector{T2}(undef, ncols(A) - nunits)
   for i = 1:ncols(A) - nunits
      gens[i] = m(matrix(R, 1, ncols(K),
                    T[K[i + nunits, j]
                       for j in 1:ncols(K)]))
   end
   # extract invariant factors from S
   invariant_factors = T[S[i + nunits, i + nunits] for i in 1:num_gens]
   for i = num_gens + 1:ncols(A) - nunits
      push!(invariant_factors, zero(R))
   end
   # make matrix from gens
   mat = matrix(R, ncols(A) - nunits, ncols(K),
        T[gens[i].v[1, j] for i in 1:ncols(A) - nunits for j in 1:ncols(K)])
   M = SNFModule{T}(m, gens, invariant_factors)
   f = ModuleHomomorphism(M, m, mat)
   M.map = f 
   return M, f
end

function snf(m::SNFModule{T}) where T <: RingElement
   return m
end

@doc Markdown.doc"""
    invariant_factors(m::AbstractAlgebra.FPModule{T}) where T <: RingElement
> Return a vector of the invariant factors of the module $M$.
"""
function invariant_factors(m::AbstractAlgebra.FPModule{T}) where T <: RingElement
   R = base_ring(m)
   old_rels = rels(m)
   # put the relations into a matrix
   r = length(old_rels)
   s = ngens(m)
   A = matrix(R, r, s, T[old_rels[i][1, j] for i in 1:r for j in 1:s])
   # compute the snf
   S = snf(A)
   # count unit invariant factors
   nunits = 0
   while nunits < min(nrows(S), ncols(S))
      nunits += 1
      if !isunit(S[nunits, nunits])
         nunits -= 1
         break
      end
   end
   num_gens = nrows(S) - nunits
   # extract invariant factors from S
   invariant_factors = T[S[i + nunits, i + nunits] for i in 1:num_gens]
   for i = num_gens + 1:ncols(A) - nunits
      push!(invariant_factors, zero(R))
   end
   return invariant_factors
end

