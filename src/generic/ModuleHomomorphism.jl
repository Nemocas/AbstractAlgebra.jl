###############################################################################
#
#   ModuleHomomorphism.jl : Generic homomorphisms of free/sub/quotient modules
#
###############################################################################

###############################################################################
#
#   Basic manipulation
#
###############################################################################

inverse_mat(f::Map(ModuleIsomorphism)) = f.inverse_matrix

inverse_image_fn(f::Map(ModuleIsomorphism)) = f.inverse_image_fn

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, f::Map(ModuleIsomorphism))
   println(io, "Module isomorphism with")
   print(io, "Domain: ")
   print(IOContext(io, :compact => true), domain(f))
   println(io, "")
   print(io, "Codomain: ")
   print(IOContext(io, :compact => true), codomain(f))
end

###############################################################################
#
#   Composition
#
###############################################################################

function compose(f::Map(ModuleIsomorphism), g::Map(ModuleIsomorphism))
   check_composable(f, g)
   T = elem_type(base_ring(domain(f)))
   return ModuleIsomorphism{T}(domain(f), codomain(g), f.matrix*g.matrix,
                                             g.inverse_matrix*f.inverse_matrix)
end

###############################################################################
#
#   Inverse
#
###############################################################################

@doc Markdown.doc"""
    Base.inv(f::Map(ModuleIsomorphism))

Return the inverse map of the given module isomorphism. This is computed
cheaply.
"""
function Base.inv(f::Map(ModuleIsomorphism))
   T = elem_type(base_ring(domain(f)))
   return ModuleIsomorphism{T}(codomain(f), domain(f), inverse_mat(f), mat(f))
end

###############################################################################
#
#   Call overload
#
###############################################################################

function (f::ModuleHomomorphism{T})(a::AbstractAlgebra.FPModuleElem{T}) where T <: RingElement
   parent(a) !== domain(f) && error("Incompatible module element")
   return image_fn(f)(a)
end

function (f::ModuleIsomorphism{T})(a::AbstractAlgebra.FPModuleElem{T}) where T <: RingElement
   parent(a) !== domain(f) && error("Incompatible module element")
   return image_fn(f)(a)
end

###############################################################################
#
#   ModuleHomomorphism constructor
#
###############################################################################

function ModuleHomomorphism(M1::AbstractAlgebra.FPModule{T},
     M2::AbstractAlgebra.FPModule{T}, m::AbstractAlgebra.MatElem{T}) where
                                                               T <: RingElement
   (nrows(m) == ngens(M1) && ncols(m) == ngens(M2)) ||
                                                    error("dimension mismatch")
   return ModuleHomomorphism{T}(M1, M2, m)
end

function ModuleHomomorphism(M1::AbstractAlgebra.FPModule{T},
                 M2::AbstractAlgebra.FPModule{T}, v::Array{S, 1}) where
                         {T <: RingElement, S<:AbstractAlgebra.FPModuleElem{T}}
   return ModuleHomomorphism(M1, M2, vcat([Generic._matrix(x) for x = v]...))
end

function ModuleIsomorphism(M1::AbstractAlgebra.FPModule{T},
     M2::AbstractAlgebra.FPModule{T}, M::AbstractAlgebra.MatElem{T}) where
                                                               T <: RingElement
   # Put rows of m and target relations into a matrix
   R = base_ring(M1)
   R != base_ring(M2) && error("Incompatible modules")
   trels = rels(M2)
   q = length(trels)
   m = nrows(M)
   n = ncols(M)
   (ngens(M1) != m || ngens(M2) != n) &&
                                        error("Matrix of the wrong dimensions")
   if m == 0 || n == 0
       M_inv = matrix(R, n, m, T[])
   else
      # Put matrix M and target relations in a matrix
      mat = zero_matrix(R, m + q, n)
      for i = 1:m
         for j = 1:n
            mat[i, j] = M[i, j]
         end
      end
      for i = 1:q
         for j = 1:n
            mat[m + i, j] = trels[i].v[1, j]
         end
      end
      # Find left inverse of mat
      N = identity_matrix(R, n)
      X = solve_left(mat, N)
      # Construct matrix of inverse homomorphism from first m columns of X
      M_inv = X[:, 1:m]
   end
   return ModuleIsomorphism{T}(M1, M2, M, M_inv)
end
