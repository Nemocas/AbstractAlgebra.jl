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

domain(f::ModuleHomomorphism) = f.domain
codomain(f::ModuleHomomorphism) = f.codomain
image_fn(f::ModuleHomomorphism) = f.image_fn
matrix(f::ModuleHomomorphism) = f.matrix

domain(f::ModuleIsomorphism) = f.domain
codomain(f::ModuleIsomorphism) = f.codomain
image_fn(f::ModuleIsomorphism) = f.image_fn
inverse_image_fn(f::ModuleIsomorphism) = f.inverse_image_fn
matrix(f::ModuleIsomorphism) = f.matrix
inverse_mat(f::ModuleIsomorphism) = f.inverse_matrix

###############################################################################
#
#   Unary operators
#
###############################################################################

Base.:-(a::ModuleHomomorphism) = hom(domain(a), codomain(a), -matrix(a))

###############################################################################
#
#   Binary operators
#
###############################################################################

Base.:*(a::T, b::ModuleHomomorphism{T}) where {T <: RingElement} = hom(domain(b), codomain(b), a * matrix(b))
Base.:*(a::T, b::ModuleIsomorphism{T}) where {T <: RingElement} = hom(domain(b), codomain(b), a * matrix(b))
Base.:+(a::ModuleHomomorphism, b::ModuleHomomorphism) = hom(domain(a), codomain(a), matrix(a) + matrix(b))
Base.:-(a::ModuleHomomorphism, b::ModuleHomomorphism) = hom(domain(a), codomain(a), matrix(a) - matrix(b))

###############################################################################
#
#   Comparison
#
###############################################################################

function Base.:(==)(a::Union{ModuleHomomorphism, ModuleIsomorphism}, b::Union{ModuleHomomorphism, ModuleIsomorphism})
  domain(a) === domain(b) || return false
  codomain(a) === codomain(b) || return false
  return matrix(a) == matrix(b)
end

function Base.hash(a::Union{ModuleHomomorphism, ModuleIsomorphism}, h::UInt)
  h = hash(domain(a), h)
  h = hash(codomain(a), h)
  h = hash(matrix(a), h)
  return h
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, f::Map(ModuleIsomorphism))
  if is_terse(io)
    print(io, "Module isomorphism")
  else
    io = pretty(io)
    io = terse(io)
    print(io, "Hom: ", Lowercase(), domain(f))
    print(io, " -> ", Lowercase(), codomain(f))
  end
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
#   Preimage
#
###############################################################################

function AbstractAlgebra.solve_ctx(f::ModuleHomomorphism)
  # cache solver context
  if !isdefined(f, :solve_ctx)
    f.solve_ctx = AbstractAlgebra.compute_solve_ctx(f)
  end
  return f.solve_ctx::AbstractAlgebra.Solve.solve_context_type(base_ring(codomain(f)))
end

###############################################################################
#
#   Inverse
#
###############################################################################

@doc raw"""
    Base.inv(f::Map(ModuleIsomorphism))

Return the inverse map of the given module isomorphism. This is computed
cheaply.
"""
function Base.inv(f::Map(ModuleIsomorphism))
   T = elem_type(base_ring(domain(f)))
   return ModuleIsomorphism{T}(codomain(f), domain(f), inverse_mat(f), matrix(f))
end

function Base.inv(f::ModuleHomomorphism)
  return hom(codomain(f), domain(f), inv(matrix(f)))
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
                 M2::AbstractAlgebra.FPModule{T}, v::Vector{S}) where
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
            mat[m + i, j] = trels[i][1, j]
         end
      end
      # Find left inverse of mat
      N = identity_matrix(R, n)
      X = solve(mat, N)
      # Construct matrix of inverse homomorphism from first m columns of X
      M_inv = X[:, 1:m]
   end
   return ModuleIsomorphism{T}(M1, M2, M, M_inv)
end

function hom(V::AbstractAlgebra.Module, W::AbstractAlgebra.Module, v::Vector{<:ModuleElem}; check::Bool = true)
  if ngens(V) == 0
    return ModuleHomomorphism(V, W, zero_matrix(base_ring(V), ngens(V), ngens(W)))
  end
  return ModuleHomomorphism(V, W, reduce(vcat, [x.v for x = v]))
end

function hom(V::AbstractAlgebra.Module, W::AbstractAlgebra.Module, v::MatElem; check::Bool = true)
  return ModuleHomomorphism(V, W, v)
end
