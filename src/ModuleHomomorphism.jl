###############################################################################
#
#   ModuleHomomorphism.jl : Homomorphisms of free/sub/quotient modules
#
###############################################################################

###############################################################################
#
#   Basic manipulation
#
###############################################################################

matrix(f::Map(FPModuleHomomorphism)) = throw(NotImplementedError(:matrix, f))

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, f::Map(FPModuleHomomorphism))
  if is_terse(io)
    print(io, "Module homomorphism")
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

function compose(f::Map(FPModuleHomomorphism), g::Map(FPModuleHomomorphism))
   check_composable(f, g)
   return ModuleHomomorphism(domain(f), codomain(g), matrix(f)*matrix(g))
end

###############################################################################
#
#   Kernel
#
###############################################################################

@doc raw"""
    kernel(f::ModuleHomomorphism{T}) where T <: RingElement

Return a pair `K, g` consisting of the kernel object $K$ of the given module
homomorphism $f$ (as a submodule of its domain) and the canonical injection
from the kernel into the domain of $f$.
"""
function kernel(f::Map(FPModuleHomomorphism))
   D = domain(f)
   C = codomain(f)
   R = base_ring(D)
   crels = rels(C)
   M = matrix(f)
   # put domain relations and M in a big matrix
   nr = nrows(M) + length(crels)
   N = M
   if length(crels) != 0
      NN = reduce(vcat, crels)
      N = vcat(N, NN)
   end
   # compute the kernel
   K = kernel(N)
   V = [D(K[j:j, 1:nrows(M)]) for j in 1:nrows(K)]
   return sub(D, V)
end

###############################################################################
#
#   Image
#
###############################################################################

@doc raw"""
    image(f::Map(FPModuleHomomorphism))

Return a pair `I, g` consisting of the image object $I$ of the given module
homomorphism $f$ (as a submodule of its codomain) and the canonical injection
from the image into the codomain of $f$
"""
function image(f::Map(FPModuleHomomorphism))
   D = domain(f)
   C = codomain(f)
   R = base_ring(D)
   G = gens(D)
   V = elem_type(C)[f(v) for v in G]
   return sub(C, V)
end

###############################################################################
#
#   Preimage
#
###############################################################################

@doc raw"""
    preimage(f::Map(FPModuleHomomorphism),
             v::FPModuleElem{T}) where T <: RingElement
Return a preimage of $v$ under the homomorphism $f$, i.e. an element of the
domain of $f$ that maps to $v$ under $f$. Note that this has no special
mathematical properties. It is an element of the set theoretical preimage of
the map $f$ as a map of sets, if one exists. The preimage is neither
unique nor chosen in a canonical way in general. When no such element exists,
an exception is raised.
"""
function preimage(f::Map(FPModuleHomomorphism), v::FPModuleElem{T}) where {T<:RingElement}
  return preimage(f, [v])[1]
end

function preimage(
  f::Map(FPModuleHomomorphism), v::Vector{<:FPModuleElem{T}}
) where {T<:RingElement}
  fl, b = has_preimage_with_preimage(f, v)
  !fl && error("Element has no preimage")
  return b
end

@doc raw"""
    has_preimage_with_preimage(f::Map(FPModuleHomomorphism),
      v::FPModuleElem{T}) where T <: RingElement

Check if $v$ has a preimage under the homomorphism $f$.
If it does, return a tuple (true, y) for $y$ in domain(f) such that $f(y) = x$ holds,
otherwise, return (false, id) where id is the identity of domain(f).
"""
function has_preimage_with_preimage(
  f::Map(FPModuleHomomorphism), v::FPModuleElem{T}
) where {T<:RingElement}
  fl, b = has_preimage_with_preimage(f, [v])
  return fl, b[1]
end

function has_preimage_with_preimage(
  f::Map(FPModuleHomomorphism), v::Vector{<:FPModuleElem{T}}
) where {T<:RingElement}
  D = domain(f)
  C = codomain(f)
  R = base_ring(C)
  if length(v) == 0
    return true, elem_type(domain(f))[]
  end
  parent(v[1]) !== C && error("Incompatible element")
  M = matrix(f)
  trels = rels(C)
  # Put rows of M and target relations into a matrix
  q = length(trels)
  m = nrows(M)
  n = ncols(M)
  if n == 0 #target module is trivial, so v is all zero
    return true, elem_type(D)[D(zero_matrix(R, 1, m)) for x in v]
  elseif m == 0 #source module is trivial
    return all(iszero, v), elem_type(D)[D(zero_matrix(R, 1, m)) for x in v]
  else
    if !isdefined(f, :solve_ctx)
      # Put matrix M and target relations in a matrix
      matr = zero_matrix(R, m + q, n)
      matr[1:m, 1:n] = M
      for i in 1:q
        matr[m + i, :] = trels[i]
      end
      # Find left inverse of mat
      f.solve_ctx = solve_init(matr)
    end
    inmat = reduce(vcat, Generic._matrix.(v))
    fl, x = can_solve_with_solution(f.solve_ctx, inmat)
    if fl
      return true, elem_type(D)[D(x[i:i, 1:m]) for i in 1:length(v)]
    else
      return false, elem_type(D)[D(zero_matrix(R, 1, m)) for x in v]
    end
  end
end

###############################################################################
#
#   ModuleHomomorphism constructor
#
###############################################################################

@doc raw"""
    ModuleHomomorphism(M1::FPModule{T},
                       M2::FPModule{T}, m::MatElem{T}) where T <: RingElement
Create the homomorphism $f : M_1 \to M_2$ represented by the matrix $m$.
"""
function ModuleHomomorphism(M1::FPModule{T},
                         M2::FPModule{T}, m::MatElem{T}) where T <: RingElement
   return Generic.ModuleHomomorphism{T}(M1, M2, m)
end

function ModuleHomomorphism(M1::FPModule{T},
               M2::FPModule{T}, v::Vector{S}) where
                        {T <: RingElement, S<:FPModuleElem{T}}
   return Generic.ModuleHomomorphism(M1, M2, v)
end

function ModuleHomomorphism(M1::Module, M2::Module, A...)
   Generic.ModuleHomomorphism(M1, M2, A...)
end

function module_homomorphism(M1::Module, M2::Module, m::MatElem)
   Generic.ModuleHomomorphism(M1, M2, m)
end

@doc raw"""
    ModuleIsomorphism(M1::FPModule{T}, M2::FPModule{T}, M::MatElem{T},
                      minv::MatElem{T}) where T <: RingElement
Create the isomorphism $f : M_1 \to M_2$ represented by the matrix $M$. The
inverse morphism is automatically computed.
"""
function ModuleIsomorphism(M1::FPModule{T},
                          M2::FPModule{T}, M::MatElem{T}) where T <: RingElement
   return Generic.ModuleIsomorphism(M1, M2, M)
end

function ModuleIsomorphism(M1::Module, M2::Module, m::MatElem)
   Generic.ModuleIsomorphism(M1, M2, m)
end

function module_isomorphism(M1::Module, M2::Module, m::MatElem)
   Generic.ModuleIsomorphism(M1, M2, m)
end
