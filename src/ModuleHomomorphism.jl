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

matrix(f::Map(FPModuleHomomorphism)) = f.matrix

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, f::Map(FPModuleHomomorphism))
  if get(io, :supercompact, false)
    print(io, "Module homomorphism")
  else
    io = pretty(io)
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
   return ModuleHomomorphism(domain(f), codomain(g), f.matrix*g.matrix)
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
from the kernel into the domain of $f$
"""
function kernel(f::Map(FPModuleHomomorphism))
   D = domain(f)
   C = codomain(f)
   R = base_ring(D)
   crels = rels(C)
   M = matrix(f)
   # put domain relations and M in a big matrix
   # swap rows so we can get upper triangular wrt original data
   nr = nrows(M) + length(crels)
   local N
   if length(crels) == 0
     N = M
   else
     N = zero_matrix(R, nr, ncols(M))
     for i = 1:nrows(M)
        N[nr - i + 1, :] = view(M, i:i, :)
     end
     for i = 1:length(crels)
        N[nr - i - nrows(M) + 1, :] = crels[i]
     end
   end
   # compute the kernel
   num_gens, K = AbstractAlgebra._left_kernel(N)
   # Construct generators of kernel submodule, reversing rows
   # and columns so they're correct wrt to original data and
   # in upper triangular form
   #TODO use other primitives, don't ever use getindex, ...
   V = Vector{elem_type(D)}(undef, num_gens)
   if length(crels) == 0
     for j=1:num_gens
       V[j] = D(view(K, j:j, :))
     end
   else
     for j = 1:num_gens
        V[j] = D([K[num_gens - j + 1, nr - k + 1] for k = 1:nrows(M)])
     end
   end
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
function preimage(f::Map(FPModuleHomomorphism), v::FPModuleElem{T}) where
                                                               T <: RingElement
   return preimage(f, [v])[1]
end

function preimage(f::Map(FPModuleHomomorphism), v::Vector{<:FPModuleElem{T}}) where
                                                               T <: RingElement
   D = domain(f)
   C = codomain(f)
   R = base_ring(C)
   if length(v) == 0
     return elem_type(domain(f))[]
   end
   parent(v[1]) !== C && error("Incompatible element")
   M = matrix(f)
   trels = rels(C)
   # Put rows of M and target relations into a matrix
   q = length(trels)
   m = nrows(M)
   n = ncols(M)
   if m == 0 || n == 0
     return [D(zero_matrix(R, 1, m)) for x = v]
   else
      # Put matrix M and target relations in a matrix
      matr = zero_matrix(R, m + q, n)
      matr[1:m, 1:n] = M
      for i = 1:q
        matr[m + i, :] = trels[i]
      end
      # Find left inverse of mat
      inmat = zero_matrix(R, length(v), n)
      for i=1:length(v)
        inmat[i, :] = Generic._matrix(v[i])
      end
      x = solve(matr, inmat)
      return [D(view(x, i:i,  1:m)) for i=1:length(v)]
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
