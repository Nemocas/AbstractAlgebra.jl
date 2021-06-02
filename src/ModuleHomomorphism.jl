###############################################################################
#
#   ModuleHomomorphism.jl : Homomorphisms of free/sub/quotient modules
#
###############################################################################

export ModuleHomomorphism, image, inverse_image_fn, inverse_mat, mat,
       preimage

###############################################################################
#
#   Basic manipulation
#
###############################################################################

mat(f::Map(FPModuleHomomorphism)) = f.matrix

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, f::Map(FPModuleHomomorphism))
   println(io, "Module homomorphism with")
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

function compose(f::Map(FPModuleHomomorphism), g::Map(FPModuleHomomorphism))
   check_composable(f, g)
   return ModuleHomomorphism(domain(f), codomain(g), f.matrix*g.matrix)
end

###############################################################################
#
#   Kernel
#
###############################################################################

@doc Markdown.doc"""
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
   M = mat(f)
   # put domain relations and M in a big matrix
   # swap rows so we can get upper triangular wrt original data
   nr = nrows(M) + length(crels)
   N = zero_matrix(R, nr, ncols(M))
   for i = 1:nrows(M)
      for j = 1:ncols(M)
         N[nr - i + 1, j] = M[i, j]
      end
   end
   for i = 1:length(crels)
      for j = 1:ncols(M)
         N[nr - i - nrows(M) + 1, j] = crels[i][1, j]
      end
   end
   # compute the kernel
   num_gens, K = left_kernel(N)
   # Construct generators of kernel submodule, reversing rows
   # and columns so they're correct wrt to original data and
   # in upper triangular form
   V = Vector{elem_type(D)}(undef, num_gens)
   for j = 1:num_gens
      V[j] = D([K[num_gens - j + 1, nr - k + 1] for k = 1:nrows(M)])
   end
   return sub(D, V)
end

###############################################################################
#
#   Image
#
###############################################################################

@doc Markdown.doc"""
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

@doc Markdown.doc"""
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
   D = domain(f)
   C = codomain(f)
   R = base_ring(C)
   parent(v) !== C && error("Incompatible element")
   M = mat(f)
   trels = rels(C)
   # Put rows of M and target relations into a matrix
   q = length(trels)
   m = nrows(M)
   n = ncols(M)
   ncols(Generic._matrix(v)) != n && error("Incompatible element")
   if m == 0 || n == 0
       return D(zero_matrix(R, 1, m))
   else
      # Put matrix M and target relations in a matrix
      matr = zero_matrix(R, m + q, n)
      for i = 1:m
         for j = 1:n
            matr[i, j] = M[i, j]
         end
      end
      for i = 1:q
         for j = 1:n
            matr[m + i, j] = trels[i][1, j]
         end
      end
      # Find left inverse of mat
      x = solve_left(matr, Generic._matrix(v))
      if q != 0
         x = matrix(R, 1, m, T[x[1, i] for i in 1:m])
      end
      return D(x)
   end
end

###############################################################################
#
#   ModuleHomomorphism constructor
#
###############################################################################

@doc Markdown.doc"""
    ModuleHomomorphism(M1::FPModule{T},
	                      M2::FPModule{T}, m::MatElem{T}) where T <: RingElement
Create the homomorphism $f : M_1 \to M_2$ represented by the matrix $m$.
"""
function ModuleHomomorphism(M1::FPModule{T},
                         M2::FPModule{T}, m::MatElem{T}) where T <: RingElement
   return Generic.ModuleHomomorphism{T}(M1, M2, m)
end

function ModuleHomomorphism(M1::FPModule{T},
               M2::FPModule{T}, v::Array{S, 1}) where
                        {T <: RingElement, S<:FPModuleElem{T}}
   return Generic.ModuleHomomorphism(M1, M2, v)
end

@doc Markdown.doc"""
    ModuleIsomorphism(M1::FPModule{T}, M2::FPModule{T}, M::MatElem{T},
	                                     minv::MatElem{T}) where T <: RingElement
Create the isomorphism $f : M_1 \to M_2$ represented by the matrix $M$. The
inverse morphism is automatically computed.
"""
function ModuleIsomorphism(M1::FPModule{T},
                          M2::FPModule{T}, M::MatElem{T}) where T <: RingElement
   return Generic.ModuleIsomorphism(M1, M2, M)
end
