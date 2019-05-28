###############################################################################
#
#   ModuleHomomorphism.jl : Homomorphisms of free/sub/quotient modules
#
###############################################################################

export ModuleHomomorphism, image

###############################################################################
#
#   Basic manipulation
#
###############################################################################

mat(f::ModuleHomomorphism{T}) where T <: RingElement = f.matrix

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, f::ModuleHomomorphism)
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

function compose(f::ModuleHomomorphism{T}, g::ModuleHomomorphism{T}) where T <: RingElement
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
> Returns a pair `K, g` consisting of the kernel object $K$ of the given module
> homomorphism $f$ (as a submodule of its domain) and the canonical injection
> from the kernel into the domain of $f$
"""
function kernel(f::ModuleHomomorphism{T}) where T <: RingElement
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
   return Submodule(D, V)
end

###############################################################################
#
#   Image
#
###############################################################################

@doc Markdown.doc"""
    image(f::ModuleHomomorphism{T}) where T <: RingElement
> Returns a pair `I, g` consisting of the image object $I$ of the given module
> homomorphism $f$ (as a submodule of its codomain) and the canonical injection
> from the image into the codomain of $f$
"""
function image(f::ModuleHomomorphism{T}) where T <: RingElement
   D = domain(f)
   C = codomain(f)
   R = base_ring(D)
   G = gens(D)
   V = elem_type(C)[f(v) for v in G]
   return Submodule(C, V)
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

###############################################################################
#
#   ModuleHomomorphism constructor
#
###############################################################################

@doc Markdown.doc"""
    ModuleHomomorphism(M1::AbstractAlgebra.FPModule{T}, M2::AbstractAlgebra.FPModule{T}, m::AbstractAlgebra.MatElem{T}) where T <: RingElement
> Create the homomorphism $f : M_1 \to M_2$ represented by the matrix $m$.
"""
function ModuleHomomorphism(M1::AbstractAlgebra.FPModule{T}, M2::AbstractAlgebra.FPModule{T}, m::AbstractAlgebra.MatElem{T}) where T <: RingElement
   return ModuleHomomorphism{T}(M1, M2, m)
end

