###############################################################################
#
#   ModuleHomomorphism.jl : Homomorphisms of free/sub/quotient modules
#
###############################################################################

export ModuleHomomorphism

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, f::ModuleHomomorphism)
   println(io, "Module homomorphism with")
   print(io, "Domain: ")
   show(io, domain(f))
   println(io, "")
   print(io, "Codomain: ")
   show(io, codomain(f))
end

###############################################################################
#
#   Call overload
#
###############################################################################

function (f::ModuleHomomorphism{T})(a::AbstractAlgebra.ModuleElem{T}) where T <: RingElement
   parent(a) != domain(f) && error("Incompatible module element")
   return image_fn(f)(a)
end

###############################################################################
#
#   ModuleHomomorphism constructor
#
###############################################################################

@doc Markdown.doc"""
    ModuleHomomorphism(M1::AbstractAlgebra.Module{T}, M2::AbstractAlgebra.Module{T}, m::AbstractAlgebra.MatElem{T}) where T <: RingElement
> Create the homomorphism $f : M_1 \to M_2$ represented by the matrix $m$.
"""
function ModuleHomomorphism(M1::AbstractAlgebra.Module{T}, M2::AbstractAlgebra.Module{T}, m::AbstractAlgebra.MatElem{T}) where T <: RingElement
   return ModuleHomomorphism{T}(M1, M2, m)
end

