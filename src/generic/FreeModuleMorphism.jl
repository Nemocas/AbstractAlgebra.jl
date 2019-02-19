###############################################################################
#
#   FreeModuleMorphism.jl : Morphisms of free modules
#
###############################################################################

export FreeModuleMorphism

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, f::FreeModuleMorphism)
   println(io, "Free module morphism with")
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

function (f::FreeModuleMorphism{T})(a::free_module_elem{T}) where T <: Union{RingElement, NCRingElem}
   parent(a) != domain(f) && error("Incompatible module element")
   return image_fn(f)(a)
end

###############################################################################
#
#   FreeModuleMorphism constructor
#
###############################################################################

@doc Markdown.doc"""
    FreeModuleMorphism(M1::FreeModule{T}, M2::FreeModule{T}, m::AbstractAlgebra.MatElem{T}) where T <: Union{RingElement, NCRingElem}
> Create the morphism $f : M_1 \to M_2$ represented by the matrix $m$.
"""
function FreeModuleMorphism(M1::FreeModule{T}, M2::FreeModule{T}, m::AbstractAlgebra.MatElem{T}) where T <: Union{RingElement, NCRingElem}
   return FreeModuleMorphism{T}(M1, M2, m)
end

