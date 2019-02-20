###############################################################################
#
#   ModuleMorphism.jl : Morphisms of sub/quotient modules
#
###############################################################################

export ModuleMorphism

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, f::ModuleMorphism)
   println(io, "Module morphism with")
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

function (f::ModuleMorphism{T})(a::AbstractAlgebra.ModuleElem{T}) where T <: RingElement
   parent(a) != domain(f) && error("Incompatible module element")
   return image_fn(f)(a)
end

###############################################################################
#
#   ModuleMorphism constructor
#
###############################################################################

@doc Markdown.doc"""
    ModuleMorphism(M1::AbstractAlgebra.Module{T}, M2::AbstractAlgebra.Module{T}, m::AbstractAlgebra.MatElem{T}) where T <: RingElement
> Create the morphism $f : M_1 \to M_2$ represented by the matrix $m$.
"""
function ModuleMorphism(M1::AbstractAlgebra.Module{T}, M2::AbstractAlgebra.Module{T}, m::AbstractAlgebra.MatElem{T}) where T <: RingElement
   return ModuleMorphism{T}(M1, M2, m)
end

