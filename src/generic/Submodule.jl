export Submodule, submodule_elem, ngens

base_ring(N::Submodule{T}) where T <: RingElement = N.base_ring

generators(N::Submodule{T}) where T <: RingElement = N.gens::Vector{elem_type(N.m)}

function Submodule(m::AbstractAlgebra.FPModule{T}, gens::Vector{S}) where {T <: RingElement, S <: AbstractAlgebra.FPModuleElem{T}}
   M = Submodule{T}(m, gens)
   return M, x-> x
end

function Submodule(m::AbstractAlgebra.FPModule{T}, gens::Vector{Any}) where T <: RingElement
   return Submodule(m, elem_type(m)[])
end

function Submodule(m::AbstractAlgebra.FPModule{T}, subs::Vector{Submodule{T}}) where T <: RingElement
   gens = vcat((generators(s) for s in subs)...)
   return Submodule(m, gens)
end

