mutable struct FreeModule{T <: Union{RingElement, NCRingElem}} <: AbstractAlgebra.FPModule{T}
   base_ring::NCRing
   rank::Int
end

mutable struct free_module_elem{T <: Union{RingElement, NCRingElem}} <: AbstractAlgebra.FPModuleElem{T}
    parent::FreeModule{T}
    v::Vector{T}
end

mutable struct Submodule{T <: RingElement} <: AbstractAlgebra.FPModule{T}
   m::AbstractAlgebra.FPModule{T}
   gens::Vector{<:AbstractAlgebra.FPModuleElem{T}}

   function Submodule{T}(M::AbstractAlgebra.FPModule{T}, gens::Vector{<:AbstractAlgebra.FPModuleElem{T}}) where T<: RingElement
      z = new{T}(M, gens)
   end
end

