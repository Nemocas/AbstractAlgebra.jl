module AbstractAlgebra

export ZZ, FreeModule, Submodule, Generic, test_submodule

abstract type RingElem end
abstract type Ring end
abstract type FPModule{T} end
abstract type FPModuleElem{T} end

mutable struct Integers{T <: Integer} <: Ring
end

const RingElement = Union{RingElem, Integer, Rational, AbstractFloat}

#### Generic submodule ####

module Generic

import ..AbstractAlgebra: Ring, RingElem, RingElement

using ..AbstractAlgebra

export Submodule, FreeModule, free_module_elem

mutable struct FreeModule{T <: RingElement} <: AbstractAlgebra.FPModule{T}
   base_ring::Ring
   rank::Int
end

mutable struct free_module_elem{T <: RingElement} <: AbstractAlgebra.FPModuleElem{T}
    parent::FreeModule{T}
    v::Vector{T}
end

mutable struct Submodule{T <: RingElement} <: AbstractAlgebra.FPModule{T}
   m::AbstractAlgebra.FPModule{T}
   gens::Vector{<:AbstractAlgebra.FPModuleElem{T}}
end

elem_type(::FreeModule{T}) where T <: RingElement = free_module_elem{T}

function (M::FreeModule{T})(a::Vector{T}) where T <: RingElement
   return free_module_elem{T}(M, a)
end

function FreeModule(R::Ring, rank::Int)
   T = elem_type(R)
   return FreeModule{T}(R, rank)
end

generators(N::Submodule{T}) where T <: RingElement = N.gens::Vector{elem_type(N.m)}

function Submodule(m::AbstractAlgebra.FPModule{T}, gens::Vector{Any}) where T <: RingElement
   return Submodule(m, elem_type(m)[])
end

function Submodule(m::AbstractAlgebra.FPModule{T}, subs::Vector{Submodule{T}}) where T <: RingElement
   gens = vcat((generators(s) for s in subs)...)
   return Submodule(m, gens)
end

end #### Generic submodule ####

import .Generic: elem_type

function FreeModule(R::Ring, rank::Int)
   Generic.FreeModule(R, rank)
end

function Submodule(m::FPModule{T}, gens::Vector{<:FPModuleElem{T}}) where T <: RingElement
   Generic.Submodule(m, gens)
end

function Submodule(m::FPModule{T}, subs::Vector{<:Generic.Submodule{T}}) where T <: RingElement
   Generic.Submodule(m, subs)
end

function Submodule(m::FPModule{T}, subs::Vector{S}) where {T <: RingElement, S}
   return Submodule(m, elem_type(m)[])
end

elem_type(::Integers{T}) where T <: Integer = T

ZZ = Integers{BigInt}()

function test_submodule()
#   M = FreeModule(ZZ, 2)  # <-----------------

   M = FreeModule(ZZ, 3)
   nsubs = 0
   subs = [Submodule(M, [M(BigInt[1, 2, 3])]) for i in 1:nsubs]
   N = Submodule(M, subs)
end

end # module
