module AbstractAlgebra

export QQ, FreeModule, VectorSpace, Submodule, Generic, test_submodule

abstract type NCRingElem end

abstract type RingElem <: NCRingElem end

abstract type FieldElem <: RingElem end

abstract type NCRing end

abstract type Ring <: NCRing end

abstract type Field <: Ring end

abstract type Group end

abstract type Module{T} <: Group end

abstract type FPModule{T} <: Module{T} end

abstract type ModuleElem{T} end

abstract type FPModuleElem{T} <: ModuleElem{T} end

mutable struct Integers{T <: Integer} <: Ring
end

mutable struct Rationals{T <: Integer} <: Field
end

const RingElement = Union{RingElem, Integer, Rational, AbstractFloat}

const FieldElement = Union{FieldElem, Rational, AbstractFloat}

###############################################################################
#
#   Generic submodule
#
###############################################################################

module Generic

import Base: rand

import AbstractAlgebra: Rationals, NCRing, NCRingElem, Ring, RingElem,
       RingElement

using AbstractAlgebra

export Submodule, submodule_elem, ngens, FreeModule, free_module_elem

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
end

parent_type(::Type{free_module_elem{T}}) where T <: RingElement = FreeModule{T}

base_ring(M::FreeModule{T}) where T <: Union{RingElement, NCRingElem} = M.base_ring::parent_type(T)

elem_type(::FreeModule{T}) where T <: Union{RingElement, NCRingElem} = free_module_elem{T}

rank(M::FreeModule{T}) where T <: Union{RingElement, NCRingElem} = M.rank

ngens(M::FreeModule{T}) where T <: Union{RingElement, NCRingElem} = M.rank

function (M::FreeModule{T})(a::Vector{T}) where T <: Union{RingElement, NCRingElem}
   return free_module_elem{T}(M, a)
end

function FreeModule(R::NCRing, rank::Int)
   T = elem_type(R)
   return FreeModule{T}(R, rank)
end

base_ring(N::Submodule{T}) where T <: RingElement = N.base_ring

generators(N::Submodule{T}) where T <: RingElement = N.gens::Vector{elem_type(N.m)}

function Submodule(m::AbstractAlgebra.FPModule{T}, gens::Vector{Any}) where T <: RingElement
   return Submodule(m, elem_type(m)[])
end

function Submodule(m::AbstractAlgebra.FPModule{T}, subs::Vector{Submodule{T}}) where T <: RingElement
   gens = vcat((generators(s) for s in subs)...)
   return Submodule(m, gens)
end

end # generic

import .Generic: elem_type, parent_type

function FreeModule(R::NCRing, rank::Int; cached::Bool = true)
   Generic.FreeModule(R, rank)
end

function Submodule(m::Module{T}, gens::Vector{<:ModuleElem{T}}) where T <: RingElement
   Generic.Submodule(m, gens)
end

function Submodule(m::Module{T}, subs::Vector{<:Generic.Submodule{T}}) where T <: RingElement
   Generic.Submodule(m, subs)
end

elem_type(::Rationals{T}) where T <: Integer = Rational{T}

parent_type(::Type{Rational{T}}) where T <: Integer = Rationals{T}

QQ = Rationals{BigInt}()

function test_submodule()
#   M = FreeModule(QQ, 2)

   M = FreeModule(QQ, 5)
   nsubs = 0
   subs = [Submodule(M, [M(BigInt[1, 2, 3, 4, 5])]) for i in 1:nsubs]
   N = Submodule(M, subs)
end

end # module
