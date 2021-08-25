mutable struct PolyRing{T <: RingElement} <: AbstractAlgebra.PolyRing{T}
   base_ring::Ring
   S::Symbol

   function PolyRing{T}(R::Ring, s::Symbol, cached::Bool = true) where T <: RingElement
      return get_cached!(PolyID, (R, s), cached) do
         new{T}(R, s)
      end::PolyRing{T}
   end
end

const PolyID = CacheDictType{Tuple{Ring, Symbol}, Ring}()

mutable struct Poly{T <: RingElement} <: AbstractAlgebra.PolyElem{T}
   coeffs::Vector{T}
   length::Int
   parent::PolyRing{T}

   Poly{T}() where T <: RingElement = new{T}(Array{T}(undef, 0), 0)

   function Poly{T}(b::Vector{T}) where T <: RingElement
      z = new{T}(b)
      z.length = normalise(z, length(b))
      return z
   end

   Poly{T}(a::T) where T <: RingElement = iszero(a) ? new{T}(Array{T}(undef, 0), 0) : new{T}([a], 1)
end

mutable struct FracField{T <: RingElem} <: AbstractAlgebra.FracField{T}
   base_ring::Ring

   function FracField{T}(R::Ring, cached::Bool = true) where T <: RingElem
      return get_cached!(FracDict, R, cached) do
         new{T}(R)
      end::FracField{T}
   end
end

const FracDict = CacheDictType{Ring, Ring}()

mutable struct Frac{T <: RingElem} <: AbstractAlgebra.FracElem{T}
   num::T
   den::T
   parent::FracField{T}

   Frac{T}(num::T, den::T) where T <: RingElem = new{T}(num, den)
end

struct RationalFunctionField{T <: FieldElement} <: AbstractAlgebra.Field
   S::Symbol
   fraction_field::FracField{<:PolyElem{T}}
   base_ring::Field

   function RationalFunctionField{T}(k::Field, frac_field::FracField{<:PolyElem{T}}, sym::Symbol, cached::Bool = true) where T <: FieldElement
      return get_cached!(RationalFunctionFieldDict, (k, sym), cached) do
         U = elem_type(k)
         new{U}(sym, frac_field, k)
      end::RationalFunctionField{T}
   end
end

const RationalFunctionFieldDict = CacheDictType{Tuple{Field, Symbol}, Field}()

mutable struct Rat{T <: FieldElement} <: AbstractAlgebra.FieldElem
   d::Frac{<:PolyElem{T}}
   parent::RationalFunctionField{T}

   Rat{T}(f::Frac{<:PolyElem{T}}) where T <: FieldElement = new{T}(f)
end

