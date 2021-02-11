###############################################################################
#
#   Integers / Integer
#
###############################################################################

struct Integers{T <: Integer} <: Ring
end

###############################################################################
#
#   Rationals / Rational
#
###############################################################################

struct Rationals{T <: Integer} <: Field
end

###############################################################################
#
#   Floats / AbstractFloat
#
###############################################################################

struct Floats{T <: AbstractFloat} <: Field
end

###############################################################################
#
#   GFField/GFElem
#
###############################################################################

mutable struct GFField{T <: Integer} <: FinField
   p::T

   function GFField{T}(p::T) where T <: Integer
      if haskey(GFFieldID, (T, p))
         z = GFFieldID[T, p]::GFField{T}
      else
         z = new{T}(p)
         GFFieldID[T, p] = z
      end
      return z
   end
end

const GFFieldID = Dict{Tuple{DataType, Integer}, Field}()

struct GFElem{T <: Integer} <: FinFieldElem
   d::T
   parent::GFField{T}
end

###############################################################################
#
#   Unions of AbstactAlgebra abstract types and Julia types
#
###############################################################################

const RingElement   = Union{RingElem,   Integer, Rational, AbstractFloat}
const NCRingElement = Union{NCRingElem, Integer, Rational, AbstractFloat}

const FieldElement = Union{FieldElem, Rational, AbstractFloat}
