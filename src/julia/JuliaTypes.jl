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

   function GFField{T}(p::T; cached::Bool = true) where T <: Integer
      if cached && haskey(GFFieldID, (T, p))
         z = GFFieldID[T, p]::GFField{T}
      else
         z = new{T}(p)
         if cached
           GFFieldID[T, p] = z
         end
      end
      return z
   end
end

const GFFieldID = Dict{Tuple{DataType, Integer}, Field}()

struct GFElem{T <: Integer} <: FinFieldElem
   d::T
   parent::GFField{T}
end
