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

@attributes mutable struct GFField{T <: Integer} <: FinField
  p::T

  function GFField{T}(p::T; cached::Bool = true) where T <: Integer
    return get_cached!(GFFieldID, (T, p), cached) do
      new{T}(p)
    end::GFField{T}
  end
end

const GFFieldID = CacheDictType{Tuple{DataType, Integer}, Field}()

struct GFElem{T <: Integer} <: FinFieldElem
  d::T
  parent::GFField{T}
end
