###############################################################################
#
#   Mat space
#
###############################################################################

struct MatSpace{T <: NCRingElement} <: Module{T}
  base_ring::NCRing
  nrows::Int
  ncols::Int

  function MatSpace{T}(R::NCRing, r::Int, c::Int, cached::Bool = true) where T <: NCRingElement
     # TODO/FIXME: `cached` is ignored and only exists for backwards compatibility
     @assert elem_type(R) === T
     (r < 0 || c < 0) && error("Dimensions must be non-negative")
     return new{T}(R, r, c)
  end
end

###############################################################################
#
#   Universal ring
#
###############################################################################

@attributes mutable struct UniversalRing{T <: RingElem, U <: RingElement} <: Ring
  base_ring::Ring

  function UniversalRing{T, U}(R::Ring) where {T <: RingElem, U <: RingElement}
    @assert elem_type(R) == T
    @assert elem_type(coefficient_ring(R)) == U
    return new{T, U}(R)
  end
end

mutable struct UniversalRingElem{T <: RingElem, U <: RingElement} <: RingElem
  p::T
  parent::UniversalRing{T, U}
end

###############################################################################
#
#   Universal polynomial ring
#
###############################################################################

struct UnivPolyCoeffs{T <: RingElem}
  poly::T
end

struct UnivPolyExponentVectors{T <: RingElem}
  poly::T
end

struct UnivPolyTerms{T <: RingElem}
  poly::T
end

struct UnivPolyMonomials{T <: RingElem}
  poly::T
end
