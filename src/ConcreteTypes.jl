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

@attributes mutable struct UniversalRing{T <: RingElement} <: Ring
  base_ring::Ring

  function UniversalRing{T}(R::Ring) where {T<:RingElement}
    @assert elem_type(R) == T
    return new{T}(R)
  end
end

mutable struct UniversalRingElem{T <: RingElement} <: RingElem
  p::T
  parent::UniversalRing{T}
end

###############################################################################
#
#   Universal polynomial ring
#
###############################################################################

const UniversalPolyRing{T} = UniversalRing{<:MPolyRingElem{T}}

const UniversalPolyRingElem{T} = UniversalRingElem{<:MPolyRingElem{T}}

const UnivPoly{T} = UniversalPolyRingElem{T}

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
