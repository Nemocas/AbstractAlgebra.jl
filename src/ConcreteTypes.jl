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
#   Laurent polynomial ring
#
###############################################################################

@attributes mutable struct LaurentPolyRing{T <: RingElement, PR <: PolyRing{T}} <: Ring
  polyring::PR

  function LaurentPolyRing(pr::PR, cached::Bool = true) where {T <: RingElement, PR <: PolyRing{T}}
    return get_cached!(LaurentPolyRingID, pr, cached) do
      new{T, PR}(pr)
    end::LaurentPolyRing{T, PR}
  end
end

const LaurentPolyRingID = CacheDictType{Ring, LaurentPolyRing}()

mutable struct LaurentPolyRingElem{T <: RingElement, PE <: PolyRingElem{T}, LR <: LaurentPolyRing{T}} <: RingElem
  parent::LR
  poly::PE
  mindeg::Int

  # Represents `poly * x^mindeg`, where `x` is the generator of the parent.
  # Not normalized: `LaurentPolyRingElem(parent, poly*x^i, mindeg-i)` is the
  # same Laurent polynomial for every integer `i`.

  function LaurentPolyRingElem(parent::LR, poly::PE, mindeg::Int = 0) where {
                               T <: RingElement, PE <: PolyRingElem{T}, LR <: LaurentPolyRing{T}}
    return new{T, PE, LR}(parent, poly, mindeg)
  end
end

###############################################################################
#
#   Multivariate Laurent polynomial ring
#
###############################################################################

@attributes mutable struct LaurentMPolyRing{T <: RingElement, PR <: MPolyRing{T}} <: Ring
  mpolyring::PR

  function LaurentMPolyRing(pr::PR, cached::Bool = true) where {T <: RingElement, PR <: MPolyRing{T}}
    return get_cached!(LaurentMPolyRingID, pr, cached) do
      new{T, PR}(pr)
    end::LaurentMPolyRing{T, PR}
  end
end

const LaurentMPolyRingID = CacheDictType{Ring, LaurentMPolyRing}()

mutable struct LaurentMPolyRingElem{T <: RingElement, PE <: MPolyRingElem{T}, LR <: LaurentMPolyRing{T}} <: RingElem
  parent::LR
  mpoly::PE            # not necessarily owned by object
  mindegs::Vector{Int} # meaning ditto, vector not necessarily owned by object

  function LaurentMPolyRingElem(parent::LR, poly::PE, mindegs::Vector{Int} = zeros(Int, nvars(parent))) where {
                                T <: RingElement, PE <: MPolyRingElem{T}, LR <: LaurentMPolyRing{T}}
    return new{T, PE, LR}(parent, poly, mindegs)
  end
end

###############################################################################
#
#   Universal ring
#
###############################################################################

@attributes mutable struct UniversalRing{T <: RingElem, U <: RingElement} <: Ring
  base_ring::Ring

  function UniversalRing(R::Ring)
    @req allunique(symbols(R)) "Variable names must be unique"
    return new{elem_type(R), coefficient_type(R)}(R)
  end
end

mutable struct UniversalRingElem{T <: RingElem, U <: RingElement} <: RingElem
  data::T
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
