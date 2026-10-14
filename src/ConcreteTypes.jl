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
#   MatRing / MatRingElem
#
###############################################################################

struct MatRing{T <: NCRingElement} <: NCRing
   n::Int
   base_ring::NCRing

   function MatRing{T}(R::NCRing, n::Int) where T <: NCRingElement
      new{T}(n, R)
   end
end

struct MatRingElem{T <: NCRingElement} <: NCRingElem
   data::MatElem{T}

  function MatRingElem{TT}(A::MatElem{TT}) where TT <: NCRingElement
    nrows(A) == ncols(A) || error("Matrix must be square")
    return new{TT}(A)
   end
end

MatRingElem(A::MatElem{TT}) where TT <: NCRingElement = MatRingElem{TT}(A)

function MatRingElem(R::NCRing, n::Int, A::Vector{T}) where T <: NCRingElement
   t = matrix(R, n, n, A)
   return MatRingElem(t)
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
