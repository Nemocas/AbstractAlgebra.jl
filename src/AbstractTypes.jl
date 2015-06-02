###############################################################################
#
#   Abstract types
#
###############################################################################

abstract Pari

abstract Flint

abstract Antic

abstract Generic

abstract Set{T}

abstract Ring{T} <: Set{T}

abstract Field{T} <: Ring{T}

abstract RingElem

abstract FieldElem <: RingElem

abstract PolyElem{T} <: RingElem

abstract PowerSeriesElem <: RingElem

# not always mathematical ring elements
abstract MatElem <: RingElem

abstract FiniteFieldElem <: FieldElem

abstract NumberFieldElem <: FieldElem

abstract MaximalOrderElem <: RingElem

