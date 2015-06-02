###############################################################################
#
#   Abstract types
#
###############################################################################

abstract Ring

abstract Field <: Ring

abstract RingElem

abstract FieldElem <: RingElem

abstract PolyElem{T} <: RingElem

abstract PowerSeriesElem <: RingElem

# not always mathematical ring elements
abstract MatElem <: RingElem

abstract FiniteFieldElem <: FieldElem

abstract NumberFieldElem <: FieldElem

abstract MaximalOrderElem <: RingElem

