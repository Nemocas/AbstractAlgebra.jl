###############################################################################
#
#   Abstract types
#
###############################################################################

abstract Pari

abstract Flint

abstract Antic

abstract Generic

abstract Collection{T}

abstract Ring{T} <: Collection{T}

abstract Field{T} <: Ring{T}

abstract CollectionElem

abstract RingElem <: CollectionElem

abstract FieldElem <: RingElem

abstract PolyElem{T} <: RingElem

abstract ResidueElem{T} <: RingElem

abstract FractionElem{T} <: FieldElem

abstract SeriesElem{T} <: RingElem

# not always mathematical ring elements
abstract MatElem <: RingElem

abstract FiniteFieldElem <: FieldElem

abstract NumberFieldElem <: FieldElem

abstract MaximalOrderElem <: RingElem

