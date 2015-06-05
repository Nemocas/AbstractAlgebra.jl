###############################################################################
#
#   Abstract types
#
###############################################################################

# libraries

abstract Pari

abstract Flint

abstract Antic

abstract Generic

# mathematical domains, parameterised by a library
# these are the type classes of parent objects

   abstract Collection{T}

   abstract Ring{T} <: Collection{T}

   abstract Field{T} <: Ring{T}

# elements of mathematical domains

   abstract CollectionElem

   abstract RingElem <: CollectionElem

   abstract FieldElem <: RingElem

# mathematical objects parameterised by an element type
# these are the type classes of mathematical objects
# that have some kind of base ring, and a generic 
# implementation is meaningful over that base ring

   abstract PolyElem{T} <: RingElem

   abstract ResidueElem{T} <: RingElem

   abstract FractionElem{T} <: FieldElem

   abstract SeriesElem{T} <: RingElem

   # not always mathematical ring elements
   # later we'll maybe distinguish MatAlgebraElem, MatModuleElem
   abstract MatElem{T} <: RingElem

# leaf objects, with no parameterisation
# these are also type classes of mathematical objects
# usually provided by a C library and not by generic

   abstract IntegerRingElem <: RingElem

   abstract FiniteFieldElem <: FieldElem

   abstract NumberFieldElem <: FieldElem

   abstract MaximalOrderElem <: RingElem

   abstract PadicFieldElem <: FieldElem

