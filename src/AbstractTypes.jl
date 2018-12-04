###############################################################################
#
#   Abstract types
#
###############################################################################

# broad mathematical domains
# these contain the type classes of parent objects

abstract type Set end

abstract type Group <: Set end

abstract type NCRing <: Set end

abstract type Ring <: NCRing end

abstract type Field <: Ring end

# elements of mathematical domains

abstract type SetElem end

abstract type GroupElem <: SetElem end

abstract type NCRingElem <: SetElem end

abstract type RingElem <: NCRingElem end

abstract type FieldElem <: RingElem end

# parameterized domains

abstract type Module{T} <: Group end

abstract type Ideal{T} <: Set end

# elements of parameterised domains

abstract type ModuleElem{T} <: GroupElem end

abstract type IdealElem{T} <: SetElem end

abstract type Map{D, C, S, T} <: SetElem end

abstract type SetMap end

abstract type FunctionalMap <: SetMap end

abstract type IdentityMap <: SetMap end

# rings, fields etc, parameterised by an element type
# these are the type classes of different kinds of
# mathematical rings/fields/etc, which have a base ring,
# and for which a generic implementation is possible
# over that base ring

abstract type PolyRing{T} <: Ring end

abstract type NCPolyRing{T} <: NCRing end

abstract type MPolyRing{T} <: Ring end

abstract type SeriesRing{T} <: Ring end

abstract type ResRing{T} <: Ring end

abstract type ResField{T} <: Field end

abstract type FracField{T} <: Field end

abstract type MatSpace{T} <: Module{T} end

abstract type MatAlgebra{T} <: NCRing end

# mathematical objects parameterised by an element type
# these are the type classes of mathematical objects
# that have some kind of base ring, and a generic
# implementation is meaningful over that base ring

abstract type PolyElem{T} <: RingElem end

abstract type NCPolyElem{T} <: NCRingElem end

abstract type MPolyElem{T} <: RingElem end

abstract type ResElem{T} <: RingElem end

abstract type ResFieldElem{T} <: FieldElem end

abstract type FracElem{T} <: FieldElem end

abstract type SeriesElem{T} <: RingElem end

abstract type RelSeriesElem{T} <: SeriesElem{T} end

abstract type AbsSeriesElem{T} <: SeriesElem{T} end

abstract type MatElem{T} <: ModuleElem{T} end

abstract type MatAlgElem{T} <: NCRingElem end

# additional abstract types for parents, added ad hoc to form
# collections of types as needed by applications

abstract type FinField <: Field end    # for fq, fq_nmod, etc

# additional abstract types for elements, added ad hoc to form
# collections of types as needed by applications

abstract type FinFieldElem <: FieldElem end # for fq, fq_nmod, etc

################################################################################
#
#   Promotion system
#
# The promote_rule functions are not extending Base.promote_rule. The
# AbstractAlgebra promotion system is orthogonal to the built-in julia promotion
# system. The julia system assumes that whenever you have a method signature of
# the form Base.promote_rule(::Type{T}, ::Type{S}) = R, then there is also a
# corresponding Base.convert(::Type{R}, ::T) and similar for S. Since we
# cannot use the julia convert system (we need an instance of the type and not
# the type), we cannot use the julia promotion system.
#
# The AbstractAlgebra promotion system is used to define catch all functions
# for arithmetic between arbitrary ring elements.
#
################################################################################

promote_rule(T, U) = Union{}
