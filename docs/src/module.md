# Module Interface

All module types in AbstractAlgebra follow the following interface.

Modules can be built over both commutative and noncommutative rings.

## Types and parents

AbstractAlgebra provides two abstract types for modules and their elements:

  * `Module{T}` is the abstract type for module parent types
  * `ModuleElem{T}` is the abstract type for module element types

Note that the abstract types are parameterised. The type `T` should usually be the type
of elements of the ring the module is over.

## Required functionality for modules

In addition to the required functionality for the Ring/NCRing interface (and in the case
of polynomials over a field, the Euclidean Ring interface), the Polynomial Ring interface
has the following required functions.

We suppose that `R` is a fictitious base ring and that `S` is a module over `R` with
parent object `S` of type `MyModule{T}`. We also assume the elements in the module have
type `MyModuleElem{T}`, where `T` is the type of elements of the ring the module is
over.

Of course, in practice these types may not be parameterised, but we use parameterised
types here to make the interface clearer.

Note that the type `T` must (transitively) belong to the abstract type `RingElem` or
`NCRingElem`.

We describe the functionality below for modules over commutative rings, i.e. with
element type belonging to `RingElem`, however similar constructors should be available
for element types belonging to `NCRingElem` instead, if the basse ring is
noncommutative.
