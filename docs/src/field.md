# Field functionality

## Abstract types for rings

All field types in AbstractAlgebra belong to the `Field` abstract
type and field elements belong to the `FieldElem` abstract type.

As Julia types cannot belong to our `FieldElem` type hierarchy, we also
provide the union type `FieldElement` which includes `FieldElem` in union with
the Julia types `Rational` and `AbstractFloat`.

Note that

```julia
Field <: Ring
FieldElem <: RingElem
FieldElement <: RingElement
```

Of course all `Ring` functionality is available for AbstractAlgebra fields and
their elements.

## Functions for types and parents of fields

```julia
characteristic(R::MyParent)
```

Return the characteristic of the field. If the characteristic is not known, an
exception is raised.

## Basic functions

```julia
isunit(f::MyElem)
```

Return `true` if the given element is invertible, i.e. nonzero in the field.

