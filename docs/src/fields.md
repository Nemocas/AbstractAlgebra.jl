# Field Interface

AbstractAlgebra.jl generic code makes use of a standardised set of functions which it
expects to be implemented for all fields. Here we document this interface. All libraries
which want to make use of the generic capabilities of AbstractAlgebra.jl must supply
all of the required functionality for their fields.

## Types

Most fields must supply two types:
  - a type for the parent object (representing the field itself)
  - a type for elements of that field

For example, the generic fraction field type in AbstractAlgebra.jl provides two 
types in `generic/GenericTypes.jl`: 

  - `Generic.FracField{T}` for the parent objects
  - `Generic.Frac{T}` for the actual fractions

The parent type must belong to `AbstractAlgebra.Field` and the element type must belong
to `AbstractAlgebra.FieldElem`. Of course, the types may belong to these abstract types
transitively.

For parameterised fields, we advise that the types of both the parent objects and
element objects to be parameterised by the types of the elements of the base ring.

There can be variations on this theme: e.g. in some areas of mathematics there is a
notion of a coefficient domain, in which case it may make sense to parameterise all
types by the type of elements of this coefficient domain. But note that this may have
implications for the ad hoc operators one might like to explicitly implement.

## Parent object caches

In many cases, it is desirable to have only one object in the system to represent each
field. This means that if the same field is constructed twice, elements of the two fields
will be compatible as far as arithmetic is concerned.

In order to facilitate this, global caches of fields are stored in AbstractAlgebra.jl,
usually implemented using dictionaries. For example, the `Generic.FracField` parent
objects are looked up in a dictionary `FracDict` to see if they have been previously
defined.

Whether these global caches are provided or not, depends on both mathematical and
algorithmic considerations. E.g. in the case of number fields, it isn't desirable to
identify all number fields with the same defining polynomial, as they may be considered
with distinct embeddings into one another. In other cases, identifying whether two
fields are the same may be prohibitively expensive. Generally, it may only make sense
algorithmically to identify two fields if they were constructed from identical data.

If a global cache is provided, it must be optionally possible to construct the parent
objects without caching. This is done by passing a boolean value `cached` to the inner
constructor of the parent object. See generic/GenericTypes.jl` for examples of how to
construct and handle such caches.

## Required functions for all fields

In the following, we list all the functions that are required to be provided for fields
in AbstractAlgebra.jl or by external libraries wanting to use AbstractAlgebra.jl.

We give this interface for fictitious types `MyParent` for the type of the field parent
object `R` and `MyElem` for the type of the elements of the field.

Note that generic functions in AbstractAlgebra.jl may not rely on the existence of
functions that are not documented here. If they do, those functions will only be
available for fields that implement that additional functionality, and should be
documented as such.

In the first place, all fields are rings and therefore any field type must implement
all of the Ring interface. The functionality below is in addition to this basic
functionality.

### Data type and parent object methods

```julia
characteristic(R::MyParent)
```

Return the characteristic of the field.

### Basic manipulation of rings and elements

```julia
isunit(f::MyElem)
```

Return `true` if the given element is invertible, i.e. nonzero in the field.

### Inversion

```julia
inv(f::MyElem)
```

Return the inverse of the given element in the field. If $f = 0$, an error is thrown.

