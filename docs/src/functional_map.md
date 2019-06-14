```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Functional maps

A functional map in AbstractAlgebra is a map which can be applied by evaluating a Julia
function or closure. It is represented by a map object that contains such a
function/closure, usually in a field called `image_fn`.

All functional maps belong to the map class `AbstractAlgebra.FunctionalMap`.

A generic concrete type `Generic.FunctionalMap` is provided by the Generic module to
implement a generic functional map type. This allows for functional maps that contain no
extra data, other than a Julia function/closure.

Custom map types can also be defined which have map class `AbstractAlgebra.FunctionalMap`.

## Functional map interface

All functional map types must define their supertypes as in the following example:

```julia
mutable struct MyFunctionalMap{D, C} <: Map{D, C, FunctionalMap, MyFunctionalMap}
   # some fields
   image_fn::Function
end
```
Of course `MyFunctionalMap` need not be parameterised if the types `D` and `C` of the
domain and codomain objects are known.

### Required functions for functional maps

The following functions must be defined for all functional map types or classes:

```julia
image_fn(M::Map(MyFunctionalMap))
```

Return the Julia function or closure that corresponds to application of the map $M$. This
function only needs to be provided if this function is not stored in an `image_fn`
field of the `MyFunctionalMap` type.

## Generic functional maps

The Generic module provides a concrete type `FunctionalMap` which merely keeps track of
a Julia function/closure implementing the map.

Such maps can be constructed using the following function:

```julia
map_from_func(R, S, f::Function)
```

Construct the generic functional map with domain and codomain given by the parent objects
$R$ and $S$ corresponding to the Julia function $f$.

**Examples**

```jldoctest
julia> f = map_from_func(ZZ, ZZ, x -> x + 1)
Map with the following data

Domain:
=======
Integers

Codomain:
========
Integers


julia> f(ZZ(2))
3

```

