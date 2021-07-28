# Maps

Maps in AbstractAlgebra can be constructed from Julia functions, or they can be
represented by some other kind of data, e.g. a matrix, or built up from other maps.

In the following, we will always use the word "function" to mean a Julia function, and
reserve the word "map" for a map on sets, whether mathematically, or as an object in the
system.

## Parent objects

Maps in AbstractAlgebra currently don't have parents.

## Map classes

Each map in AbstractAlgebra has a map class. Currently the supported map classes are
`SetMap` which is at the root of the tree of map classes and `FunctionalMap` for maps
which are built on Julia functions and `IdentityMap`.

Other map classes can be added by inheriting from one of these. A map class is in
practice nothing more than a Julia abstract type.

One might naturally assume that map types belong directly to these classes in the way
that types of other objects in the system belong to abstract types in the AbstractAlgebra
type hierarchy. However, in order to provide an extensible system, this is not the case.

Instead, a map type `MyMap` say will belong to an abstract type of the form
`Map{D, C, T, MyMap}`, where `D` is the type of the object representing the domain of
the map type (this can also be an abstract type, such as `Group`), `C` is the type of
the object representing the codomain of the map type and `T` is the map class.

Firstly note that all maps in the system belong to the abstract type `Map`. Moreover,
maps from a domain of type `D` to a codomain of type `C` all belong to `Map{D, C}`.

However, because a four parameter type system becomes quite cumbersome to use, we
provide a number of functions for referring to maps by their map class or by their
specific map type.

Maps with map class `FunctionalMap` all belong to `Map(FunctionalMap))`. One can also
restrict the domain and codomain by writing `Map(FunctionalMap){D, C}`.

Finally, if a function should only work for a map of a the very specific map type
`MyMap`, one writes `M::Map(MyMap))` or `M::Map(MyMap){D, C}` if one wishes to restrict
the domain and codomain types.

## Implementing new map types

There are two common kinds of map type that developers will need to write. The first
has a fixed domain and codomain, and the second is a type parameterised by the types of
the domain and codomain. We give two simple examples here of how this might look.

In the case of fixed domain and codomain, e.g. `Integers{BigInt}`, we would write it as
follows:

```julia
mutable struct MyMap <: Map{Integers{BigInt}, Integers{BigInt}, SetMap, MyMap}
   # some data fields
end
```

In the case of parameterisation by the type of the domain and codomain:

```julia
mutable struct MyMap{D, C} <: Map{D, C, SetMap, MyMap}
   # some data fields
end
```

## Getters and setters

When writing new map types, it is very important to define getters and setters of the
fields of the new map type, rather than to access them directly.

Let us suppose that the `MyMap` type has a field called `foo`. Rather than access this
field by writing `M.foo`, one must access it using `foo(M)`.

If such a getter only needs to access the field `foo` of `M`, there is a standard way of
defining such a getter and setter when defining a new map type.

```julia
foo(M::Map(MyMap)) = get_field(M, :foo)
```

To set a field of a map, one needs a setter, which can be implemented as follows:

```julia
set_foo!(M::Map(MyMap), a) = set_field(M, :foo, a)
```

## Map functions

```julia
domain(M::Map(MyMap))
codomain(M::Map(MyMap))
```

Return the domain and codomain parent objects respectively, for the map $M$.

All maps can be applied to elements in the domain.

```julia
(M::Map(MyMap)(a))
```

Apply the map `M` to the element `a` of the domain of `M`.

### Identity maps

There is a concrete map type `Generic.IdentityMap{D}` for the identity map on a given
domain. Here `D` is the type of the object representing that domain.

`Generic.IdentityMap` belongs to the map type
`Map{D, C, AbstractAlgebra.IdentityMap, IdentityMap}`.

An identity map has the property that when composed with any map whose domain or codomain
is compatible, that map will be returned as the composition. Identity maps can therefore
serve as a starting point when building up a composition of maps, starting an
identity map.

To construct an identity map for a given domain, specified by a parent object `R`, say,
we have the following function.

```julia
identity_map(R::Set)
```

### Composition of maps

Any two compatible maps in AbstractAlgebra can be composed and any composition can be
applied.

To construct a composition map from two existing maps, we have the following function:

```julia
compose(f::Map{D, U}, g::Map{U, C}) where {D, U, C}
```

As a shortcut for this function we have the following operator:

```julia
*(f::Map{D, U}, g::Map{U, C}) where {D, U, C} = compose(f, g)
```

!!! note

    Observe the order of composition. If we have maps $f : X \to Y$, $g : Y \to Z$ the correct
    order of the maps in this operator is `f*g`, so that `(f*g)(x) = g(f(x))`.

    This is chosen so that for left $R$-module morphisms represented by a matrix, the order
    of matrix multiplication will match the order of composition of the corresponding
    morphisms.

