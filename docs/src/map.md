# Map interface

Maps in Nemo model maps on sets $f : D \to C$ for some domain $D$ and codomain $C$, which
have no real limitations except that elements of the codomain and domain be represented
by element objects in the system..

Maps $f : D \to C$ in AbstractAlgebra are modeled by Julia objects that are able to be
called on a single element $d \in D$ of the domain to yield an element $f(d) \in C$ of the
codomain. We say that the map is being applied.

Maps can be constructed from Julia functions, or they can be represented by some other
kind of data, e.g. a matrix, or built up from other maps.

In the following, we will always use the word "function" to mean a Julia function, and
reserve the word "map" for a map on sets, whether mathematically, or as an object in the
system..

Maps in AbstractAlgebra have a domain and codomain, can be applied, composed and composed
with the identity map (assuming its domain is compatible). Various special kinds of map
provide more functionality.

## Parent objects

Maps in AbstractAlgebra currently don't have parents. This will change later when
AbstractAlgebra has a category system, so that the parent of a map can be some sort of
Hom set.

## Map classes

All maps in AbstractAlgebra belong to a class of maps. The classes are modeled as
abstract types that lie in a hierarchy, inheriting from `SetMap` at the top of the
hierarchy. Other classes that inherit from `SetMap` are `FunctionalMap` for maps that
are constructed from a Julia function (or closure), and `IdentityMap` for the class of
the identity maps within the system.

One might naturally assume that map types belong directly to these classes in the way
that types of other objects in the system belong to abstract types in the AbstractAlgebra
type hierarchy. However, in order to provide an extensible system, this is not the case.

Instead, a map type `MyMap` will belong to an abstract type of the form
`Map{D, C, T, MyMap}`, where `D` is the type of the object representing the domain of
the map type (this can also be an abstract type, such as `Group`), `C` is the type of
the object representing the codomain of the map type and `T` is the map class that `MyMap`
belongs to, e.g. `SetMap` or `FunctionalMap`.

Because a four parameter type system becomes quite cumbersome to use, we provide a number
of functions for referring to collections of map types.

If writing a function that accepts any map type, one makes the type of its argument
belong to `Map`. For example `f(M::Map) = 1`.

If writing a function that accepts any map from a domain of type `D` to a codomain of
type `C`, one makes writes for example `f(M::Map{D, C}) = 2`. Note that `D` and `C` can
be abstract types, such as `Group`, but otherwise must be the types of the parent objects
representing the domain and codomain.

A function that accepts any map belonging to a given map class might be written as
`f(M::Map(FunctionalMap)) = 3` or `f(M::Map(FunctionalMap){D, C}) = 4` for example, where
`D` and `C` are the types of the parent objects for the domain and codomain.

Finally, if a function should only work for a map of a given map type `MyMap`, say, one
writes this `f(M::Map(MyMap))` or `f(M::Map(MyMap){D, C}`, where as usual `D` and `C`
are the types of the domain and codomain parent objects.

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

As mentioned above, to write a function that only accepts maps of type `MyMap`, one
writes the functions as follows:

```julia
function my_fun(M::Map(MyMap))
```

The `Map` function then computes the correct type to use, which is actually not `MyMap`
if all features of the generic Map infrastructure are required. It is bad practice to
write functions for `MyMap` directly instead of `Map(MyMap)`, since other users will be
unable to use generic constructions over the map type `MyMap`.

## Getters and setters

When writing new map types, it is very important to define getters and setters of the
fields of the new map type, rather than to access them directly.

Let us suppose that the `MyMap` type has a field called `foo`. Rather than access this
field by writing `M.foo`, one must access it using `foo(M)` (at least until Julia 1.1).

If such a getter only needs to access the field `foo` of `M`, there is a standard way of
defining such a getter and setter when defining a new map type.

```julia
foo(M::Map(MyMap)) = get_field(M, :foo)
```

To set a field of a map, one needs a setter, which can be implemented as follows:

```julia
set_foo!(M::Map(MyMap), a) = set_field(M, :foo, a)
```

In general, setters should be used rarely for map types.

Note that by providing getter and setter functions, map types need not even contain
fields with the given name. For example, for a `MyMap` map type for maps between integers, one does not wish to explicitly store the domain and codomain in `MyMap`. Instead, we
can define the getter functions `domain` and `codomain` to return `JuliaZZ` for any
`MyMap` object.

```julia
domain(M::Map(MyMap)) = JuliaZZ
codomain(M::Map(MyMap)) = JuliaZZ
```

## Required functionality for maps

All map types must implement a standard interface, which we specify here.

We will define this interface for a custom map type `MyMap` belonging to `Map(SetMap)`,
`SetMap` being the map class that all maps types belong to.

Note that map types do not need to contain any specific fields, but must provide accessor
functions (getters and setters) in the manner described above.

The required accessors for map types of class `SetMap` are as follows.

```julia
domain(M::Map(MyMap))
codomain(M::Map(MyMap))
```

Return the domain and codomain parent objects respectively, for the map $M$. It is only
necessary to define these functions if the map type `MyMap` does not contain fields
`domain` and `codomain` containing these parent objects.

It is also necessary to be able to apply a map. This amounts to overloading the call
method for objects belonging to `Map(MyMap)`.

```julia
(M::Map(MyMap)(a))
```

Apply the map `M` to the element `a` of the domain of `M`. Note that it is usual to
add a type assertion to the return value of this function, asserting that the return
value has type `elem_type(C)` where `C` is the type of the codomain parent object.

## Optional functionality for maps

The Generic module in AbstractAlgebra automatically provides certain functionality for
map types, assuming that they satisfy the full interface described above.

However, certain map types or map classes might like to provide their own implementation
of this functionality, overriding the generic functionality.

We describe this optional functionality here.

### Show method

Custom map types may like to provide a custom `show` method if the default of displaying
the domain and codomain of the map is not sufficient.

```julia
show(io::IO, M::Map(MyMap))
```

### Identity maps

There is a concrete map type `Generic.IdentityMap{D}` for the identity map on a given
domain. Here `D` is the type of the object representing that domain.

`Generic.IdentityMap` belongs to the supertype
`Map{D, C, AbstractAlgebra.IdentityMap, IdentityMap}`.

Note that the map class is also called `IdentityMap`. It is an abstract type, whereas
`Generic.IdentityMap` is a concrete type in the Generic module.

An identity map has the property that when composed with any map whose domain or codomain
is compatible, that map will be returned as the composition. Identity maps can therefore
serve as a starting point when building up a composition of maps, starting an
identity map.

We do not cached identity maps in the system, so that if more than one is created on the
same domain, there will be more than one such map in the system. This underscores the
fact that there is in general no way for the system to know if two maps compose to give
an identity map, and therefore the only two maps that can be composed to give an identity
map are identity maps on the same domain.

To construct an identity map for a given domain, specified by a parent object `R`, say,
we have the following function.

```julia
identity_map(R::D) where D <: Set
```

Return an identity map on the domain $R$.

Of course there is nothing stopping a map type or class from implementing its own
identity map type, and defining composition of maps of the same kind with such an
identity map. In such a case, the class of such an identity map type must belong to
`AbstractAlgebra.IdentityMap` so that composition with other map types still works.

### Composition of maps

Any two compatible maps in AbstractAlgebra can be composed and any composition can be
applied.

In order to facilitate this, the Generic module provides a type `CompositionMap{D, C}`,
which contains two maps `map1` and `map2`, corresponding to the two maps to be applied
in a composition, in the order they should be applied.

To construct a composition map from two existing maps, we have the following function:

```julia
compose(f::Map{D, U}, g::Map{U, C}) where {D, U, C}
```

Compose the two maps $f$ and $g$, i.e. return the map $h$ such that $h(x) = g(f(x))$.

As a shortcut for this function we have the following operator:

```julia
*(f::Map{D, U}, g::Map{U, C}) where {D, U, C} = compose(f, g)
```

Note the order of composition. If we have maps $f : X \to Y$, $g : Y \to Z$ the correct
order of the maps in this operator is `f*g`, so that `(f*g)(x) = g(f(x))`.

This is chosen so that for left $R$-module morphisms represented by a matrix, the order
of matrix multiplication will match the order of composition of the corresponding
morphisms.

Of course, a custom map type or class of maps can implement its own composition type
and compose function.

This is the case with the `FunctionalMap` class for example, which caches the Julia
function/closure corresponding to the composition of two functional maps. As this cached
function needs to be stored inside the composition, a special type is necessary for the
composition of two functional maps.

By default, `compose` will check that the two maps are composable, i.e. the codomain
of the first map matches the domain of the second map. This is implemented by the
following function:

```julia
check_composable(f::Map{D, U}, g::Map{U, C})
```
Raise an exception if the codomain of $f$ doesn't match the domain of $g$.

Note that composite maps should keep track of the two maps they were constructed from. To
access these maps, the following functions are provided:

```julia
map1(f::CompositeMap)
map2(f::CompositeMap)
```

Any custom composite map type must also provide these functions for that map type, even
if there exist fields with those names. This is because there is no common map class for
all composite map types. Therefore the Generic system cannot provide fallbacks for all
such composite map types.
