```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Map with inverse

It is not possible to provide generic functionality to invert a map. However, sometimes
one knows an inverse map explicitly and would like to keep track of this.

Recall that as map composition is not commutative, there is a notion of a left inverse
and a right inverse for maps.

To keep track of such inverse maps, AbstractAlgebra provides data types
`Generic.MapWithRetraction` and `Generic.MapWithSection`.

Given a map $f : X \to Y$, a retraction of $f$ is a map $g : Y \to X$ such that
$g(f(x)) = x$ for all $x \in X$.

Given a map $f : X \to Y$, a section of $f$ is a map $g : Y \to X$ such that
$f(g(x)) = x$ for all $y \in Y$.

In AbstractAlgebra, a map with retraction/section is an object containing a pair of maps,
the second of which is a retraction/section of the first.

Maps with retraction/section can be composed, and we also define the inverse of such a
pair to be the map with the pair swapped. Thus the inverse of a map with retraction is
a map with section.

## Map with inverse constructors

To construct a map with retraction/section from a pair of maps, we have the following
functions:

```@docs
map_with_retraction
map_with_section
```

For convenience we allow construction of maps with retraction/section from a pair of
Julia functions/closures.

```@docs
map_with_retraction_from_func
map_with_section_from_func
```

## Functionality for maps with inverses

The following functionality is provided for maps with inverses.

```julia
inv(M::Generic.MapWithRetraction)
inv(M::Generic.MapWithSection)
```

Return the map with the two maps contained in $M$ swapped. In the first case, a
`MapWithSection` is returned. In the second case a `MapWithRetraction` is returned.

To access the two maps stored in a map with retraction/section, we have the following:

```@docs
image_map
retraction_map
section_map
```

**Examples**

```jldoctest
julia> f = map_with_retraction_from_func(x -> x + 1, x -> x - 1, ZZ, ZZ)
Map with retraction
  from integers
  to integers

julia> g = inv(f)
Map with section
  from integers
  to integers

julia> h = f*g
Composite map
  from integers
  to integers
which is the composite of
  Map: integers -> integers
  Map: integers -> integers

julia> a = h(ZZ(1))
1

```
