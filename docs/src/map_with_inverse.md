# Map with inverse

It is not possible to provide generic functionality to invert a map. However, sometimes
one knows an inverse map explicitly and would like to keep track of this.

Recall that as map composition is not commutative, there is a notion of a left inverse
and a right inverse for maps.

To keep track of such inverse maps, AbstractAlgebra provides data types
`Generic.MapWithRetraction` and `GenericMapWithSection`.

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

```julia
map_with_retraction(m::Map{D, C}, r::Map{C, D}) where {D, C}
map_with_section(m::Map{D, C}, s::Map{C, D}) where {D, C}
```

Construct the map with retraction/section given a known retraction/section $r$ or $s$
respectively, of $m$.

For convenience we allow construction of maps with retraction/section from a pair of
Julia functions/closures.

```julia
map_with_retraction_from_func(R, S, f::Function, r::Function)
map_with_section_from_func(R, S, f::Function, s::Function)
```

Construct the map with retraction/section such that the map is given by the function $f$
and the retraction/section is given by the function $r$ or $s$ respectively. Here $R$ is
the parent object representing the domain and $S$ is the parent object representing the
codomain of $f$.

**Examples**

```julia
f = map_with_retraction_from_func(ZZ, ZZ, x -> x + 1, x -> x - 1)

a = f(ZZ(1))
```

## Functionality for maps with inverses

The following functionality is provided for maps with inverses.

```julia
inv(M::MapWithRetraction)
inv(M::MapWithSection)
```

Return the map with the two maps contained in $M$ swapped. In the first case, a
`MapWithSection` is returned. In the second case a `MapWithRetraction` is returned.

To access the two maps stored in a map with retraction/section, we have the following:

```julia
image_map(M::MapWithRetraction)
image_map(M::MapWithSection)
retraction_map(M::MapWithRetraction)
section_map(M::MapWithSection)
```

The first two of these functions return the first map in a map with retraction/section,
the second two functions return the corresponding second maps.

**Examples**

```julia
f = map_with_retraction_from_func(ZZ, ZZ, x -> x + 1, x -> x - 1)
g = inv(f)
h = f*g

a = h(ZZ(1))
```

