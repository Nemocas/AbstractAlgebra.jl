# Extending the interface of AbstractAlgebra.jl

In this section we will discuss on how to extend the interface of
AbstractAlgebra.jl.

## Elements and parents

Any implementation with elements and parents should implement the following
interface:

```@docs
parent
elem_type
parent_type
```

### Aquiring associated elements and parents

Further, if one has a base ring, like polynomials over the integers
$\mathbb{Z}[x]$, then one should implement

```@docs
base_ring
```

## Special elements

For rings, one has to extend the following methods:

```@docs
one
zero
```

Groups should only extend at least one of these. The one that is required
depends on if the group is additive (commutative) or multiplicative.

## Basic manipulation

If one would like to implement a ring, these are the basic manipulation methods
that all rings should extend:

```@docs
isone
iszero
isunit
```

With the same logic as earlier, groups only need to extend one of the methods
`isone` and `iszero`.
