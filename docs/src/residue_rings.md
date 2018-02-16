# Residue Ring Interface

Residue rings (currently a quotient ring modulo a principal ideal) are supported in
AbstractAlgebra.jl, at least for Euclidean base rings. In addition to the standard Ring
interface, some additional functions are required to be present for residue rings.

## Types and parents

AbstractAlgebra provides four abstract types for residue rings and their elements:

  * `ResRing{T}` is the abstract type for residue ring parent types
  * `ResField{T}` is the abstract type for residue rings known to be fields
  * `ResElem{T}` is the abstract type for types of elements of residue rings (residues)
  * `ResFieldElem{T}` is the abstract type for types of elements of residue fields

We have that `ResRing{T} <: AbstractAlgebra.Ring` and 
`ResElem{T} <: AbstractAlgebra.RingElem`.

Note that these abstract types are parameterised. The type `T` should usually be the type
of elements of the base ring of the residue ring/field.

If the parent object for a residue ring has type `MyResRing` and residues in that ring
have type `MyRes` then one would have:

  * `MyResRing <: ResRing{BigInt}`
  * `MyRes <: ResElem{BigInt}`

Residue rings should be made unique on the system by caching parent objects (unless
an optional `cache` parameter is set to `false`). Residue rings should at least be
distinguished based on their base ring and modulus (the principal ideal one is taking
a quotient of the base ring by).

See `src/generic/GenericTypes.jl` for an example of how to implement such a cache (which
usually makes use of a dictionary).

## Required functionality for residue rings

In addition to the required functionality for the Ring interface the Residue Ring
interface has the following required functions.

We suppose that `R` is a fictitious base ring, $m$ is an element of that ring, and that
`S` is the residue ring (quotient ring) $R/(m)$ with parent object `S` of type
`MyResRing{T}`. We also assume the residues $r \pmod{m}$ in the residue ring have type
`MyRes{T}`, where `T` is the type of elements of the base ring.

Of course, in practice these types may not be parameterised, but we use parameterised
types here to make the interface clearer.

Note that the type `T` must (transitively) belong to the abstract type `RingElem`.

### Data type and parent object methods

```julia
modulus(S::MyResRing{T}) where T <: AbstractAlgebra.RingElem
```

Return the modulus of the given residue ring, i.e. if the residue ring $S$ was specified
to be $R/(m)$, return $m$.

**Examples**

```julia
R, x = PolynomialRing(QQ, "x")
S = ResidueRing(R, x^3 + 3x + 1)

m = modulus(S)
```

### Basic manipulation of rings and elements

```julia
data(f::MyRes{T}) where T <: AbstractAlgebra.RingElem
```

Given a residue $r \pmod{m}$, represented as such, return $r$.

**Examples**

```julia
R, x = PolynomialRing(QQ, "x")
S = ResidueRing(R, x^3 + 3x + 1)

f = S(x^2 + 2)

d = data(f)
```

