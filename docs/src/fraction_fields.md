# Fraction Field Interface

Fraction fields are supported in AbstractAlgebra.jl, at least for gcd domains.
In addition to the standard Ring interface, some additional functions are required to be
present for fraction fields.

## Types and parents

AbstractAlgebra provides two abstract types for fraction fields and their elements:

  * `FracField{T}` is the abstract type for fraction field parent types
  * `FracElem{T}` is the abstract type for types of fractions

We have that `FracField{T} <: AbstractAlgebra.Field` and 
`FracElem{T} <: AbstractAlgebra.FieldElem`.

Note that both abstract types are parameterised. The type `T` should usually be the type
of elements of the base ring of the fraction field.

Fraction fields should be made unique on the system by caching parent objects (unless
an optional `cache` parameter is set to `false`). Fraction fields should at least be
distinguished based on their base ring.

See `src/generic/GenericTypes.jl` for an example of how to implement such a cache (which
usually makes use of a dictionary).

## Required functionality for fraction fields

In addition to the required functionality for the Field interface the Fraction Field
interface has the following required functions.

We suppose that `R` is a fictitious base ring, and that `S` is the fraction field with 
parent object `S` of type `MyFracField{T}`. We also assume the fractions in the field 
have type `MyFrac{T}`, where `T` is the type of elements of the base ring.

Of course, in practice these types may not be parameterised, but we use parameterised
types here to make the interface clearer.

Note that the type `T` must (transitively) belong to the abstract type `RingElem`.

### Constructors

We provide the following constructors. Note that these constructors don't require
construction of the parent object first. This is easier to achieve if the fraction
element type doesn't contain a reference to the parent object, but merely contains a
reference to the base ring. The parent object can then be constructed on demand.

```julia
//(x::T, y::T) where T <: AbstractAlgebra.RingElem
```

Return the fraction $x/y$.

```julia
//(x::T, y::AbstractAlgebra.FracElem{T}) where T <: AbstractAlgebra.RingElem
```

Return $x/y$ where $x$ is in the base ring of $y$.

```julia
//(x::AbstractAlgebra.FracElem{T}, y::T) where T <: AbstractAlgebra.RingElem
```

Return $x/y$ where $y$ is in the base ring of $x$.

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")

f = (x^2 + x + 1)//(x^3 + 3x + 1)
g = f//x
h = x//f
```

### Basic manipulation of fields and elements

```julia
numerator(d::MyFrac{T}) where T <: AbstractAlgebra.RingElem
```

Given a fraction $d = a/b$ return $a$, where $a/b$ is in lowest terms with respect to
the `canonical_unit` and `gcd` functions on the base ring.

```julia
denominator(d::MyFrac{T}) where T <: AbstractAlgebra.RingElem
```

Given a fraction $d = a/b$ return $b$, where $a/b$ is in lowest terms with respect to
the `canonical_unit` and `gcd` functions on the base ring.

**Examples**

```julia
R, x = PolynomialRing(QQ, "x")

f = (x^2 + x + 1)//(x^3 + 3x + 1)

n = numerator(f)
d = denominator(f)
```

