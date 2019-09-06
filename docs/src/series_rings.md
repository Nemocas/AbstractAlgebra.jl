```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Series Ring Interface

Univariate power series rings are supported in AbstractAlgebra in a variety of different
forms, including absolute and relative precision models and Laurent series.

In addition to the standard Ring interface, numerous additional functions are required
to be present for power series rings.

## Types and parents

AbstractAlgebra provides two abstract types for power series rings and their elements:

  * `SeriesRing{T}` is the abstract type for all power series ring parent types
  * `SeriesElem{T}` is the abstract type for all power series types

We have that `SeriesRing{T} <: AbstractAlgebra.Ring` and
`SeriesElem{T} <: AbstractAlgebra.RingElem`.

Note that both abstract types are parameterised. The type `T` should usually be the type
of elements of the coefficient ring of the power series ring. For example, in the case
of $\mathbb{Z}[[x]]$ the type `T` would be the type of an integer, e.g. `BigInt`.

Within the `SeriesElem{T}` abstract type is the abstract type `RelSeriesElem{T}` for
relative power series, and `AbsSeriesElem{T}` for absolute power series.

Relative series are typically stored with a valuation and a series that is either
zero or that has nonzero constant term. Absolute series are stored starting from the
constant term, even if it is zero.

If the parent object for a relative series ring over the bignum integers has type
`MySeriesRing` and series in that ring have type `MySeries` then one would have:

  * `MySeriesRing <: SeriesRing{BigInt}`
  * `MySeries <: RelSeriesElem{BigInt}`

Series rings should be made unique on the system by caching parent objects (unless
an optional `cache` parameter is set to `false`). Series rings should at least be
distinguished based on their base (coefficient) ring. But if they have the same base
ring and symbol (for their variable/generator) and same default precision, they should
certainly have the same parent object.

See `src/generic/GenericTypes.jl` for an example of how to implement such a cache (which
usually makes use of a dictionary).

## Required functionality for series

In addition to the required functionality for the Ring interface the Series Ring
interface has the following required functions.

We suppose that `R` is a fictitious base ring (coefficient ring) and that `S` is a
series ring over `R` (e.g. $S = R[[x]]$) with parent object `S` of type
`MySeriesRing{T}`. We also assume the series in the ring have type `MySeries{T}`, where
`T` is the type of elements of the base (coefficient) ring.

Of course, in practice these types may not be parameterised, but we use parameterised
types here to make the interface clearer.

Note that the type `T` must (transitively) belong to the abstract type `RingElem`.

### Constructors

In addition to the standard constructors, the following constructors, taking an array of
coefficients, must be available.

For relative power series and Laurent series we have:

```julia
(S::MySeriesRing{T})(A::Array{T, 1}, len::Int, prec::Int, val::Int) where T <: AbstractAlgebra.RingElem
```

Create the series in the given ring whose valuation is `val`, whose absolute precision
is given by `prec` and the coefficients of which are given by `A`, starting from the
first nonzero term. Only `len` terms of the array are used, the remaining terms being
ignored. The value `len` cannot exceed the length of the supplied array.

It is permitted to have trailing zeros in the array, but it is not needed, even if the
precision minus the valuation is bigger than the length of the array.

**Examples**

```jldoctest
julia> S, x = PowerSeriesRing(QQ, 10, "x"; model=:capped_relative)
(Univariate power series ring in x over Rationals, x+O(x^11))

julia> T, y = LaurentSeriesRing(ZZ, 10, "y")
(Laurent series ring in y over Integers, y+O(y^11))

julia> U, z = LaurentSeriesField(QQ, 10, "z")
(Laurent series field in z over Rationals, z+O(z^11))
```

For absolute power series we have:

```julia
(S::MySeriesRing{T})(A::Array{T, 1}, len::Int, prec::Int) where T <: AbstractAlgebra.RingElem
```

Create the series in the given ring whose absolute precision is given by `prec` and the
coefficients of which are given by `A`, starting from the constant term. Only `len`
terms of the array are used, the remaining terms being ignored.

Note that `len` is usually maintained separately of any polynomial that is underlying
the power series. This allows for easy trucation of a power series without actually
modifying the polynomial underlying it.

It is permitted to have trailing zeros in the array, but it is not needed, even if the
precision is bigger than the length of the array.

**Examples**

```jldoctest
julia> S, x = PowerSeriesRing(QQ, 10, "x"; model=:capped_absolute)
(Univariate power series ring in x over Rationals, x+O(x^10))

julia> f = S(Rational{BigInt}[0, 2, 3, 1], 4, 6)
(2//1)*x+(3//1)*x^2+x^3+O(x^6)

```

### Data type and parent object methods

```julia
var(S::MySeriesRing{T}) where T <: AbstractAlgebra.RingElem
```

Return a `Symbol` representing the variable (generator) of the series ring. Note
that this is a `Symbol` not a `String`, though its string value will usually be used
when printing series.

```julia
max_precision(S::MySeriesRing{T}) where T <: AbstractAlgebra.RingElem
```

Return the (default) maximum precision of the power series ring. This is the precision
that the output of an operation will be if it cannot be represented to full precision
(e.g. because it mathematically has infinite precision).

This value is usually supplied upon creation of the series ring and stored in the ring.
It is independent of the precision which each series in the ring actually has. Those
are stored on a per element basis in the actual series elements.

**Examples**

```jldoctest
julia> S, x = PowerSeriesRing(QQ, 10, "x")
(Univariate power series ring in x over Rationals, x+O(x^11))

julia> vsym = var(S)
:x

julia> max_precision(S) == 10
true

```

### Basic manipulation of rings and elements

```julia
pol_length(f::MySeries{T}) where T <: AbstractAlgebra.RingElem
```

Return the length of the polynomial underlying the given power series. This is not
generally useful to the user, but is used internally.

```julia
set_length!(f::MySeries{T}, n::Int) where T <: AbstractAlgebra.RingElem
```

This function sets the effective length of the polynomial underlying the given series.
The function doesn't modify the actual polynomial, but simply changes the number of
terms of the polynomial which are considered to belong to the power series. The
remaining terms are ignored.

This function cannot set the length to a value greater than the length of any underlying
polynomial.

The function mutates the series in-place but does not return the mutated series.

```julia
precision(f::MySeries{T})
```

Return the absolute precision of $f$.

```julia
set_prec!(f::MySeries{T}, prec::Int)
```

Set the absolute precision of the given series to the given value.

This function mutates the series in-place but does not return the mutated series.

```julia
valuation(f::MySeries{T})
```

Return the valuation of the given series.

```julia
set_val!(f::MySeries{T}, val::Int)
```

For relative series and Laurent series only, this function alters the valuation of the
given series to the given value.

The series is mutated in-place but does not return the mutated series.

```julia
polcoeff(f::MySeries{T}, n::Int)
```

Return the coefficient of degree `n` of the polynomial underlying the series. If `n`
is larger than the degree of this polynomial, zero is returned. This function is not
generally of use to the user but is used internally.


```julia
setcoeff!(f::MySeries{T}, n::Int, a::T) where T <: AbstractAlgebra.RingElem
```

Set the degree $n$ coefficient of the polynomial underlying $f$ to $a$. This mutates
the polynomial in-place if possible and returns the mutated series (so that immutable
types can also be supported). The function must not assume that the polynomial already
has space for $n + 1$ coefficients. The polynomial must be resized if this is not the
case.

Note that this function is not required to normalise the polynomial and is not
necessarily useful to the user, but is used extensively by the generic functionality in
AbstractAlgebra.jl. It is for setting raw coefficients in the representation.

```julia
normalise(f::MySeries{T}, n::Int)
```

Given a series $f$ represented by a polynomial of at least the given length, return the
normalised length of the underlying polynomial assuming it has length at most $n$. This
function does not actually normalise the polynomial and is not particularly useful to
the user. It is used internally.

```julia
renormalize!(f::MySeries{T}) where T <: AbstractAlgebra.RingElem
```

Given a relative series or Laurent series whose underlying polynomial has zero constant
term, say as the result of some internal computation, renormalise the series so that the
polynomial has nonzero constant term. The precision and valuation of the series are
adjusted to compensate. This function is not intended to be useful to the user, but is
used internally.

```julia
fit!(f::MySeries{T}, n::Int) where T <: AbstractAlgebra.RingElem
```

Ensure that the polynomial underlying $f$ internally has space for $n$ coefficients.
This function must mutate the series in-place if it is mutable. It does not return the
mutated series. Immutable types can still be supported by defining this function to do
nothing.

Some interfaces for C polynomial types automatically manage the internal allocation of
polynomials in every function that can be called on them. Explicit adjustment by
the generic code in AbstractAlgebra.jl is not required. In such cases, this function
can also be defined to do nothing.

```julia
gen(R::MySeriesRing{T}) where T <: AbstractAlgebra.RingElem
```

Return the generator `x` of the series ring.

**Examples**

```jldoctest
julia> S, x = PowerSeriesRing(ZZ, 10, "x")
(Univariate power series ring in x over Integers, x+O(x^11))

julia> f = 1 + 3x + x^3 + O(x^5)
1+3*x+x^3+O(x^5)

julia> g = S(BigInt[1, 2, 0, 1, 0, 0, 0], 4, 10, 3);

julia> n = pol_length(f)
4

julia> c = polcoeff(f, 1)
3

julia> set_length!(g, 3)
3

julia> g = setcoeff!(g, 2, BigInt(11))
x^3+2*x^4+11*x^5+O(x^10)

julia> fit!(g, 8)

julia> g = setcoeff!(g, 7, BigInt(4))
x^3+2*x^4+11*x^5+O(x^10)

julia> w = gen(S)
x+O(x^11)

julia> isgen(w)
true

```
