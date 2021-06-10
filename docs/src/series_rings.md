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

We have that `SeriesRing{T} <: Ring` and
`SeriesElem{T} <: RingElem`.

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
(S::MySeriesRing{T})(A::Array{T, 1}, len::Int, prec::Int, val::Int) where T <: RingElem
```

Create the series in the given ring whose valuation is `val`, whose absolute precision
is given by `prec` and the coefficients of which are given by `A`, starting from the
first nonzero term. Only `len` terms of the array are used, the remaining terms being
ignored. The value `len` cannot exceed the length of the supplied array.

It is permitted to have trailing zeros in the array, but it is not needed, even if the
precision minus the valuation is bigger than the length of the array.

```julia
(S::MySeriesRing{T})(A::Array{U, 1}, len::Int, prec::Int, val::Int) where {T <: RingElem, U <: RingElem}
```

As above, but where the array is an array of coefficient that can be coerced into the
base ring of the series ring.

```julia
(S::MySeriesRing{T})(A::Array{U, 1}, len::Int, prec::Int, val::Int) where {T <: RingElem, U <: Integer}
```

As above, but where the array is an array of integers that can be coerced into the
base ring of the series ring.

It may be desirable to implement an addition version which accepts an array of Julia
`Int` values if this can be done more efficiently.

**Examples**

```jldoctest
julia> S, x = PowerSeriesRing(QQ, 10, "x"; model=:capped_relative)
(Univariate power series ring in x over Rationals, x + O(x^11))

julia> T, y = LaurentSeriesRing(ZZ, 10, "y")
(Laurent series ring in y over Integers, y + O(y^11))

julia> U, z = LaurentSeriesField(QQ, 10, "z")
(Laurent series field in z over Rationals, z + O(z^11))
```

For absolute power series we have:

```julia
(S::MySeriesRing{T})(A::Array{T, 1}, len::Int, prec::Int) where T <: RingElem
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
(Univariate power series ring in x over Rationals, x + O(x^10))

julia> f = S(Rational{BigInt}[0, 2, 3, 1], 4, 6)
2*x + 3*x^2 + x^3 + O(x^6)

```

It is also possible to create series directly without having to create the
corresponding series ring.

```julia
abs_series(R::Ring, arr::Vector{T}, len::Int, prec::Int, var::AbstractString="x"; max_precision::Int=prec, cached::Bool=true) where T
rel_series(R::Ring, arr::Vector{T}, len::Int, prec::Int, val::Int, var::AbstractString="x"; max_precision::Int=prec, cached::Bool=true) where T
```

Create the power series over the given base ring `R` with coefficients
specified by `arr` with the given absolute precision `prec` and in the case
of relative series with the given valuation `val`.

Note that more coefficients may be specified than are actually used. Only
the first `len` coefficients are made part of the series, the remainder being
stored internally but ignored.

In the case of absolute series one must have `prec >= len` and in the case of
relative series one must have `prec >= len + val`.

By default the series are created in a ring with variable `x` and
`max_precision` equal to `prec`, however one may specify these directly to
override the defaults. Note that series are only compatible if they have the
same coefficient ring `R`, `max_precision` and variable string `var`.

Also by default any parent ring created is cached. If this behaviour is not
desired, set `cached=false`. However, this means that subsequent series created
in the same way will not be compatible. Instead, one should use the parent
object of the first series to create subsequent series instead of calling this
function repeatedly with cached=false.

**Examples**

```jldoctest
julia> f = abs_series(ZZ, [1, 2, 3], 3, 5, "y")
1 + 2*y + 3*y^2 + O(y^5)

julia> g = rel_series(ZZ, [1, 2, 3], 3, 7, 4)
x^4 + 2*x^5 + 3*x^6 + O(x^7)

julia> k = abs_series(ZZ, [1, 2, 3], 1, 6, cached=false)
1 + O(x^6)

julia> p = rel_series(ZZ, BigInt[], 0, 3, 1)
O(x^3)

julia> q = abs_series(ZZ, [], 0, 6)
O(x^6)

julia> s = abs_series(ZZ, [1, 2, 3], 3, 5; max_precision=10)
1 + 2*x + 3*x^2 + O(x^5)
```

### Data type and parent object methods

```julia
var(S::MySeriesRing{T}) where T <: RingElem
```

Return a `Symbol` representing the variable (generator) of the series ring. Note
that this is a `Symbol` not a `String`, though its string value will usually be used
when printing series.

Custom series types over a given ring should define one of the following functions
which return the type of an absolute or relative series object over that ring.

```julia
abs_series_type(::Type{T}) where T <: RingElement
rel_series_type(::Type{T}) where T <: RingElement
```

Return the type of a series whose coefficients have the given type.

This function is defined for generic series and only needs to be defined for
custom series rings, e.g. ones defined by a C implementation.

```julia
max_precision(S::MySeriesRing{T}) where T <: RingElem
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
(Univariate power series ring in x over Rationals, x + O(x^11))

julia> vsym = var(S)
:x

julia> max_precision(S) == 10
true

```

### Basic manipulation of rings and elements

```julia
pol_length(f::MySeries{T}) where T <: RingElem
```

Return the length of the polynomial underlying the given power series. This is not
generally useful to the user, but is used internally.

```julia
set_length!(f::MySeries{T}, n::Int) where T <: RingElem
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
set_precision!(f::MySeries{T}, prec::Int)
```

Set the absolute precision of the given series to the given value.

This return the updated series.

```julia
valuation(f::MySeries{T})
```

Return the valuation of the given series.

```julia
set_valuation!(f::MySeries{T}, val::Int)
```

For relative series and Laurent series only, this function alters the valuation of the
given series to the given value.

This function returns the updated series.

```julia
polcoeff(f::MySeries{T}, n::Int)
```

Return the coefficient of degree `n` of the polynomial underlying the series. If `n`
is larger than the degree of this polynomial, zero is returned. This function is not
generally of use to the user but is used internally.


```julia
setcoeff!(f::MySeries{T}, n::Int, a::T) where T <: RingElem
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
renormalize!(f::MySeries{T}) where T <: RingElem
```

Given a relative series or Laurent series whose underlying polynomial has zero constant
term, say as the result of some internal computation, renormalise the series so that the
polynomial has nonzero constant term. The precision and valuation of the series are
adjusted to compensate. This function is not intended to be useful to the user, but is
used internally.

```julia
fit!(f::MySeries{T}, n::Int) where T <: RingElem
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
gen(R::MySeriesRing{T}) where T <: RingElem
```

Return the generator `x` of the series ring.

**Examples**

```jldoctest
julia> S, x = PowerSeriesRing(ZZ, 10, "x")
(Univariate power series ring in x over Integers, x + O(x^11))

julia> f = 1 + 3x + x^3 + O(x^5)
1 + 3*x + x^3 + O(x^5)

julia> g = S(BigInt[1, 2, 0, 1, 0, 0, 0], 4, 10, 3);

julia> n = pol_length(f)
4

julia> c = polcoeff(f, 1)
3

julia> set_length!(g, 3)
x^3 + 2*x^4 + O(x^10)

julia> g = setcoeff!(g, 2, BigInt(11))
x^3 + 2*x^4 + 11*x^5 + O(x^10)

julia> fit!(g, 8)

julia> g = setcoeff!(g, 7, BigInt(4))
x^3 + 2*x^4 + 11*x^5 + O(x^10)

julia> w = gen(S)
x + O(x^11)

julia> isgen(w)
true

```

## Optional functionality for series

### Similar and zero

The following functions are available for all absolute and relative series
types. The functions `similar` and `zero` do the same thing, but are provided
for uniformity with other parts of the interface.

```julia
similar(x::MySeries, R::Ring, max_prec::Int, var::Symbol=var(parent(x)); cached::Bool=true)
zero(a::MySeries, R::Ring, max_prec::Int, var::Symbol=var(parent(a)); cached::Bool=true)
```

Construct the zero series with the given variable (if specified), coefficients
in the specified coefficient ring and with relative/absolute precision cap on
its parent ring as given by `max_prec`.

```julia
similar(x::MySeries, R::Ring, var::Symbol=var(parent(x)); cached::Bool=true)
similar(x::MySeries, max_prec::Int, var::Symbol=var(parent(x)); cached::Bool=true)
similar(x::MySeries, var::Symbol=var(parent(x)); cached::Bool=true)
similar(x::MySeries, R::Ring, max_prec::Int, var::String; cached::Bool=true)
similar(x::MySeries, R::Ring, var::String; cached::Bool=true)
similar(x::MySeries, max_prec::Int, var::String; cached::Bool=true)
similar(x::MySeries, var::String; cached::Bool=true)
zero(x::MySeries, R::Ring, var::Symbol=var(parent(x)); cached::Bool=true)
zero(x::MySeries, max_prec::Int, var::Symbol=var(parent(x)); cached::Bool=true)
zero(x::MySeries, var::Symbol=var(parent(x)); cached::Bool=true)
zero(x::MySeries, R::Ring, max_prec::Int, var::String; cached::Bool=true)
zero(x::MySeries, R::Ring, var::String; cached::Bool=true)
zero(x::MySeries, max_prec::Int, var::String; cached::Bool=true)
zero(x::MySeries, var::String; cached::Bool=true)
```

As above, but use the precision cap of the parent ring of `x` and the
`base_ring` of `x` if these are not specified.

Custom series rings may choose which series type is best-suited to return for
the given coefficient ring, precision cap and variable, however they should
return a series with the same model as `x`, i.e. relative or series.

If custom implementations don't specialise these function the default return
type is a `Generic.AbsSeries` or `Generic.RelSeries`.

The default implementation of zero calls out to similar, so it's generally
sufficient to specialise only similar. For both similar and zero only the most
general method has to be implemented as all other methods call out to this more
general method.

