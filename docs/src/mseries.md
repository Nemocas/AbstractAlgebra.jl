```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Multivariate series

AbstractAlgebra.jl provide multivariate series over a commutative ring.

Series with capped absolute precision are provided with and without
weights.

For the unweighted case precision in each variable can be set per series,
but is capped at some maximum precision which is set when defining the
ring.

For the weighted case, a single precision is set on the ring only. Terms
are truncated at that precision (after applying weights).

## Generic multivariate series

Generic multivariate series over a commutative ring,
`AbsMSeries{T}` is implemented in `src/generic/AbsMSeries.jl`.

Such series are capped absolute series and have type `Generic.AbsMSeries{T}`
where `T` is the type of elements of the coefficient ring.

Internally they consist of a multivariate polynomial. For unweighted series
they also contain a vector of precisions, one for each variable.

For weighted series weights and a precision are stored on the ring only.
The vector of precisions in the series objects is ignored.

See the file `src/generic/GenericTypes.jl` for details of the type.

The series are implemented in terms of multivariate polynomials which are
used internally to keep track of the coefficients of the series.

Only lex ordering is provided at present both weighted and unweighted,
though series print in reverse order to what multivariate polynomials would
print, i.e. least significant term first, as would be expected for series.

Parent objects of such series have type `Generic.AbsMSeriesRing{T}`.

The symbol representation of the variables and the multivariate polynomial
ring is stored in the parent object.

## Abstract types

Multivariate series element types belong to the abstract type
`MSeriesElem{T}` and the multivariate series ring types
belong to the abstract type `MSeriesRing{T}`. This enables
one to write generic functions that can accept any AbstractAlgebra
multivariate series type.

## Multivariate series ring constructors

In order to construct multivariate series in AbstractAlgebra.jl, one must first
construct the series ring itself. This is accomplished with the following
constructors.

For the unweighted case:

```julia
power_series_ring(R::Ring, prec::Vector{Int}, s::AbstractVector{<:VarName}; cached::Bool = true)
```

Given a base ring `R` and a vector of strings `s` specifying how the generators
(variables) should be printed, along with a vector of precisions, one for each
variable, return a tuple `U, (x, y, ...)` representing the new series ring $S$
and the generators $x, y, \ldots$ of the ring as a tuple. By default the parent
object `S` will depend on `R`, the precision vector and the variable names
`x, y, ...` and will be cached. Setting the optional argument `cached` to
`false` will prevent the parent object `S` from being cached.

In the weighted case:

```
power_series_ring(R::Ring, weights::Vector{Int}, s::AbstractVector{<:VarName}, prec::Int; cached::Bool = true)
```

Given a base ring `R` and a vector of strings `s` specifying how the generators
(variables) should be printed, along with a vector of weights, one for each
variable and a bound on the (weighted) precision, return a tuple
`U, (x, y, ...)` representing the new series ring $S$ and the generators
$x, y, \ldots$ of the ring as a tuple. By default the parent object `S` will
depend on `R`, the precision, the vector of weights and the variable names
`x, y, ...` and will be cached. Setting the optional argument `cached` to
`false` will prevent the parent object `S` from being cached.

Here are some examples of creating multivariate series rings and making use of
the resulting parent objects to coerce various elements into the series ring.

Note that one can also use the function call `O(x^n)` with unweighted series to
specify the precision in the variable `x` of a given series expression should be
precision `n`.

!!! note

    It is not possible to use `x^0` in the `O()` function, since there is no
    distinction between `x^0` and `y^0` as far as the system is concerned. If one
    wishes to set the precision of a variable to precision `0`, one must use the
    `set_precision!` function described below.

If one wants a series with the same precision in all variables, one can use
`O(R, n)` where `R` is the series ring and `n` is the desired precision.

If all the precisions are to be the same, the vector of integers for the
precisions can be replaced by a single integer in the constructor.

**Examples**

```jldoctest
julia> R, (x, y) = power_series_ring(ZZ, [2, 3], ["x", "y"])
(Multivariate power series ring in 2 variables over integers, AbstractAlgebra.Generic.AbsMSeries{BigInt, AbstractAlgebra.Generic.MPoly{BigInt}}[x + O(y^3) + O(x^2), y + O(y^3) + O(x^2)])

julia> f = R()
O(y^3) + O(x^2)

julia> g = R(123)
123 + O(y^3) + O(x^2)

julia> h = R(BigInt(1234))
1234 + O(y^3) + O(x^2)

julia> k = R(x + 1)
1 + x + O(y^3) + O(x^2)

julia> m = x + y + O(y^2)
y + x + O(y^2) + O(x^2)

julia> R, (x, y) = power_series_ring(ZZ, 3, ["x", "y"])
(Multivariate power series ring in 2 variables over integers, AbstractAlgebra.Generic.AbsMSeries{BigInt, AbstractAlgebra.Generic.MPoly{BigInt}}[x + O(y^3) + O(x^3), y + O(y^3) + O(x^3)])

julia> n = x + y + O(R, 2)
y + x + O(y^2) + O(x^2)

julia> R, (x, y) = power_series_ring(ZZ, [2, 3], 10, ["x", "y"])
(Multivariate power series ring in 2 variables over integers, AbstractAlgebra.Generic.AbsMSeries{BigInt, AbstractAlgebra.Generic.MPoly{BigInt}}[x + O(10), y + O(10)])

julia> R()
O(10)

julia> R(x)
x + O(10)
```

## Basic ring functionality

Once a multivariate series ring is constructed, there are various ways to
construct series in that ring.

The easiest way is simply using the generators returned by the
`power_series_ring` constructor and build up the power series using basic
arithmetic, as described in the Ring interface.

The power series rings in AbstractAlgebra.jl implement the full Ring interface.

We give some examples of such functionality. 

!!! note

    The divexact function can currently only divide by unit series (i.e. whose
    constant coefficient is invertible).

**Examples**

```jldoctest
julia> R, (x,) = power_series_ring(ZZ, [5], ["x"])
(Multivariate power series ring in 1 variable over integers, AbstractAlgebra.Generic.AbsMSeries{BigInt, AbstractAlgebra.Generic.MPoly{BigInt}}[x + O(x^5)])

julia> f = x^3 + 3x + 21
21 + 3*x + x^3 + O(x^5)

julia> h = zero(R)
O(x^5)

julia> k = one(R)
1 + O(x^5)

julia> isone(k)
true

julia> iszero(f)
false

julia> n = length(f)
3

julia> U = base_ring(R)
Integers

julia> v = symbols(R)
1-element Vector{Symbol}:
 :x

julia> T = parent(x + 1)
Multivariate power series ring in 1 variable x
  over integers

julia> f == deepcopy(f)
true

julia> t = divexact(f*x, 1 + x)
21*x - 18*x^2 + 18*x^3 - 17*x^4 + O(x^5)

julia> R, (x, y) = power_series_ring(ZZ, [2, 3], 10, ["x", "y"])
(Multivariate power series ring in 2 variables over integers, AbstractAlgebra.Generic.AbsMSeries{BigInt, AbstractAlgebra.Generic.MPoly{BigInt}}[x + O(10), y + O(10)])

julia> f = 3x^2*y + 1
1 + 3*y*x^2 + O(10)

julia> one(R)
1 + O(10)
```

## Power series functionality provided by AbstractAlgebra.jl

The functionality listed below is automatically provided by AbstractAlgebra.jl
for absolute series over any commutative ring.

### Basic functionality

The following are provided for weighted and unweighted series:

```@docs
number_of_variables(::Generic.AbsMSeriesRing)
```

```@docs
symbols(::MSeriesRing)
```

```@docs
precision(::Generic.AbsMSeries)
```

```@docs
coeff(::Generic.AbsMSeries, ::Int)
```

```@docs
characteristic(::Generic.AbsMSeries)
```

```@docs
gen(::Generic.AbsMSeriesRing, ::Int)
```

```@docs
gens(::Generic.AbsMSeriesRing)
```

```@docs
is_gen(::Generic.AbsMSeries)
```

```@docs
is_unit(::Generic.AbsMSeries)
```

```@docs
length(::Generic.AbsMSeries)
```

The following are only available for unweighted series.

```@docs
max_precision(::Generic.AbsMSeriesRing)
```

```@docs
valuation(::Generic.AbsMSeries)
```

### Iteration

```@docs
coefficients(::Generic.AbsMSeries)
```

```@docs
exponent_vectors(::Generic.AbsMSeries)
```

### Truncation

```@docs
truncate(::Generic.AbsMSeries, ::Vector{Int})
truncate(::Generic.AbsMSeries, ::Int)
```

### Exact division

```@docs
divexact(::Generic.AbsMSeries{T}, ::Generic.AbsMSeries{T}) where T <: RingElem
```

### Evaluation


```@docs
evaluate(::U, ::Vector{Int}, ::Vector{U}) where {T <: RingElement, U <: Generic.AbsMSeries{T}}
```

```@docs
evaluate(::U, ::Vector{U}, ::Vector{U}) where {T <: RingElement, U <: Generic.AbsMSeries{T}}
```

```@docs
evaluate(::U, ::Vector{U}) where {T <: RingElement, U <: Generic.AbsMSeries{T}}
```

### Random generation

```@docs
rand(::MSeriesRing, term_range, v...)
```


