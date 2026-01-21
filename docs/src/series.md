```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Power series

AbstractAlgebra.jl allows the creation of capped relative and absolute power series over
any computable commutative ring $R$.

Capped relative power series are power series of the form
$a_jx^j + a_{j+1}x^{j+1} + \cdots + a_{k-1}x^{k-1} + O(x^k)$
where $a_j \in R$ and the relative precision $k - j$ is at most
equal to some specified precision $n$.

Capped absolute power series are power series of the form
$a_jx^j + a_{j+1}x^{j+1} + \cdots + a_{n-1}x^{n-1} + O(x^n)$
where $j \geq 0$, $a_j \in R$ and the precision $n$ is fixed.

There are two implementations of relative series: relative power series, implemented in
`src/RelSeries.jl` for which $j > 0$ in the above description, and Laurent
series where $j$ can be negative, implemented in `src/Laurent.jl`.
Note that there are two implementations for Laurent series, one over rings and one
over fields, though in practice most of the implementation uses the same code in both
cases.

There is a single implementation of absolute series: absolute power series, implemented
in `src/AbsSeries.jl`.

## Generic power series types

AbstractAlgebra.jl provides generic series types implemented in
`src/generic/AbsSeries.jl`, `src/generic/RelSeries.jl` and
`src/generic/LaurentSeries.jl` which implement the Series interface.

These generic series have types `Generic.RelSeries{T}`, `Generic.AbsSeries{T}`,
`Generic.LaurentSeriesRingElem{T}` and `Generic.LaurentSeriesFieldElem{T}`. See
the file `src/generic/GenericTypes.jl` for details.

The parent objects have types `Generic.AbsPowerSeriesRing{T}`
and `Generic.RelPowerSeriesRing{T}` and `Generic.LaurentSeriesRing{T}` respectively.

The default precision, string representation of the variable and base ring $R$
of a generic power series are stored in its parent object.

## Abstract types
Relative power series elements belong to the abstract type `RelPowerSeriesRingElem`.

Laurent series elements belong directly to either `RingElem` or
`FieldElem` since it is more useful to be able to distinguish whether
they belong to a ring or field than it is to distinguish that they are relative
series.

Absolute power series elements belong to `AbsPowerSeriesRingElem`.

The parent types for relative and absolute power series,
`Generic.RelPowerSeriesRing{T}` and `Generic.AbsPowerSeriesRing{T}` respectively,
belong to `SeriesRing{T}`.

The parent types of Laurent series belong directly to `Ring` and `Field`
respectively.

## Series ring constructors

In order to construct series in AbstractAlgebra.jl, one must first construct the ring
itself. This is accomplished with any of the following constructors.

```julia
power_series_ring(R::Ring, prec_max::Int, s::VarName; cached::Bool = true, model::Symbol=:capped_relative)
```

```julia
laurent_series_ring(R::Ring, prec_max::Int, s::VarName; cached::Bool = true)
```

```julia
laurent_series_ring(R::Field, prec_max::Int, s::VarName; cached::Bool = true)
```

Given a base ring `R`, a maximum precision (relative or absolute, depending on the
model) and a string `s` specifying how the generator (variable) should be
printed, return a tuple `S, x` representing the series ring and its generator.

By default, `S` will depend only on `S`, `x` and the maximum precision and will be
cached. Setting the optional argument `cached` to `false` will prevent this.

In the case of power series, the optional argument `model` can be set to either
`:capped_absolute` or `:capped_relative`, depending on which power series model is
required.

It is also possible to construct absolute and relative power series with a
default variable. These are lightweight constructors and should be used in
generic algorithms wherever possible when creating series rings where the
symbol does not matter.

```julia
AbsPowerSeriesRing(R::Ring, prec::Int)
RelPowerSeriesRing(R::Ring, prec::Int)
```

Return the absolute or relative power series ring over the given base ring $R$
and with precision cap given by `prec`. Note that a tuple is not returned, only
the power series ring itself, not a generator.

Here are some examples of constructing various kinds of series rings and coercing
various elements into those rings.

**Examples**

```jldoctest
julia> R, x = power_series_ring(ZZ, 10, :x)
(Univariate power series ring over integers, x + O(x^11))

julia> S, y = power_series_ring(ZZ, 10, :y; model=:capped_absolute)
(Univariate power series ring over integers, y + O(y^10))

julia> T, z = laurent_series_ring(ZZ, 10, :z)
(Laurent series ring in z over integers, z + O(z^11))

julia> U, w = laurent_series_field(QQ, 10, :w)
(Laurent series field in w over rationals, w + O(w^11))

julia> f = R()
O(x^10)

julia> g = S(123)
123 + O(y^10)

julia> h = U(BigInt(1234))
1234 + O(w^10)

julia> k = T(z + 1)
1 + z + O(z^10)

julia> V = AbsPowerSeriesRing(ZZ, 10)
Univariate power series ring in x with precision 10
  over integers
```

## Power series constructors

Series can be constructed using arithmetic operators using the generator of the
series. Also see the big-oh notation below for specifying the precision.

All of the standard ring constructors can also be used to construct power series.

```julia
(R::SeriesRing)() # constructs zero
(R::SeriesRing)(c::Integer)
(R::SeriesRing)(c::elem_type(R))
(R::SeriesRing{T})(a::T) where T <: RingElement
```

In addition, the following constructors that are specific to power series are
provided. They take an array of coefficients, a length, precision and
valuation. Coefficients will be coerced into the coefficient ring if they are
not already in that ring.

For relative series we have:

```julia
(S::SeriesRing{T})(A::Vector{T}, len::Int, prec::Int, val::Int) where T <: RingElem
(S::SeriesRing{T})(A::Vector{U}, len::Int, prec::Int, val::Int) where {T <: RingElem, U <: RingElem}
(S::SeriesRing{T})(A::Vector{U}, len::Int, prec::Int, val::Int) where {T <: RingElem, U <: Integer}
```

And for absolute series:

```julia
(S::SeriesRing{T})(A::Vector{T}, len::Int, prec::Int) where T <: RingElem
```

It is also possible to create series directly without having to create the
corresponding series ring.

```julia
abs_series(R::Ring, arr::Vector{T}, len::Int, prec::Int, var::VarName=:x; max_precision::Int=prec, cached::Bool=true) where T
rel_series(R::Ring, arr::Vector{T}, len::Int, prec::Int, val::Int, var::VarName=:x; max_precision::Int=prec, cached::Bool=true) where T
laurent_series(R::Ring, arr::Vector{T}, len::Int, prec::Int, val::Int, scale::Int, var::VarName=:x; max_precision::Int=prec, cached::Bool=true) where T
```

**Examples**

```jldoctest
julia> S, x = power_series_ring(QQ, 10, :x; model=:capped_absolute)
(Univariate power series ring over rationals, x + O(x^10))

julia> f = S(Rational{BigInt}[0, 2, 3, 1], 4, 6)
2*x + 3*x^2 + x^3 + O(x^6)

julia> f = abs_series(ZZ, [1, 2, 3], 3, 5, :y)
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

julia> s = laurent_series(ZZ, [1, 2, 3], 3, 5, 0, 2; max_precision=10)
1 + 2*x^2 + 3*x^4 + O(x^5)
```

## Big-oh notation

Series elements can be given a precision using the big-oh notation. This is provided
by a function of the following form, (or something equivalent for Laurent series):

```julia
O(x::SeriesElem)
```

**Examples**

```jldoctest
julia> R, x = power_series_ring(ZZ, 10, :x)
(Univariate power series ring over integers, x + O(x^11))

julia> S, y = laurent_series_ring(ZZ, 10, :y)
(Laurent series ring in y over integers, y + O(y^11))

julia> f = 1 + 2x + O(x^5)
1 + 2*x + O(x^5)

julia> g = 2y + 7y^2 + O(y^7)
2*y + 7*y^2 + O(y^7)
```

What is happening here in practice is that `O(x^n)` is creating the series `0 + O(x^n)`
and the rules for addition of series dictate that if this is added to a series of
greater precision, then the lower of the two precisions must be used.

Of course it may be that the precision of the series that `O(x^n)` is added to is
already lower than `n`, in which case adding `O(x^n)` has no effect. This is the case
if the default precision is too low, since `x` on its own has the default precision.

## Power series models

Capped relative power series have their maximum relative precision capped at
some value `prec_max`. This means that if the leading term of a nonzero
power series element is $c_ax^a$ and the precision is $b$ then the power series
is of the form  $c_ax^a + c_{a+1}x^{a+1} + \ldots + O(x^{a + b})$.

The zero power series is simply taken to be $0 + O(x^b)$.

The capped relative model has the advantage that power series are stable
multiplicatively. In other words, for nonzero power series $f$ and $g$ we
have that `divexact(f*g), g) == f`.

However, capped relative power series are not additively stable, i.e. we
do not always have $(f + g) - g = f$.

Similar comments apply to Laurent series.

On the other hand, capped absolute power series have their absolute precision capped.
This means that if the leading term of a nonzero power series element is $c_ax^a$ and
the precision is $b$ then the power series is of the form
$c_ax^a + c_{a+1}x^{a+1} + \ldots + O(x^b)$.

Capped absolute series are additively stable, but not necessarily multiplicatively
stable.

For all models, the maximum precision is also used as a default precision in the case
of coercing coefficients into the ring and for any computation where the result could
mathematically be given to infinite precision.

In all models we say that two power series are equal if they agree up to the minimum
**absolute** precision of the two power series.

Thus, for example, $x^5 + O(x^{10}) == 0 + O(x^5)$, since the minimum absolute
precision is $5$.

During computations, it is possible for power series to lose relative
precision due to cancellation. For example if $f = x^3 + x^5 + O(x^8)$ and
$g = x^3 + x^6 + O(x^8)$ then $f - g = x^5 - x^6 + O(x^8)$ which now has
relative precision $3$ instead of relative precision $5$.

Amongst other things, this means that equality is not transitive. For example
$x^6 + O(x^{11}) == 0 + O(x^5)$ and $x^7 + O(x^{12}) == 0 + O(x^5)$ but
$x^6 + O(x^{11}) \neq x^7 + O(x^{12})$.

Sometimes it is necessary to compare power series not just for arithmetic
equality, as above, but to see if they have precisely the same precision and
terms. For this purpose we introduce the `isequal` function.

For example, if $f = x^2 + O(x^7)$ and $g = x^2 + O(x^8)$ and $h = 0 + O(x^2)$
then $f == g$, $f == h$ and $g == h$, but `isequal(f, g)`, `isequal(f, h)` and
`isequal(g, h)` would all return `false`. However, if $k = x^2 + O(x^7)$ then
`isequal(f, k)` would return `true`.

There are further difficulties if we construct polynomial over power series.
For example, consider the polynomial in $y$ over the power series ring in $x$
over the rationals. Normalisation of such polynomials is problematic. For
instance, what is the leading coefficient of $(0 + O(x^{10}))y + (1 + O(x^{10}))$?

If one takes it to be $(0 + O(x^{10}))$ then some functions may not terminate
due to the fact that algorithms may require the degree of polynomials to
decrease with each iteration. Instead, the degree may remain constant and
simply accumulate leading terms which are arithmetically zero but not
identically zero.

On the other hand, when constructing power series over other power series, if
we simply throw away terms which are arithmetically equal to zero, our
computations may have different output depending on the order in which the
power series are added!

One should be aware of these difficulties when working with power series.
Power series, as represented on a computer, simply don't satisfy the axioms
of a ring. They must be used with care in order to approximate operations in
a mathematical power series ring.

Simply increasing the precision will not necessarily give a "more correct"
answer and some computations may not even terminate due to the presence of
arithmetic zeroes!

An absolute power series ring over a ring $R$ with precision $p$ behaves
very much like the quotient $R[x]/(x^p)$ of the polynomial ring over $R$. Therefore
one can often treat absolute power series rings as though they were rings. However,
this depends on all series being given a precision equal to the specified maximum
precision and not a lower precision.

## Functions for types and parents of series rings

```julia
base_ring(R::SeriesRing)
base_ring(a::SeriesElem)
```

Return the coefficient ring of the given series ring or series.


```julia
parent(a::SeriesElem)
```

Return the parent of the given series.


## Series functions

Unless otherwise noted, the functions below are available for all series
models, including Laurent series. We denote this by using the abstract type
`RelPowerSeriesRingElem`, even though absolute series and Laurent series types
do not belong to this abstract type.

### Basic functionality

Series implement the Ring Interface

```julia
zero(R::SeriesRing)
one(R::SeriesRing)
iszero(a::SeriesElem)
isone(a::SeriesElem)
```

```julia
divexact(a::T, b::T) where T <: SeriesElem
inv(a::SeriesElem)
```

Series also implement the Series Interface, the most important basic
functions being the following.

```julia
var(S::SeriesRing)
```

Return a symbol for the variable of the given series ring.

```julia
max_precision(S::SeriesRing)
```

Return the precision cap of the given series ring.


```julia
precision(f::SeriesElem)
valuation(f::SeriesElem)
```

```julia
gen(R::SeriesRing)
```

The following functions are also provided for all series.

```julia
coeff(a::SeriesElem, n::Int)
```

Return the degree $n$ coefficient of the given power series. Note coefficients
are numbered from $n = 0$ for the constant coefficient. If $n$ exceeds the
current precision of the power series, the function returns a zero coefficient.

For power series types, $n$ must be non-negative. Laurent series do not have this
restriction.

```@docs
modulus(::SeriesElem{T}) where T <: ResElem
```

```@docs
is_gen(::RelPowerSeriesRingElem)
```

**Examples**

```jldoctest
julia> S, x = power_series_ring(ZZ, 10, :x)
(Univariate power series ring over integers, x + O(x^11))

julia> f = 1 + 3x + x^3 + O(x^10)
1 + 3*x + x^3 + O(x^10)

julia> g = 1 + 2x + x^2 + O(x^10)
1 + 2*x + x^2 + O(x^10)

julia> h = zero(S)
O(x^10)

julia> k = one(S)
1 + O(x^10)

julia> isone(k)
true

julia> iszero(f)
false

julia> n = pol_length(f)
4

julia> c = polcoeff(f, 3)
1

julia> U = base_ring(S)
Integers

julia> v = var(S)
:x

julia> max_precision(S) == 10
true

julia> T = parent(x + 1)
Univariate power series ring in x with precision 10
  over integers

julia> g == deepcopy(g)
true

julia> t = divexact(2g, 2)
1 + 2*x + x^2 + O(x^10)

julia> p = precision(f)
10

julia> R, t = power_series_ring(QQ, 10, :t)
(Univariate power series ring over rationals, t + O(t^11))

julia> S, x = power_series_ring(R, 30, :x)
(Univariate power series ring over R, x + O(x^31))

julia> a = O(x^4)
O(x^4)

julia> b = (t + 3)*x + (t^2 + 1)*x^2 + O(x^4)
(3 + t + O(t^10))*x + (1 + t^2 + O(t^10))*x^2 + O(x^4)

julia> k = is_gen(gen(R))
true

julia> m = is_unit(-1 + x + 2x^2)
true

julia> n = valuation(a)
4

julia> p = valuation(b)
1

julia> c = coeff(b, 2)
1 + t^2 + O(t^10)

julia> S, x = power_series_ring(ZZ, 10, :x)
(Univariate power series ring over integers, x + O(x^11))

julia> f = 1 + 3x + x^3 + O(x^5)
1 + 3*x + x^3 + O(x^5)

julia> g = S(BigInt[1, 2, 0, 1, 0, 0, 0], 4, 10, 3);

julia> set_length!(g, 3)
x^3 + 2*x^4 + O(x^10)

julia> g = setcoeff!(g, 2, BigInt(11))
x^3 + 2*x^4 + 11*x^5 + O(x^10)

julia> fit!(g, 8)

julia> g = setcoeff!(g, 7, BigInt(4))
x^3 + 2*x^4 + 11*x^5 + O(x^10)
```

### Change base ring

```@docs
map_coefficients(::Any, ::AbsPowerSeriesRingElem{<:RingElem})
change_base_ring(::Ring, ::AbsPowerSeriesRingElem{<:RingElem})
```

**Examples**

```jldoctest
julia> R, x = power_series_ring(ZZ, 10, :x)
(Univariate power series ring over integers, x + O(x^11))

julia> f = 4*x^6 + x^7 + 9*x^8 + 16*x^9 + 25*x^10 + O(x^11)
4*x^6 + x^7 + 9*x^8 + 16*x^9 + 25*x^10 + O(x^11)

julia> map_coefficients(AbstractAlgebra.sqrt, f)
2*x^6 + x^7 + 3*x^8 + 4*x^9 + 5*x^10 + O(x^11)

julia> change_base_ring(QQ, f)
4*x^6 + x^7 + 9*x^8 + 16*x^9 + 25*x^10 + O(x^11)
```

### Shifting

```@docs
shift_left(::RelPowerSeriesRingElem{T}, ::Int) where T <: RingElem
```

```@docs
shift_right(::RelPowerSeriesRingElem{T}, ::Int) where T <: RingElem
```

**Examples**

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S, x = power_series_ring(R, 30, :x)
(Univariate power series ring over R, x + O(x^31))

julia> a = 2x + x^3
2*x + x^3 + O(x^31)

julia> b = O(x^4)
O(x^4)

julia> c = 1 + x + 2x^2 + O(x^5)
1 + x + 2*x^2 + O(x^5)

julia> d = 2x + x^3 + O(x^4)
2*x + x^3 + O(x^4)

julia> f = shift_left(a, 2)
2*x^3 + x^5 + O(x^33)

julia> g = shift_left(b, 2)
O(x^6)

julia> h = shift_right(c, 1)
1 + 2*x + O(x^4)

julia> k = shift_right(d, 3)
1 + O(x^1)

```

### Truncation

```@docs
truncate(::RelPowerSeriesRingElem{T}, ::Int) where T <: RingElem
```

**Examples**

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S, x = power_series_ring(R, 30, :x)
(Univariate power series ring over R, x + O(x^31))

julia> a = 2x + x^3
2*x + x^3 + O(x^31)

julia> b = O(x^4)
O(x^4)

julia> c = 1 + x + 2x^2 + O(x^5)
1 + x + 2*x^2 + O(x^5)

julia> d = 2x + x^3 + O(x^4)
2*x + x^3 + O(x^4)

julia> f = truncate(a, 3)
2*x + O(x^3)

julia> g = truncate(b, 2)
O(x^2)

julia> h = truncate(c, 7)
1 + x + 2*x^2 + O(x^5)

julia> k = truncate(d, 5)
2*x + x^3 + O(x^4)

```

### Division

```@docs
Base.inv(::RelPowerSeriesRingElem)
```

**Examples**

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S, x = power_series_ring(R, 30, :x)
(Univariate power series ring over R, x + O(x^31))

julia> a = 1 + x + 2x^2 + O(x^5)
1 + x + 2*x^2 + O(x^5)

julia> b = S(-1)
-1 + O(x^30)

julia> c = inv(a)
1 - x - x^2 + 3*x^3 - x^4 + O(x^5)

julia> d = inv(b)
-1 + O(x^30)

```

### Composition

```@docs
compose(a::RelPowerSeriesRingElem, b::RelPowerSeriesRingElem)
```

Note that `subst` can be used instead of `compose`, however the provided
functionality is the same. General series substitution is not well-defined.

### Derivative and integral

```@docs
derivative(a::RelPowerSeriesRingElem)
```

```@docs
integral(a::RelPowerSeriesRingElem)
```

### Special functions

```@docs
Base.log(a::SeriesElem{T}) where T <: FieldElem
```

```@docs
Base.exp(a::RelPowerSeriesRingElem)
```
Methods for `is_square` and `sqrt` are provided for inputs of type `RelPowerSeriesRingElem`.


**Examples**

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S, x = power_series_ring(R, 30, :x)
(Univariate power series ring over R, x + O(x^31))

julia> T, z = power_series_ring(QQ, 30, :z)
(Univariate power series ring over rationals, z + O(z^31))

julia> a = 1 + z + 3z^2 + O(z^5)
1 + z + 3*z^2 + O(z^5)

julia> b = z + 2z^2 + 5z^3 + O(z^5)
z + 2*z^2 + 5*z^3 + O(z^5)

julia> c = exp(x + O(x^40))
1 + x + 1//2*x^2 + 1//6*x^3 + 1//24*x^4 + 1//120*x^5 + 1//720*x^6 + 1//5040*x^7 + 1//40320*x^8 + 1//362880*x^9 + 1//3628800*x^10 + 1//39916800*x^11 + 1//479001600*x^12 + 1//6227020800*x^13 + 1//87178291200*x^14 + 1//1307674368000*x^15 + 1//20922789888000*x^16 + 1//355687428096000*x^17 + 1//6402373705728000*x^18 + 1//121645100408832000*x^19 + 1//2432902008176640000*x^20 + 1//51090942171709440000*x^21 + 1//1124000727777607680000*x^22 + 1//25852016738884976640000*x^23 + 1//620448401733239439360000*x^24 + 1//15511210043330985984000000*x^25 + 1//403291461126605635584000000*x^26 + 1//10888869450418352160768000000*x^27 + 1//304888344611713860501504000000*x^28 + 1//8841761993739701954543616000000*x^29 + 1//265252859812191058636308480000000*x^30 + O(x^31)

julia> d = divexact(x, exp(x + O(x^40)) - 1)
1 - 1//2*x + 1//12*x^2 - 1//720*x^4 + 1//30240*x^6 - 1//1209600*x^8 + 1//47900160*x^10 - 691//1307674368000*x^12 + 1//74724249600*x^14 - 3617//10670622842880000*x^16 + 43867//5109094217170944000*x^18 - 174611//802857662698291200000*x^20 + 77683//14101100039391805440000*x^22 - 236364091//1693824136731743669452800000*x^24 + 657931//186134520519971831808000000*x^26 - 3392780147//37893265687455865519472640000000*x^28 + O(x^29)

julia> f = exp(b)
1 + z + 5//2*z^2 + 43//6*z^3 + 193//24*z^4 + O(z^5)

julia> log(exp(b)) == b
true

julia> h = sqrt(a)
1 + 1//2*z + 11//8*z^2 - 11//16*z^3 - 77//128*z^4 + O(z^5)

```

### Random generation

Random series can be constructed using the `rand` function. A range of possible
valuations is provided. The maximum precision of the ring is used as a bound on
the precision. Other parameters are used to construct random coefficients.

```julia
rand(R::SeriesRing, val_range::AbstractUnitRange{Int}, v...)
```

**Examples**

```jldoctest; setup = :(import Random; Random.seed!(42))
julia> R, x = power_series_ring(ZZ, 10, :x)
(Univariate power series ring over integers, x + O(x^11))

julia> f = rand(R, 3:5, -10:10)
3*x^4 - x^5 + 4*x^7 + 4*x^8 - 7*x^9 + 2*x^10 + 4*x^11 - x^12 - 4*x^13 + O(x^14)
```
