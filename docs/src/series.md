```@meta
CurrentModule = AbstractAlgebra
```

# Generic power series

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
`src/generic/RelSeries.jl` for which $j > 0$ in the above description, and Laurent
series where $j$ can be negative, implemented in `src/generic/Laurent.jl`.
Note that there are two implementations for Laurent series, one over rings and one
over fields, though in practice most of the implementation uses the same code in both
cases.

There is a single implementation of absolute series: absolute power series, implemented
in `src/generic/AbsSeries.jl`.

As well as implementing the Series Ring interface, the series modules in
AbstractAlgebra.jl implement the generic algorithms described below.

All of the generic functionality is part of the `Generic` submodule of
AbstractAlgebra.jl. This is exported by default so that it is not necessary to qualify
function names.

## Types and parent objects

The types of generic series implemented by AbstractAlgebra.jl are
`Generic.RelSeries{T}`, `Generic.AbsSeries{T}`, `Generic.LaurentSeriesRingElem{T}` and
`Generic.LaurentSeriesFieldElem{T}`.

Relative power series elements belong to the abstract type
`AbstractAlgebra.RelSeriesElem`.

Laurent series elements belong directly to either `AbstractAlgebra.RingElem` or
`AbstractAlgebra.FieldElem` since it is more useful to be able to distinguish whether
they belong to a ring or field than it is to distinguish that they are relative series.

Absolute power series elements belong to `AbstractAlgebra.AbsSeriesElem`.

The parent types for relative and absolute power series, `Generic.RelSeriesRing{T}` 
and `Generic.AbsSeriesRing{T}` respectively, belong to `AbstractAlgebra.SeriesRing{T}`.

The parent types for Laurent series rings and fields, `Generic.LaurentSeriesRing{T}`
and `Generic.LaurentSeriesField{T}` respectively, belong directly to 
`AbstractAlgebra.Ring` and `AbstractAlgebra.Field` respectively.

The default precision, string representation of the variable and base ring $R$ of a
generic power series are stored in its parent object. 

## Series ring constructors

In order to construct series in AbstractAlgebra.jl, one must first construct the ring
itself. This is accomplished with any of the following constructors.

```julia
PowerSeriesRing(R::AbstractAlgebra.Ring, prec_max::Int, s::AbstractString; cached::Bool = true, model=:capped_relative)
```

```julia
LaurentSeriesRing(R::AbstractAlgebra.Ring, prec_max::Int, s::AbstractString; cached::Bool = true)
```

```julia
LaurentSeriesRing(R::AbstractAlgebra.Field, prec_max::Int, s::AbstractString; cached::Bool = true)
```

Given a base ring `R`, a maximum precision (relative or absolute, depending on the
model) and a string `s` specifying how the generator (variable) should be
printed, return a typle `S, x` representing the series ring and its generator.

By default, `S` will depend only on `S`, `x` and the maximum precision and will be
cached. Setting the optional argument `cached` to `false` will prevent this.

In the case of power series, the optional argument `model` can be set to either
`:capped_absolute` or `capped_relative`, depending on which power series model is
required.

Here are some examples of constructing various kinds of series rings and coercing
various elements into those rings.

**Examples**

```julia
R, x = PowerSeriesRing(ZZ, 10, "x")
S, y = PowerSeriesRing(ZZ, 10, "y"; model=:capped_absolute)
T, z = LaurentSeriesRing(ZZ, 10, "z")
U, w = LaurentSeriesField(QQ, 10, "w")

f = R()
g = S(123)
h = U(BigInt(1234))
k = T(z + 1)
```

## Big-oh notation

Series elements can be given a precision using the big-oh notation. This is provided
by a function of the following form, (or something equivalent for Laurent series):

```julia
O(x::SeriesElem)
```

**Examples**

```julia
R, x = PowerSeriesRing(ZZ, 10, "x")
S, y = LaurentSeriesRing(ZZ, 10, "y")

f = 1 + 2x + O(x^5)
g = 2y + 7y^2 + O(y^7)
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

## Basic ring functionality

All power series models provide the functionality described in the Ring and Series Ring
interfaces.

**Examples**

```julia
S, x = PowerSeriesRing(ZZ, 10, "x")

f = 1 + 3x + x^3 + O(x^10)
g = 1 + 2x + x^2 + O(x^10)

h = zero(S)
k = one(S)
isone(k) == true
iszero(f) == false
n = pol_length(f)
c = polcoeff(f, 3)
U = base_ring(S)
v = var(S)
T = parent(x + 1)
g == deepcopy(g)
t = divexact(2g, 2)
p = precision(f)
```

## Series functionality provided by AbstractAlgebra.jl


The functionality below is automatically provided by AbstractAlgebra.jl for any series
module that implements the full Series Ring interface. This includes AbstractAlgebra's
own generic series rings.

Of course, modules are encouraged to provide specific implementations of the functions
described here, that override the generic implementation.

Unless otherwise noted, the functions are available for all series models, including
Laurent series. We denote this by using the abstract type
`AbstractAlgebra.RelSeriesElem`, even though absolute series and Laurent series types
do not belong to this abstract type.

### Basic functionality

```
coeff(a::AbstractAlgebra.SeriesElem, n::Int)
```

Return the degree $n$ coefficient of the given power series. Note coefficients
are numbered from $n = 0$ for the constant coefficient. If $n$ exceeds the
current precision of the power series, the function returns a zero coefficient.

For power series types, $n$ must be non-negative. Laurent series do not have this
restriction.

```@docs
modulus{T <: ResElem}(::SeriesElem{T})
```

```@docs
isgen(::RelSeriesElem)
```

```@docs
isunit(::RelSeriesElem)
```

**Examples**

```julia
R, t = PowerSeriesRing(QQ, 10, "t")
S, x = PowerSeriesRing(R, 30, "x")

a = O(x^4)
b = (t + 3)*x + (t^2 + 1)*x^2 + O(x^4)

k = isgen(gen(R))
m = isunit(-1 + x + 2x^2)
n = valuation(a)
p = valuation(b)
c = coeff(b, 2)
```

### Shifting

```@docs
shift_left{T <: RingElem}(::RelSeriesElem{T}, ::Int)
```

```@docs
shift_right{T <: RingElem}(::RelSeriesElem{T}, ::Int)
```

**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S, x = PowerSeriesRing(R, 30, "x")

a = 2x + x^3
b = O(x^4)
c = 1 + x + 2x^2 + O(x^5)
d = 2x + x^3 + O(x^4)

f = shift_left(a, 2)
g = shift_left(b, 2)
h = shift_right(c, 1)
k = shift_right(d, 3)
```

### Truncation

```@docs
truncate{T <: RingElem}(::RelSeriesElem{T}, ::Int)
```

**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S, x = PowerSeriesRing(R, 30, "x")

a = 2x + x^3
b = O(x^4)
c = 1 + x + 2x^2 + O(x^5)
d = 2x + x^3 + O(x^4)

f = truncate(a, 3)
g = truncate(b, 2)
h = truncate(c, 7)
k = truncate(d, 5)
```

### Division

```@docs
inv(::RelSeriesElem)
```

**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S, x = PowerSeriesRing(R, 30, "x")

a = 1 + x + 2x^2 + O(x^5)
b = S(-1)

c = inv(a)
d = inv(b)
```

### Special functions

```@docs
Base.exp(a::RelSeriesElem)
```

```@docs
Base.sqrt(a::RelSeriesElem)
```


**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S, x = PowerSeriesRing(R, 30, "x")
T, z = PowerSeriesRing(QQ, 30, "z")

a = 1 + z + 3z^2 + O(z^5)
b = z + 2z^2 + 5z^3 + O(z^5)

c = exp(x + O(x^40))
d = divexact(x, exp(x + O(x^40)) - 1)
f = exp(b)
h = sqrt(a)
```
