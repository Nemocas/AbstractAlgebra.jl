```@meta
CurrentModule = AbstractAlgebra
```

# Generic Puiseux series

AbstractAlgebra.jl allows the creation of Puiseux series over any computable commutative
ring $R$.

Puiseux series are power series of the form
$a_jx^{j/m} + a_{j+1}x^{(j+1)/m} + \cdots + a_{k-1}x^{(k-1)/m} + O(x^{k/m})$
for some integer $m > 0$ where $i \geq 0$, $a_i \in R$ and the relative precision
$k - j$ is at most equal to some specified precision $n$.

The generic Puiseux series module is implemented in `src/generic/PuiseuxSeries.jl`.

As well as implementing the Series Ring interface, the Puiseux series module in
AbstractAlgebra.jl implements the generic algorithms described below.

All of the generic functionality is part of the `Generic` submodule of
AbstractAlgebra.jl. This is exported by default so that it is not necessary to qualify
function names.

## Types and parent objects

The types of generic polynomials implemented by AbstractAlgebra.jl are
`Generic.PuiseuxSeriesRingElem{T}` and `Generic.PuiseuxSeriesFieldElem{T}`.

Both series element types belong to the union type
`Generic.PuiseuxSeriesElem`.

Puiseux series elements belong directly to either `AbstractAlgebra.RingElem` or
`AbstractAlgebra.FieldElem` since it is more useful to be able to distinguish whether
they belong to a ring or field than it is to distinguish that they are Puiseux series.

The parent types for Puiseux series, `Generic.PuiseuxSeriesRing{T}` 
and `Generic.PuiseuxSeriesField{T}` respectively, belong to `AbstractAlgebra.Ring`
and `AbstractAlgebra.Field` respectively.

The default precision, string representation of the variable and base ring $R$ of a
generic Puiseux series are stored in its parent object. 

## Puisuex series ring constructors

In order to construct Puiseux series in AbstractAlgebra.jl, one must first construct the
ring itself. This is accomplished with any of the following constructors.

```julia
PuiseuxSeriesRing(R::AbstractAlgebra.Ring, prec_max::Int, s::AbstractString; cached::Bool = true)
```

```julia
PuiseuxSeriesRing(R::AbstractAlgebra.Field, prec_max::Int, s::AbstractString; cached::Bool = true)
```

```julia
PuiseuxSeriesField(R::AbstractAlgebra.Field, prec_max::Int, s::AbstractString; cached::Bool = true)
```

Given a base ring `R`, a maximum relative precision and a string `s` specifying how the
generator (variable) should be printed, return a typle `S, x` representing the Puiseux
series ring and its generator.

By default, `S` will depend only on `S`, `x` and the maximum precision and will be
cached. Setting the optional argument `cached` to `false` will prevent this.

Here are some examples of constructing various kinds of Puiseux series rings and coercing
various elements into those rings.

**Examples**

```julia
R, x = PuiseuxSeriesRing(ZZ, 10, "x")
S, y = PuiseuxSeriesField(QQ, 10, "y")

f = R()
g = S(123)
h = R(BigInt(1234))
k = S(y + 1)
```

## Big-oh notation

Series elements can be given a precision using the big-oh notation. This is provided
by a function of the following form, (or something equivalent for Laurent series):

```julia
O(x::SeriesElem)
```

**Examples**

```julia
R, x = PuiseuxSeriesRing(ZZ, 10, "x")

f = 1 + 2x + O(x^5)
g = 2x^(1//3) + 7y^(2//3) + O(y^(7//3))
```

What is happening here in practice is that `O(x^n)` is creating the series `0 + O(x^n)`
and the rules for addition of series dictate that if this is added to a series of 
greater precision, then the lower of the two precisions must be used.

Of course it may be that the precision of the series that `O(x^n)` is added to is
already lower than `n`, in which case adding `O(x^n)` has no effect. This is the case
if the default precision is too low, since `x` on its own has the default precision.

## Puiseux series implementation

Puiseux series have their maximum relative precision capped at some value `prec_max`.
This refers to the internal Laurent series used to store the Puiseux series, i.e. the
series without denominators in the exponents.

The Puiseux series type stores such a Laurent series and a `scale` or denominator for
the exponents. For example, $f(x) = 1 + x^{1/3} + 2x^{2/3} + O(x^{7/3})$ would be
stored as a Laurent series $1 + x + 2x^2 + O(x^7)$ and a scale of $3$..

The maximum precision is also used as a default (Laurent) precision in the case
of coercing coefficients into the ring and for any computation where the result could
mathematically be given to infinite precision.

In all models we say that two Puiseux series are equal if they agree up to the minimum
**absolute** precision of the two power series.

Thus, for example, $x^5 + O(x^{10}) == 0 + O(x^5)$, since the minimum absolute
precision is $5$.

Sometimes it is necessary to compare Puiseux series not just for arithmetic
equality, as above, but to see if they have precisely the same precision and
terms. For this purpose we introduce the `isequal` function.

For example, if $f = x^2 + O(x^7)$ and $g = x^2 + O(x^8)$ and $h = 0 + O(x^2)$
then $f == g$, $f == h$ and $g == h$, but `isequal(f, g)`, `isequal(f, h)` and
`isequal(g, h)` would all return `false`. However, if $k = x^2 + O(x^7)$ then
`isequal(f, k)` would return `true`.

There are a number of technicalities that must be observed when working with Puiseux
series. As these are the same as for the other series rings in AbstractAlgebra.jl,
we refer the reader to the documentation of series rings for information about these
issues.

## Basic ring functionality

All Puiseux series provide the functionality described in the Ring and Series Ring
interfaces with the exception of the `pol_length` and `polcoeff` functions.

**Examples**

```julia
S, x = PuiseuxSeriesRing(ZZ, 10, "x")

f = 1 + 3x + x^3 + O(x^10)
g = 1 + 2x^(1//3) + x^(2//3) + O(x^(7//3))

h = zero(S)
k = one(S)
isone(k) == true
iszero(f) == false
U = base_ring(S)
v = var(S)
T = parent(x + 1)
g == deepcopy(g)
t = divexact(2g, 2)
p = precision(f)
```

## Puiseux series functionality provided by AbstractAlgebra.jl

The functionality below is automatically provided by AbstractAlgebra.jl for any Puiseux
series.

Of course, modules are encouraged to provide specific implementations of the functions
described here, that override the generic implementation.

### Basic functionality

```
coeff(a::Generic.PuiseuxSeriesElem, n::Int)
```

```
coeff(a::Generic.PuiseuxSeriesElem, n::Rational{Int})
```

Return the coefficient of the term of exponent $n$ of the given power series. If $n$
exceeds the current precision of the power series or does not correspond to a nonzero
term of the Puiseux series, the function returns a zero coefficient.

```@docs
modulus{T <: ResElem}(::Generic.PuiseuxSeriesElem{T})
```

```@docs
isgen(::Generic.PuiseuxSeriesElem)
```

```@docs
isunit(::Generic.PuiseuxSeriesElem)
```

**Examples**

```julia
R, t = PuiseuxSeriesRing(QQ, 10, "t")
S, x = PuiseuxSeriesRing(R, 30, "x")

a = O(x^4)
b = (t + 3)*x + (t^2 + 1)*x^2 + O(x^4)

k = isgen(gen(R))
m = isunit(-1 + x^(1//3) + 2x^2)
n = valuation(a)
p = valuation(b)
c = coeff(b, 2)
```

### Division

```@docs
inv(::Generic.PuiseuxSeriesElem)
```

**Examples**

```julia
R, x = PuiseuxSeriesRing(QQ, 30, "x")

a = 1 + x + 2x^2 + O(x^5)
b = R(-1)

c = inv(a)
d = inv(b)
```

### Special functions

```@docs
Base.exp(a::Generic.PuiseuxSeriesElem)
```

```@docs
Base.sqrt(a::Generic.PuiseuxSeriesElem)
```


**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S, x = PuiseuxSeriesRing(R, 30, "x")
T, z = PuiseuxSeriesRing(QQ, 30, "z")

a = 1 + z + 3z^2 + O(z^5)
b = z + 2z^2 + 5z^3 + O(z^5)

c = exp(x + O(x^40))
d = divexact(x, exp(x + O(x^40)) - 1)
f = exp(b)
h = sqrt(a)
```
