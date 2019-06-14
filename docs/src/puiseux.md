```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
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

```jldoctest
julia> R, x = PuiseuxSeriesRing(ZZ, 10, "x")
(Puiseux series ring in x over Integers, x+O(x^11))

julia> S, y = PuiseuxSeriesField(QQ, 10, "y")
(Puiseux series field in y over Rationals, y+O(y^11))

julia> f = R()
0+O(x^10)

julia> g = S(123)
123//1+O(y^10)

julia> h = R(BigInt(1234))
1234+O(x^10)

julia> k = S(y + 1)
1//1+y+O(y^10)

```

## Big-oh notation

Series elements can be given a precision using the big-oh notation. This is provided
by a function of the following form, (or something equivalent for Laurent series):

```julia
O(x::SeriesElem)
```

**Examples**

```jldoctest
julia> R, x = PuiseuxSeriesRing(ZZ, 10, "x")
(Puiseux series ring in x over Integers, x+O(x^11))

julia> f = 1 + 2x + O(x^5)
1+2*x+O(x^5)

julia> g = 2x^(1//3) + 7x^(2//3) + O(x^(7//3))
2*x^1//3+7*x^2//3+O(x^7//3)
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

```jldoctest
julia> S, x = PuiseuxSeriesRing(ZZ, 10, "x")
(Puiseux series ring in x over Integers, x+O(x^11))

julia> f = 1 + 3x + x^3 + O(x^10)
1+3*x+x^3+O(x^10)

julia> g = 1 + 2x^(1//3) + x^(2//3) + O(x^(7//3))
1+2*x^1//3+x^2//3+O(x^7//3)

julia> h = zero(S)
0+O(x^10)

julia> k = one(S)
1+O(x^10)

julia> isone(k)
true

julia> iszero(f)
false

julia> U = base_ring(S)
Integers

julia> v = var(S)
:x

julia> T = parent(x + 1)
Puiseux series ring in x over Integers

julia> g == deepcopy(g)
true

julia> t = divexact(2g, 2)
1+2*x^1//3+x^2//3+O(x^7//3)

julia> p = precision(f)
10//1

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

```jldoctest
julia> R, t = PuiseuxSeriesRing(QQ, 10, "t")
(Puiseux series field in t over Rationals, t+O(t^11))

julia> S, x = PuiseuxSeriesRing(R, 30, "x")
(Puiseux series field in x over Puiseux series field in t over Rationals, x+O(x^31))

julia> a = O(x^4)
0//1+O(t^10)+O(x^4)

julia> b = (t + 3)*x + (t^2 + 1)*x^2 + O(x^4)
(3//1+t+O(t^10))*x+(1//1+t^2+O(t^10))*x^2+O(x^4)

julia> k = isgen(gen(R))
true

julia> m = isunit(-1 + x^(1//3) + 2x^2)
true

julia> n = valuation(a)
4//1

julia> p = valuation(b)
1//1

julia> c = coeff(b, 2)
1//1+t^2+O(t^10)

```

### Division

```@docs
inv(::Generic.PuiseuxSeriesElem)
```

**Examples**

```jldoctest
julia> R, x = PuiseuxSeriesRing(QQ, 30, "x")
(Puiseux series field in x over Rationals, x+O(x^31))

julia> a = 1 + x + 2x^2 + O(x^5)
1//1+x+2//1*x^2+O(x^5)

julia> b = R(-1)
-1//1+O(x^30)

julia> c = inv(a)
1//1-x-x^2+3//1*x^3-x^4+O(x^5)

julia> d = inv(b)
-1//1+O(x^30)

```

### Special functions

```@docs
Base.exp(a::Generic.PuiseuxSeriesElem)
```

```@docs
Base.sqrt(a::Generic.PuiseuxSeriesElem)
```


**Examples**

```jldoctest
julia> R, t = PolynomialRing(QQ, "t")
(Univariate Polynomial Ring in t over Rationals, t)

julia> S, x = PuiseuxSeriesRing(R, 30, "x")
(Puiseux series ring in x over Univariate Polynomial Ring in t over Rationals, x+O(x^31))

julia> T, z = PuiseuxSeriesRing(QQ, 30, "z")
(Puiseux series field in z over Rationals, z+O(z^31))

julia> a = 1 + z + 3z^2 + O(z^5)
1//1+z+3//1*z^2+O(z^5)

julia> b = z + 2z^2 + 5z^3 + O(z^5)
z+2//1*z^2+5//1*z^3+O(z^5)

julia> c = exp(x + O(x^40))
1//1+x+1//2*x^2+1//6*x^3+1//24*x^4+1//120*x^5+1//720*x^6+1//5040*x^7+1//40320*x^8+1//362880*x^9+1//3628800*x^10+1//39916800*x^11+1//479001600*x^12+1//6227020800*x^13+1//87178291200*x^14+1//1307674368000*x^15+1//20922789888000*x^16+1//355687428096000*x^17+1//6402373705728000*x^18+1//121645100408832000*x^19+1//2432902008176640000*x^20+1//51090942171709440000*x^21+1//1124000727777607680000*x^22+1//25852016738884976640000*x^23+1//620448401733239439360000*x^24+1//15511210043330985984000000*x^25+1//403291461126605635584000000*x^26+1//10888869450418352160768000000*x^27+1//304888344611713860501504000000*x^28+1//8841761993739701954543616000000*x^29+1//265252859812191058636308480000000*x^30+O(x^31)

julia> d = divexact(x, exp(x + O(x^40)) - 1)
1//1-1//2*x+1//12*x^2-1//720*x^4+1//30240*x^6-1//1209600*x^8+1//47900160*x^10-691//1307674368000*x^12+1//74724249600*x^14-3617//10670622842880000*x^16+43867//5109094217170944000*x^18-174611//802857662698291200000*x^20+77683//14101100039391805440000*x^22-236364091//1693824136731743669452800000*x^24+657931//186134520519971831808000000*x^26-3392780147//37893265687455865519472640000000*x^28+O(x^29)

julia> f = exp(b)
1//1+z+5//2*z^2+43//6*z^3+193//24*z^4+O(z^5)

julia> h = sqrt(a)
1//1+1//2*z+11//8*z^2-11//16*z^3-77//128*z^4+O(z^5)

```
