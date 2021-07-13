```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Generic fraction fields

AbstractAlgebra.jl provides a module, implemented in `src/Fraction.jl` for
fraction fields over any gcd domain belonging to the AbstractAlgebra.jl
abstract type hierarchy.

## Generic fraction types

AbstractAlgebra.jl implements a generic fraction type `Generic.Frac{T}`
where `T` is the type of elements of the base ring. See the file
`src/generic/GenericTypes.jl` for details.

Parent objects of such fraction elements have type `Generic.FracField{T}`.

## Abstract types

All fraction element types belong to the abstract type `FracElem{T}`
and the fraction field types belong to the abstract type `FracField{T}`.
This enables one to write generic functions that can accept any AbstractAlgebra
fraction type.

Note that both the generic fraction field type `Generic.FracField{T}` and the abstract
type it belongs to, `FracField{T}` are both called `FracField`. The
former is a (parameterised) concrete type for a fraction field over a given base ring
whose elements have type `T`. The latter is an abstract type representing all
fraction field types in AbstractAlgebra.jl, whether generic or very specialised (e.g.
supplied by a C library).

## Fraction field constructors

In order to construct fractions in AbstractAlgebra.jl, one can first construct the
fraction field itself. This is accomplished with the following constructor.

```julia
FractionField(R::Ring; cached::Bool = true)
```

Given a base ring `R` return the parent object of the fraction field of $R$. By default
the parent object `S` will depend only on `R` and will be cached. Setting the optional
argument `cached` to `false` will prevent the parent object `S` from being cached.

Here are some examples of creating fraction fields and making use of the
resulting parent objects to coerce various elements into the fraction field.

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S = FractionField(R)
Fraction field of Univariate Polynomial Ring in x over Integers

julia> f = S()
0

julia> g = S(123)
123

julia> h = S(BigInt(1234))
1234

julia> k = S(x + 1)
x + 1

```

## Fraction constructors

One can construct fractions using the fraction field parent object, as for any
ring or field.

```julia
(R::FracField)() # constructs zero
(R::FracField)(c::Integer)
(R::FracField)(c::elem_type(R))
(R::FracField{T})(a::T) where T <: RingElement
```

One may also use the Julia double slash operator to construct elements of the
fraction field without constructing the fraction field parent first.

```julia
//(x::T, y::T) where T <: RingElement
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> S = FractionField(R)
Fraction field of Univariate Polynomial Ring in x over Rationals

julia> f = S(x + 1)
x + 1

julia> g = (x^2 + x + 1)//(x^3 + 3x + 1)
(x^2 + x + 1)//(x^3 + 3*x + 1)
```

## Functions for types and parents of fraction fields

Fraction fields in AbstractAlgebra.jl implement the Ring interface.

```julia
base_ring(R::FracField)
base_ring(a::FracElem)
```

Return the base ring of which the fraction field was constructed.

```julia
parent(a::FracElem)
```

Return the fraction field of the given fraction.

```julia
characteristic(R::FracField)
```

Return the characteristic of the base ring of the fraction field.


**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> S = FractionField(R)
Fraction field of Univariate Polynomial Ring in x over Rationals

julia> f = S(x + 1)
x + 1

julia> U = base_ring(S)
Univariate Polynomial Ring in x over Rationals

julia> V = base_ring(f)
Univariate Polynomial Ring in x over Rationals

julia> T = parent(f)
Fraction field of Univariate Polynomial Ring in x over Rationals

julia> m = characteristic(S)
0
```

## Fraction field functions

### Basic functions

Fraction fields implement the Ring interface.

```julia
zero(R::FracField)
one(R::FracField)
iszero(a::FracElem)
isone(a::FracElem)
```

```julia
inv(a::T) where T <: FracElem
```

They also implement the field interface.

```julia
isunit(f::FracElem)
```

And they implement the fraction field interface.

```julia
numerator(a::FracElem)
denominator(a::FracElem)
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> S = FractionField(R)
Fraction field of Univariate Polynomial Ring in x over Rationals

julia> f = S(x + 1)
x + 1

julia> g = (x^2 + x + 1)//(x^3 + 3x + 1)
(x^2 + x + 1)//(x^3 + 3*x + 1)

julia> h = zero(S)
0

julia> k = one(S)
1

julia> isone(k)
true

julia> iszero(f)
false

julia> r = deepcopy(f)
x + 1

julia> n = numerator(g)
x^2 + x + 1

julia> d = denominator(g)
x^3 + 3*x + 1
```

### Greatest common divisor

```@docs
gcd{T <: RingElem}(::FracElem{T}, ::FracElem{T})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> f = (x + 1)//(x^3 + 3x + 1)
(x + 1)//(x^3 + 3*x + 1)

julia> g = (x^2 + 2x + 1)//(x^2 + x + 1)
(x^2 + 2*x + 1)//(x^2 + x + 1)

julia> h = gcd(f, g)
(x + 1)//(x^5 + x^4 + 4*x^3 + 4*x^2 + 4*x + 1)

```

### Square root

```@docs
issquare{T <: RingElem}(::FracElem{T})
```

```@docs
Base.sqrt(::FracElem{T}) where {T <: RingElem}
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> S = FractionField(R)
Fraction field of Univariate Polynomial Ring in x over Rationals

julia> a = (21//4*x^6 - 15*x^5 + 27//14*x^4 + 9//20*x^3 + 3//7*x + 9//10)//(x + 3)
(21//4*x^6 - 15*x^5 + 27//14*x^4 + 9//20*x^3 + 3//7*x + 9//10)//(x + 3)

julia> sqrt(a^2)
(21//4*x^6 - 15*x^5 + 27//14*x^4 + 9//20*x^3 + 3//7*x + 9//10)//(x + 3)

julia> issquare(a^2)
true
```

### Remove and valuation

When working over a Euclidean domain, it is convenient to extend valuations to the
fraction field. To facilitate this, we define the following functions.

```@docs
remove{T <: RingElem}(::FracElem{T}, ::T)
```

```@docs
valuation{T <: RingElem}(::FracElem{T}, ::T)
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> f = (x + 1)//(x^3 + 3x + 1)
(x + 1)//(x^3 + 3*x + 1)

julia> g = (x^2 + 1)//(x^2 + x + 1)
(x^2 + 1)//(x^2 + x + 1)

julia> v, q = remove(f^3*g, x + 1)
(3, (x^2 + 1)//(x^11 + x^10 + 10*x^9 + 12*x^8 + 39*x^7 + 48*x^6 + 75*x^5 + 75*x^4 + 66*x^3 + 37*x^2 + 10*x + 1))

julia> v = valuation(f^3*g, x + 1)
3

```

### Random generation

Random fractions can be generated using `rand`. The parameters passed after the
fraction field tell `rand` how to generate random elements of the base ring.

```julia
rand(R::FracField, v...)
```

** Examples **

```@repl
K = FractionField(ZZ)
f = rand(K, -10:10)

R, x = PolynomialRing(ZZ, "x")
S = FractionField(R)
g = rand(S, -1:3, -10:10)
```
