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

AbstractAlgebra.jl implements a generic fraction type `Generic.FracFieldElem{T}`
where `T` is the type of elements of the base ring. See the file
`src/generic/GenericTypes.jl` for details.

Parent objects of such fraction elements have type `Generic.FracField{T}`.

## Factored fraction types

AbstractAlgebra.jl also implements a fraction type `Generic.FactoredFracFieldElem{T}`
with parent objects of such fractions having type `Generic.FactoredFracField{T}`.
As opposed to the fractions of type `Generic.FracFieldElem{T}`, which are just a
numerator and denominator, these fractions are maintained in factored form as
much as possible.

## Abstract types

All fraction element types belong to the abstract type `FracElem{T}`
and the fraction field types belong to the abstract type `FracField{T}`.
This enables one to write generic functions that can accept any AbstractAlgebra
fraction type.

!!! note

    Both the generic fraction field type `Generic.FracField{T}` and the abstract
    type it belongs to, `FracField{T}` are both called `FracField`. The
    former is a (parameterised) concrete type for a fraction field over a given base ring
    whose elements have type `T`. The latter is an abstract type representing all
    fraction field types in AbstractAlgebra.jl, whether generic or very specialised (e.g.
    supplied by a C library).

## Fraction field constructors

In order to construct fractions in AbstractAlgebra.jl, one can first construct the
fraction field itself. This is accomplished with the following constructor.

```julia
fraction_field(R::Ring; cached::Bool = true)
```

Given a base ring `R` return the parent object of the fraction field of $R$. By default
the parent object `S` will depend only on `R` and will be cached. Setting the optional
argument `cached` to `false` will prevent the parent object `S` from being cached.

Here are some examples of creating fraction fields and making use of the
resulting parent objects to coerce various elements into the fraction field.

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, "x")
(Univariate polynomial ring in x over integers, x)

julia> S = fraction_field(R)
Fraction field
  of univariate polynomial ring in x over integers

julia> f = S()
0

julia> g = S(123)
123

julia> h = S(BigInt(1234))
1234

julia> k = S(x + 1)
x + 1
```

## Factored Fraction field constructors

The corresponding factored field uses the following constructor.

```julia
FactoredFractionField(R::Ring; cached::Bool = true)
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> S = FactoredFractionField(R)
Factored fraction field of Multivariate polynomial ring in 2 variables over integers

julia> (X, Y) = (S(x), S(y))
(x, y)

julia> f = X^6*(X+Y)^2*(X^2+Y)^3*(X+2*Y)^-3*(X+3*Y)^-4
x^6*(x + y)^2*(x^2 + y)^3/((x + 2*y)^3*(x + 3*y)^4)

julia> numerator(f)
x^14 + 2*x^13*y + x^12*y^2 + 3*x^12*y + 6*x^11*y^2 + 3*x^10*y^3 + 3*x^10*y^2 + 6*x^9*y^3 + 3*x^8*y^4 + x^8*y^3 + 2*x^7*y^4 + x^6*y^5

julia> denominator(f)
x^7 + 18*x^6*y + 138*x^5*y^2 + 584*x^4*y^3 + 1473*x^3*y^4 + 2214*x^2*y^5 + 1836*x*y^6 + 648*y^7

julia> derivative(f, x)
x^5*(x + y)*(x^2 + y)^2*(7*x^5 + 58*x^4*y + 127*x^3*y^2 + x^3*y + 72*x^2*y^3 + 22*x^2*y^2 + 61*x*y^3 + 36*y^4)/((x + 2*y)^4*(x + 3*y)^5)
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
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> S = fraction_field(R)
Fraction field
  of univariate polynomial ring in x over rationals

julia> f = S(x + 1)
x + 1

julia> g = (x^2 + x + 1)//(x^3 + 3x + 1)
(x^2 + x + 1)//(x^3 + 3*x + 1)

julia> x//f
x//(x + 1)

julia> f//x
(x + 1)//x
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

Return the characteristic of the base ring of the fraction field. If the
characteristic is not known an exception is raised.


**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> S = fraction_field(R)
Fraction field
  of univariate polynomial ring in x over rationals

julia> f = S(x + 1)
x + 1

julia> U = base_ring(S)
Univariate polynomial ring in x over rationals

julia> V = base_ring(f)
Univariate polynomial ring in x over rationals

julia> T = parent(f)
Fraction field
  of univariate polynomial ring in x over rationals

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
is_unit(f::FracElem)
```

And they implement the fraction field interface.

```julia
numerator(a::FracElem)
denominator(a::FracElem)
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> S = fraction_field(R)
Fraction field
  of univariate polynomial ring in x over rationals

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
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> f = (x + 1)//(x^3 + 3x + 1)
(x + 1)//(x^3 + 3*x + 1)

julia> g = (x^2 + 2x + 1)//(x^2 + x + 1)
(x^2 + 2*x + 1)//(x^2 + x + 1)

julia> h = gcd(f, g)
(x + 1)//(x^5 + x^4 + 4*x^3 + 4*x^2 + 4*x + 1)

```

### Square root

```@docs
is_square{T <: RingElem}(::FracElem{T})
```

```@docs
Base.sqrt(::FracElem{T}) where {T <: RingElem}
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> S = fraction_field(R)
Fraction field
  of univariate polynomial ring in x over rationals

julia> a = (21//4*x^6 - 15*x^5 + 27//14*x^4 + 9//20*x^3 + 3//7*x + 9//10)//(x + 3)
(21//4*x^6 - 15*x^5 + 27//14*x^4 + 9//20*x^3 + 3//7*x + 9//10)//(x + 3)

julia> sqrt(a^2)
(21//4*x^6 - 15*x^5 + 27//14*x^4 + 9//20*x^3 + 3//7*x + 9//10)//(x + 3)

julia> is_square(a^2)
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
julia> R, x = polynomial_ring(ZZ, "x")
(Univariate polynomial ring in x over integers, x)

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

**Examples**

```jldoctest; setup = :(import Random; Random.seed!(42))
julia> K = fraction_field(ZZ)
Rationals

julia> f = rand(K, -10:10)
-1//3

julia> R, x = polynomial_ring(ZZ, "x")
(Univariate polynomial ring in x over integers, x)

julia> S = fraction_field(R)
Fraction field
  of univariate polynomial ring in x over integers

julia> g = rand(S, -1:3, -10:10)
(-4*x - 4)//(4*x^2 + x - 4)
```

### Extra functionality for factored fractions

The `Generic.FactoredFracFieldElem{T}` type implements an interface similar to that of
the `Fac{T}` type for iterating over the terms in the factorisation. There is
also the function `push_term!(a, b, e)` for efficiently performing `a *= b^e`,
and the function `normalise` returns relatively prime terms.

**Examples**

```jldoctest
julia> F = FactoredFractionField(ZZ)
Factored fraction field of Integers

julia> f = F(-1)
-1

julia> push_term!(f, 10, 10)
-10^10

julia> push_term!(f, 42, -8)
-10^10/42^8

julia> normalise(f)
-5^10*2^2/21^8

julia> unit(f)
-1

julia> collect(f)
2-element Vector{Tuple{BigInt, Int64}}:
 (10, 10)
 (42, -8)
```

