```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Rational function fields

AbstractAlgebra.jl provides a module, implemented in
`src/generic/RationalFunctionField.jl` for rational function fields $k(x)$ over
a field $k$.

All of this functionality is part of a submodule of AbstractAlgebra called
`Generic`. This is exported by default so that it is not necessary to qualify the
function names with the submodule name.

## Types and parent objects

Rational functions in $k(x)$ have type `Generic.Rat{T}` where `T` is the type of
elements of the coefficient field $k$. See the file
`src/generic/GenericTypes.jl` for details.

Parent objects corresponding to the rational function field $k$ have type
`Generic.RationalFunctionField{T}`.

The rational function types belong to the abstract type `AbstractAlgebra.Field`
and the rational function field types belong to the abstract type
`AbstractAlgebra.FieldElem`.

## Rational function field constructors

In order to construct rational functions in AbstractAlgebra.jl, one can first
construct the function field itself. This is accomplished with the following
constructor.

```julia
RationalFunctionField(k::AbstractAlgebra.Field, s::AbstractString; cached::Bool = true)
```

Given a coefficient field `k` return a tuple `(S, x)` consisting of the parent
object of the rational function field over $k$ and the generator `x`. By default
the parent object `S` will depend only on `R` and `s` and will be cached.
Setting the optional argument `cached` to `false` will prevent the parent object
`S` from being cached.

Here are some examples of creating rational function fields and making use of the
resulting parent objects to coerce various elements into the function field.

**Examples**

```jldoctest
julia> S, x = RationalFunctionField(QQ, "x")
(Rational function field over Rationals, x)

julia> f = S()
0

julia> g = S(123)
123

julia> h = S(BigInt(1234))
1234

julia> k = S(x + 1)
x + 1

```

## Basic rational function field functionality

Fraction fields in AbstractAlgebra.jl implement the full Field interface and
the entire fraction field interface.

We give some examples of such functionality.

**Examples**

```jldoctest
julia> S, x = RationalFunctionField(QQ, "x")
(Rational function field over Rationals, x)

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

julia> m = characteristic(S)
0

julia> U = base_ring(S)
Rationals

julia> V = base_ring(f)
Rationals

julia> T = parent(f)
Rational function field over Rationals

julia> r = deepcopy(f)
x + 1

julia> n = numerator(g)
x^2 + x + 1

julia> d = denominator(g)
x^3 + 3*x + 1

```

Note that numerator and denominator are returned as elements of a polynomial
ring whose variable is printed the same way as that of the generator of the
rational function field.

## Rational function field functionality provided by AbstractAlgebra.jl

The following functionality is provided for rational function fields.

### Greatest common divisor

```@docs
gcd(::Generic.Rat{T}, ::Generic.Rat{T}) where T <: FieldElement
```

**Examples**

```jldoctest
julia> R, x = RationalFunctionField(QQ, "x")
(Rational function field over Rationals, x)

julia> f = (x + 1)//(x^3 + 3x + 1)
(x + 1)//(x^3 + 3*x + 1)

julia> g = (x^2 + 2x + 1)//(x^2 + x + 1)
(x^2 + 2*x + 1)//(x^2 + x + 1)

julia> h = gcd(f, g)
(x + 1)//(x^5 + x^4 + 4*x^3 + 4*x^2 + 4*x + 1)

```

### Square root

```@docs
issquare(::Generic.Rat{T}) where T <: FieldElem
```

```@docs
Base.sqrt(::Generic.Rat{T}) where T <: FieldElem
```

**Examples**

```jldoctest
julia> R, x = RationalFunctionField(QQ, "x")
(Rational function field over Rationals, x)

julia> a = (21//4*x^6 - 15*x^5 + 27//14*x^4 + 9//20*x^3 + 3//7*x + 9//10)//(x + 3)
(21//4*x^6 - 15*x^5 + 27//14*x^4 + 9//20*x^3 + 3//7*x + 9//10)//(x + 3)

julia> sqrt(a^2)
(21//4*x^6 - 15*x^5 + 27//14*x^4 + 9//20*x^3 + 3//7*x + 9//10)//(x + 3)

julia> issquare(a^2)
true
```
