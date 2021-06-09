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

## Generic rational function field type

Rational functions in $k(x)$ have type `Generic.Rat{T}` where `T` is the type of
elements of the coefficient field $k$. See the file
`src/generic/GenericTypes.jl` for details.

Parent objects corresponding to the rational function field $k$ have type
`Generic.RationalFunctionField{T}`.

## Abstract types

The rational function types belong to the abstract type `Field`
and the rational function field types belong to the abstract type
`FieldElem`.

## Rational function field constructors

In order to construct rational functions in AbstractAlgebra.jl, one can first
construct the function field itself. This is accomplished with the following
constructor.

```julia
RationalFunctionField(k::Field, s::AbstractString; cached::Bool = true)
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

julia> m = S(numerator(x + 1, false), numerator(x + 2, false))
(x + 1)//(x + 2)
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

# Univariate function fields

Univariate function fields in AbstractAlgebra are algebraic extensions $K/k(x)$
of a rational function field $k(x)$ over a field $k$.

These are implemented in a module implemented in
`src/generic/FunctionField.jl`.

## Generic function field types

Function field objects $K/k(x)$ in AbstractAlgebra have type
`Generic.FunctionField{T}` where `T` is the type of elements of the field `k`.

Corresponding function field elements have type
`Generic.FunctionFieldElement{T}`. See the file `src/generic/GenericTypes.jl`
for details.

## Abstract types

Function field types belong to the abstract type `Field`
and their elements to the abstract type `FieldElem`.

## Function field constructors

In order to construct function fields in AbstractAlgebra.jl, one first
constructs the rational function field they are an extension of, then supplies
a polynomial over this field to the following constructor:

```julia
FunctionField(p::Poly{Rat{T}}, s::AbstractString; cached::Bool=true) where T <: FieldElement
```

Given an irreducible polynomial `p` over a rational function field return a
tuple `(S, z)` consisting of the parent object of the function field defined by
that polynomial over $k(x)$ and the generator `z`. By default the parent object
`S` will depend only on `p` and `s` and will be cached. Setting the optional
argument `cached` to `false` will prevent the parent object `S` from being
cached.

Here are some examples of creating function fields and making use of the
resulting parent objects to coerce various elements into the function field.

**Examples**

```jldoctest
julia> R1, x1 = RationalFunctionField(QQ, "x1") # characteristic 0
(Rational function field over Rationals, x1)

julia> U1, z1 = R1["z1"]
(Univariate Polynomial Ring in z1 over Rational function field over Rationals, z1)

julia> f = (x1^2 + 1)//(x1 + 1)*z1^3 + 4*z1 + 1//(x1 + 1)
(x1^2 + 1)//(x1 + 1)*z1^3 + 4*z1 + 1//(x1 + 1)

julia> S1, y1 = FunctionField(f, "y1")
(Function Field over Rationals with defining polynomial (x1^2 + 1)*y1^3 + (4*x1 + 4)*y1 + 1, y1)

julia> a = S1()
0

julia> b = S1((x1 + 1)//(x1 + 2))
(x1 + 1)//(x1 + 2)

julia> c = S1(1//3)
1//3

julia> R2, x2 = RationalFunctionField(GF(23), "x1") # characteristic p
(Rational function field over Finite field F_23, x1)

julia> U2, z2 = R2["z2"]
(Univariate Polynomial Ring in z2 over Rational function field over Finite field F_23, z2)

julia> g = z2^2 + 3z2 + 1
z2^2 + 3*z2 + 1

julia> S2, y2 = FunctionField(g, "y2")
(Function Field over Finite field F_23 with defining polynomial y2^2 + 3*y2 + 1, y2)

julia> d = S2(R2(5))
5

julia> e = S2(y2)
y2
```

## Basic function field functionality

Function fields implement the full Ring and Field interfaces. We give some
examples of such functionality.


**Examples**

```jldoctest
julia> R, x = RationalFunctionField(GF(23), "x") # characteristic p
(Rational function field over Finite field F_23, x)

julia> U, z = R["z"]
(Univariate Polynomial Ring in z over Rational function field over Finite field F_23, z)

julia> g = z^2 + 3z + 1
z^2 + 3*z + 1

julia> S, y = FunctionField(g, "y")
(Function Field over Finite field F_23 with defining polynomial y^2 + 3*y + 1, y)

julia> f = (x + 1)*y + 1
(x + 1)*y + 1

julia> base_ring(f)
Rational function field over Finite field F_23

julia> f^2
(20*x^2 + 19*x + 22)*y + 22*x^2 + 21*x

julia> f*inv(f)
1
```

## Function field functionality provided by AbstractAlgebra.jl

The following functionality is provided for function fields.

### Basic manipulation

```@docs
base_field(::Generic.FunctionField)
```

```@docs
var(::Generic.FunctionField)
```

```@docs
characteristic(S::Generic.FunctionField)
```

```@docs
defining_polynomial(R::Generic.FunctionField)
```

```@docs
Base.numerator(::Generic.FunctionField{T}, ::Bool=true) where T <: FieldElement
```

```@docs
Base.numerator(::Generic.FunctionFieldElem{T}, ::Bool=true) where T <: FieldElement
```

```@docs
degree(::Generic.FunctionField)
```

```@docs
gen(::Generic.FunctionField{T}) where T <: FieldElement
```

```@docs
isgen(a::Generic.FunctionFieldElem)
```

```@docs
coeff(::Generic.FunctionFieldElem, ::Int)
num_coeff(::Generic.FunctionFieldElem, ::Int)
```

**Examples**

```jldoctest
julia> R, x = RationalFunctionField(QQ, "x")
(Rational function field over Rationals, x)

julia> U, z = R["z"]
(Univariate Polynomial Ring in z over Rational function field over Rationals, z)

julia> g = z^2 + 3*(x + 1)//(x + 2)*z + 1
z^2 + (3*x + 3)//(x + 2)*z + 1

julia> S, y = FunctionField(g, "y")
(Function Field over Rationals with defining polynomial (x + 2)*y^2 + (3*x + 3)*y + x + 2, y)

julia> base_field(S)
Rational function field over Rationals

julia> var(S)
:y

julia> characteristic(S)
0

julia> defining_polynomial(S)
z^2 + (3*x + 3)//(x + 2)*z + 1

julia> numerator(S)
(x + 2)*y^2 + (3*x + 3)*y + x + 2

julia> denominator(S)
x + 2

julia> a = (x + 1)//(x^2 + 1)*y + 3x + 2
((x + 1)*y + 3*x^3 + 2*x^2 + 3*x + 2)//(x^2 + 1)

julia> numerator(a, false)
(x + 1)*y + 3*x^3 + 2*x^2 + 3*x + 2

julia> denominator(a, false)
x^2 + 1

julia> degree(S)
2

julia> gen(S)
y

julia> isgen(y)
true

julia> coeff(a, 1)
(x + 1)//(x^2 + 1)

julia> num_coeff(a, 1)
x + 1
```

### Trace and norm

```@docs
norm(::Generic.FunctionFieldElem)
```

```jldoctest
julia> R, x = RationalFunctionField(QQ, "x")
(Rational function field over Rationals, x)

julia> U, z = R["z"]
(Univariate Polynomial Ring in z over Rational function field over Rationals, z)

julia> g = z^2 + 3*(x + 1)//(x + 2)*z + 1
z^2 + (3*x + 3)//(x + 2)*z + 1

julia> S, y = FunctionField(g, "y")
(Function Field over Rationals with defining polynomial (x + 2)*y^2 + (3*x + 3)*y + x + 2, y)

julia> f = (-3*x - 5//3)//(x - 2)*y + (x^3 + 1//9*x^2 + 5)//(x - 2)
((-3*x - 5//3)*y + x^3 + 1//9*x^2 + 5)//(x - 2)

julia> norm(f)
(x^7 + 20//9*x^6 + 766//81*x^5 + 2027//81*x^4 + 110//3*x^3 + 682//9*x^2 + 1060//9*x + 725//9)//(x^3 - 2*x^2 - 4*x + 8)

julia> tr(f)
(2*x^4 + 38//9*x^3 + 85//9*x^2 + 24*x + 25)//(x^2 - 4)
```

