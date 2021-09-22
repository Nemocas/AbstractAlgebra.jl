```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Univariate polynomial functionality

AbstractAlgebra.jl provides a module, implemented in `src/Poly.jl` for
polynomials over any commutative ring belonging to the AbstractAlgebra abstract type
hierarchy. This functionality will work for any univariate polynomial type which
follows the Univariate Polynomial Ring interface.

## Generic univariate polynomial types

AbstractAlgebra.jl provides a generic polynomial type based on Julia arrays which
is implemented in `src/generic/Poly.jl`.

These generic polynomials have type `Generic.Poly{T}` where `T` is the type of
elements of the coefficient ring. Internally they consist of a Julia array of
coefficients and some additional fields for length and a parent object, etc. See
the file `src/generic/GenericTypes.jl` for details.

Parent objects of such polynomials have type `Generic.PolyRing{T}`.

The string representation of the variable of the polynomial ring and the
base/coefficient ring $R$ is stored in the parent object.

## Abstract types

All univariate polynomial element types belong to the abstract type
`PolyElem{T}` and the polynomial ring types belong to the abstract
type `PolyRing{T}`. This enables one to write generic functions that
can accept any AbstractAlgebra polynomial type.

!!! note

    Both the generic polynomial ring type `Generic.PolyRing{T}` and the abstract
    type it belongs to, `PolyRing{T}`, are called `PolyRing`. The
    former is a (parameterised) concrete type for a polynomial ring over a given base ring
    whose elements have type `T`. The latter is an abstract type representing all
    polynomial ring types in AbstractAlgebra.jl, whether generic or very specialised (e.g.
    supplied by a C library).

## Polynomial ring constructors

In order to construct polynomials in AbstractAlgebra.jl, one must first construct the
polynomial ring itself. This is accomplished with the following constructor.

```julia
PolynomialRing(R::Ring, s::AbstractString; cached::Bool = true)
```

Given a base ring `R` and string `s` specifying how the generator (variable) should be
printed, return a tuple `S, x` representing the new polynomial ring $S = R[x]$ and the
generator $x$ of the ring. By default the parent object `S` will depend only on `R` and
`x` and will be cached. Setting the optional argument `cached` to `false` will prevent
the parent object `S` from being cached.

A shorthand version of this function is provided: given a base ring `R`, we abbreviate
the constructor as follows.

```julia
R["x"]
```

It is also possible to create a polynomial ring with default symbol as follows.
This is a lightweight constructor and should be used in generic algorithms
wherever possible when creating polynomial rings where the symbol does not
matter.

```julia
PolyRing(R::Ring)
```

Given a base ring `R` return the polynomial ring $S = R[x]$. Note that unlike
the constructors above, the return type is not a tuple. Only the ring is
returned and not the generator. The polynomial ring is not cached.

Here are some examples of creating polynomial rings and their associated
generators.

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> T, z = QQ["z"]
(Univariate Polynomial Ring in z over Rationals, z)

julia> U = PolyRing(ZZ)
Univariate Polynomial Ring in x over Integers
```

All of the examples here are generic polynomial rings, but specialised implementations
of polynomial rings provided by external modules will also usually provide a
`PolynomialRing` constructor to allow creation of their polynomial rings.

## Polynomial constructors

Once a polynomial ring is constructed, there are various ways to construct
polynomials in that ring.

The easiest way is simply using the generator returned by the `PolynomialRing`
constructor and build up the polynomial using basic arithmetic.

The Julia language has special syntax for the construction of polynomials in terms
of a generator, e.g. we can write `2x` instead of `2*x`.

A second way is to use the polynomial ring to construct a polynomial. There are
the usual ways of constructing an element of a ring.

```julia
(R::PolyRing)() # constructs zero
(R::PolyRing)(c::Integer)
(R::PolyRing)(c::elem_type(R))
(R::PolyRing{T})(a::T) where T <: RingElement
```

For polynommials there is also the following more general constructor accepting
an array of coefficients.

```julia
(S::PolyRing{T})(A::Vector{T}) where T <: RingElem
(S::PolyRing{T})(A::Vector{U}) where T <: RingElem, U <: RingElem
(S::PolyRing{T})(A::Vector{U}) where T <: RingElem, U <: Integer
```

Construct the polynomial in the ring `S` with the given array of coefficients,
i.e. where `A[1]` is the constant coefficient.

A third way of constructing polynomials is to construct them directly without
creating the polynomial ring.

```julia
polynomial(R::Ring, arr::Vector{T}, var::String="x"; cached::Bool=true)
```

Given an array of coefficients construct the polynomial with those coefficients
over the given ring and with the given variable.

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> f = x^3 + 3x + 21
x^3 + 3*x + 21

julia> g = (x + 1)*y^2 + 2x + 1
(x + 1)*y^2 + 2*x + 1

julia> R()
0

julia> S(1)
1

julia> S(y)
y

julia> S(x)
x

julia> S, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> f = S(Rational{BigInt}[2, 3, 1])
x^2 + 3*x + 2

julia> g = S(BigInt[1, 0, 4])
4*x^2 + 1

julia> h = S([4, 7, 2, 9])
9*x^3 + 2*x^2 + 7*x + 4

julia> p = polynomial(ZZ, [1, 2, 3])
3*x^2 + 2*x + 1

julia> f = polynomial(ZZ, [1, 2, 3], "y")
3*y^2 + 2*y + 1
```

## Similar and zero

Another way of constructing polynomials is to construct one similar to an
existing polynomial using either similar or zero. 

```julia
similar(x::MyPoly{T}, R::Ring=base_ring(x)) where T <: RingElem
zero(x::MyPoly{T}, R::Ring=base_ring(x)) where T <: RingElem
```

Construct the zero polynomial with the same variable as the given polynomial
with coefficients in the given ring. Both functions behave the same way for
polynomials.

```julia
similar(x::MyPoly{T}, R::Ring, var::String=String(var(parent(x)))) where T <: RingElem
similar(x::MyPoly{T}, var::String=String(var(parent(x)))) where T <: RingElem
zero(x::MyPoly{T}, R::Ring, var::String=String(var(parent(x)))) where T <: RingElem
zero(x::MyPoly{T}, var::String=String(var(parent(x)))) where T <: RingElem
```

Construct the zero polynomial with the given variable and coefficients in the
given ring, if specified, and in the coefficient ring of the given polynomial
otherwise.

**Examples**

```@jldoctest
julia> f = 1 + 2x + 3x^2
3*x^2 + 2*x + 1

julia> g = similar(f)
0

julia> h = similar(f, QQ)
0

julia> k = similar(f, QQ, "y")
0
```

## Functions for types and parents of polynomial rings

```julia
base_ring(R::PolyRing)
base_ring(a::PolyElem)
```

Return the coefficient ring of the given polynomial ring or polynomial.

```julia
parent(a::NCRingElement)
```

Return the polynomial ring of the given polynomial..

```julia
characteristic(R::NCRing)
```

Return the characteristic of the given polynomial ring.

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> U = base_ring(S)
Univariate Polynomial Ring in x over Integers

julia> V = base_ring(y + 1)
Univariate Polynomial Ring in x over Integers

julia> T = parent(y + 1)
Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers
```

## Euclidean polynomial rings

For polynomials over a field, the Euclidean Ring interface is implemented.

```julia
mod(f::PolyElem, g::PolyElem)
divrem(f::PolyElem, g::PolyElem)
div(f::PolyElem, g::PolyElem)
```

```julia
mulmod(f::PolyElem, g::PolyElem, m::PolyElem)
powermod(f::PolyElem, e::Int, m::PolyElem)
invmod(f::PolyElem, m::PolyElem)
```

```julia
divides(f::PolyElem, g::PolyElem)
remove(f::PolyElem, p::PolyElem)
valuation(f::PolyElem, p::PolyElem)
```

```julia
gcd(f::PolyElem, g::PolyElem)
lcm(f::PolyElem, g::PolyElem)
gcdx(f::PolyElem, g::PolyElem)
gcdinv(f::PolyElem, g::PolyElem)
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> S = ResidueRing(R, x^3 + 3x + 1)
Residue ring of Univariate Polynomial Ring in x over Rationals modulo x^3 + 3*x + 1

julia> T, y = PolynomialRing(S, "y")
(Univariate Polynomial Ring in y over Residue ring of Univariate Polynomial Ring in x over Rationals modulo x^3 + 3*x + 1, y)

julia> f = (3*x^2 + x + 2)*y + x^2 + 1
(3*x^2 + x + 2)*y + x^2 + 1

julia> g = (5*x^2 + 2*x + 1)*y^2 + 2x*y + x + 1
(5*x^2 + 2*x + 1)*y^2 + 2*x*y + x + 1

julia> h = (3*x^3 + 2*x^2 + x + 7)*y^5 + 2x*y + 1
(2*x^2 - 8*x + 4)*y^5 + 2*x*y + 1

julia> invmod(f, g)
(707//3530*x^2 + 2151//1765*x + 123//3530)*y - 178//1765*x^2 - 551//3530*x + 698//1765

julia> mulmod(f, g, h)
(-30*x^2 - 43*x - 9)*y^3 + (-7*x^2 - 23*x - 7)*y^2 + (4*x^2 - 10*x - 3)*y + x^2 - 2*x

julia> powermod(f, 3, h)
(69*x^2 + 243*x + 79)*y^3 + (78*x^2 + 180*x + 63)*y^2 + (27*x^2 + 42*x + 18)*y + 3*x^2 + 3*x + 2

julia> h = mod(f, g)
(3*x^2 + x + 2)*y + x^2 + 1

julia> q, r = divrem(f, g)
(0, (3*x^2 + x + 2)*y + x^2 + 1)

julia> div(g, f)
(-5//11*x^2 + 2//11*x + 6//11)*y - 13//121*x^2 - 3//11*x - 78//121

julia> d = gcd(f*h, g*h)
y + 1//11*x^2 + 6//11

julia> k = gcdinv(f, h)
(y + 1//11*x^2 + 6//11, 0)

julia> m = lcm(f, h)
(-14*x^2 - 23*x - 2)*y - 4*x^2 - 5*x + 1

julia> flag, q = divides(g^2, g)
(true, (5*x^2 + 2*x + 1)*y^2 + 2*x*y + x + 1)

julia> valuation(3g^3, g) == 3
true

julia> val, q = remove(5g^3, g)
(3, 5)

julia> r, s, t = gcdx(g, h)
(1, 311//3530*x^2 - 2419//3530*x + 947//1765, (707//3530*x^2 + 2151//1765*x + 123//3530)*y - 178//1765*x^2 - 551//3530*x + 698//1765)

```

Functions in the Euclidean Ring interface are supported over residue rings that are
not fields, except that if an impossible inverse is encountered during the computation
an error is thrown.

## Polynomial functions

### Basic functionality

All basic ring functionality is provided for polynomials. The most important
such functions are the following.

```julia
zero(R::PolyRing)
one(R::PolyRing)
iszero(a::PolyElem)
isone(a::PolyElem)
```

```julia
divexact(a::T, b::T) where T <: PolyElem
```

All functions in the polynomial interface are provided. The most important
are the following.

```julia
var(S::PolyRing)
symbols(S::PolyRing{T}) where T <: RingElem
```

Return a symbol or length 1 array of symbols, respectively, specifying the
variable of the polynomial ring. This symbol is converted to a string when
printing polynomials in that ring.

In addition, the following basic functions are provided.

```@docs
modulus{T <: ResElem}(::PolyElem{T})
```

```@docs
leading_coefficient(::PolyElem)
trailing_coefficient(::PolyElem)
constant_coefficient(::PolynomialElem)
```

```@docs
set_coefficient!(::PolynomialElem{T}, ::Int, c::T) where T <: RingElement
```

```@docs
tail(::PolynomialElem)
```

```@docs
gen(::PolyElem)
```

```@docs
isgen(::PolyElem)
```

```@docs
ismonic(::PolyElem)
```

```@docs
issquare(::PolyElem)
```

```@docs
degree(::PolynomialElem)
```

```@docs
ismonomial(::PolyElem)
```

```@docs
ismonomial_recursive(::PolyElem)
```

```@docs
isterm(::PolyElem)
```

```@docs
isterm_recursive(::PolyElem)
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> T, z = PolynomialRing(QQ, "z")
(Univariate Polynomial Ring in z over Rationals, z)

julia> U = ResidueRing(ZZ, 17)
Residue ring of Integers modulo 17

julia> V, w = PolynomialRing(U, "w")
(Univariate Polynomial Ring in w over Residue ring of Integers modulo 17, w)

julia> var(R)
:x

julia> symbols(R)
1-element Vector{Symbol}:
 :x

julia> a = zero(S)
0

julia> b = one(S)
1

julia> isone(b)
true

julia> c = BigInt(1)//2*z^2 + BigInt(1)//3
1//2*z^2 + 1//3

julia> d = x*y^2 + (x + 1)*y + 3
x*y^2 + (x + 1)*y + 3

julia> f = leading_coefficient(d)
x

julia> y = gen(S)
y

julia> g = isgen(w)
true

julia> divexact((2x + 1)*(x + 1), (x + 1))
2*x + 1

julia> m = isunit(b)
true

julia> n = degree(d)
2

julia> r = modulus(w)
17

julia> isterm(2y^2)
true

julia> ismonomial(y^2)
true

julia> ismonomial_recursive(x*y^2)
true

julia> ismonomial(x*y^2)
false

julia> S, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> f = x^3 + 3x + 1
x^3 + 3*x + 1

julia> g = S(BigInt[1, 2, 0, 1, 0, 0, 0]);

julia> n = length(f)
4

julia> c = coeff(f, 1)
3

julia> g = set_coefficient!(g, 2, ZZ(11))
x^3 + 11*x^2 + 2*x + 1

julia> g = set_coefficient!(g, 7, ZZ(4))
4*x^7 + x^3 + 11*x^2 + 2*x + 1
```

### Iterators

An iterator is provided to return the coefficients of a univariate polynomial.
The iterator is called `coefficients` and allows iteration over the
coefficients, starting with the term of degree zero (if there is one). Note
that coefficients of each degree are given, even if they are zero. This is best
illustrated by example.

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> f = x^2 + 2
x^2 + 2

julia> C = collect(coefficients(f))
3-element Vector{BigInt}:
 2
 0
 1

julia> for c in coefficients(f)
          println(c)
       end
2
0
1
```

### Truncation

```@docs
truncate(::PolyElem, ::Int)
```

```@docs
mullow{T <: RingElem}(::PolyElem{T}, ::PolyElem{T}, ::Int)
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> f = x*y^2 + (x + 1)*y + 3
x*y^2 + (x + 1)*y + 3

julia> g = (x + 1)*y + (x^3 + 2x + 2)
(x + 1)*y + x^3 + 2*x + 2

julia> h = truncate(f, 1)
3

julia> k = mullow(f, g, 4)
(x^2 + x)*y^3 + (x^4 + 3*x^2 + 4*x + 1)*y^2 + (x^4 + x^3 + 2*x^2 + 7*x + 5)*y + 3*x^3 + 6*x + 6

```

### Reversal

```@docs
reverse(::PolyElem, ::Int)
reverse(::PolyElem)
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> f = x*y^2 + (x + 1)*y + 3
x*y^2 + (x + 1)*y + 3

julia> g = reverse(f, 7)
3*y^6 + (x + 1)*y^5 + x*y^4

julia> h = reverse(f)
3*y^2 + (x + 1)*y + x

```

### Shifting

```@docs
shift_left(::PolyElem, ::Int)
```

```@docs
shift_right(::PolyElem, ::Int)
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> f = x*y^2 + (x + 1)*y + 3
x*y^2 + (x + 1)*y + 3

julia> g = shift_left(f, 7)
x*y^9 + (x + 1)*y^8 + 3*y^7

julia> h = shift_right(f, 2)
x

```

### Inflation and deflation

```@docs
deflation(::PolyElem)
```

```@docs
inflate(::PolyElem, ::Int, ::Int)
inflate(::PolyElem, ::Int)
```

```@docs
deflate(::PolyElem, ::Int, ::Int)
deflate(::PolyElem, ::Int)
deflate(::PolyElem)
```

### Square root

```@docs
Base.sqrt(::PolyElem{T}; check::Bool) where T <: RingElement
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
g = x^2+6*x+1
sqrt(g^2)
```

### Change of base ring

```@docs
change_base_ring(::Ring, ::PolyElem{T}) where T <: RingElement
change_coefficient_ring(::Ring, ::PolyElem{T}) where T <: RingElement
map_coefficients(::Any, ::PolyElem{<:RingElement})
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
g = x^3+6*x + 1
change_base_ring(GF(2), g)
change_coefficient_ring(GF(2), g)
```

### Pseudodivision

Given two polynomials $a, b$, pseudodivision computes polynomials $q$ and $r$
with length$(r) <$ length$(b)$ such that
$$L^d a = bq + r,$$
where $d =$ length$(a) -$ length$(b) + 1$ and $L$ is the leading coefficient
of $b$.

We call $q$ the pseudoquotient and $r$ the pseudoremainder.

```@docs
pseudorem{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})
```

```@docs
pseudodivrem{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> f = x*y^2 + (x + 1)*y + 3
x*y^2 + (x + 1)*y + 3

julia> g = (x + 1)*y + (x^3 + 2x + 2)
(x + 1)*y + x^3 + 2*x + 2

julia> h = pseudorem(f, g)
x^7 + 3*x^5 + 2*x^4 + x^3 + 5*x^2 + 4*x + 1

julia> q, r = pseudodivrem(f, g)
((x^2 + x)*y - x^4 - x^2 + 1, x^7 + 3*x^5 + 2*x^4 + x^3 + 5*x^2 + 4*x + 1)

```

### Content and primitive part

```@docs
content(::PolyElem)
```

```@docs
primpart(::PolyElem)
```

**Examples**

```
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")

k = x*y^2 + (x + 1)*y + 3

n = content(k)
p = primpart(k*(x^2 + 1))
```

### Evaluation, composition and substitution

```@docs
evaluate{T <: RingElem}(::PolyElem{T}, ::T)
evaluate(::PolyElem, ::Integer)
```

```@docs
compose(::PolyElem, ::PolyElem)
```

```@docs
subst{T <: RingElem}(::PolyElem{T}, ::Any)
```

We also overload the functional notation so that the polynomial $f$ can be
evaluated at $a$ by writing $f(a)$.

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)


julia> f = x*y^2 + (x + 1)*y + 3
x*y^2 + (x + 1)*y + 3

julia> g = (x + 1)*y + (x^3 + 2x + 2)
(x + 1)*y + x^3 + 2*x + 2

julia> M = R[x + 1 2x; x - 3 2x - 1]
[x + 1       2*x]
[x - 3   2*x - 1]

julia> k = evaluate(f, 3)
12*x + 6

julia> m = evaluate(f, x^2 + 2x + 1)
x^5 + 4*x^4 + 7*x^3 + 7*x^2 + 4*x + 4

julia> n = compose(f, g)
(x^3 + 2*x^2 + x)*y^2 + (2*x^5 + 2*x^4 + 4*x^3 + 9*x^2 + 6*x + 1)*y + x^7 + 4*x^5 + 5*x^4 + 5*x^3 + 10*x^2 + 8*x + 5

julia> p = subst(f, M)
[3*x^3 - 3*x^2 + 3*x + 4       6*x^3 + 2*x^2 + 2*x]
[3*x^3 - 8*x^2 - 2*x - 3   6*x^3 - 8*x^2 + 2*x + 2]

julia> q = f(M)
[3*x^3 - 3*x^2 + 3*x + 4       6*x^3 + 2*x^2 + 2*x]
[3*x^3 - 8*x^2 - 2*x - 3   6*x^3 - 8*x^2 + 2*x + 2]

julia> r = f(23)
552*x + 26

```

### Derivative and integral

```@docs
derivative(::PolyElem)
```

```@docs
integral{T <: Union{ResElem, FieldElem}}(::PolyElem{T})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> T, z = PolynomialRing(QQ, "z")
(Univariate Polynomial Ring in z over Rationals, z)

julia> U = ResidueRing(T, z^3 + 3z + 1)
Residue ring of Univariate Polynomial Ring in z over Rationals modulo z^3 + 3*z + 1

julia> V, w = PolynomialRing(U, "w")
(Univariate Polynomial Ring in w over Residue ring of Univariate Polynomial Ring in z over Rationals modulo z^3 + 3*z + 1, w)

julia> f = x*y^2 + (x + 1)*y + 3
x*y^2 + (x + 1)*y + 3

julia> g = (z^2 + 2z + 1)*w^2 + (z + 1)*w - 2z + 4
(z^2 + 2*z + 1)*w^2 + (z + 1)*w - 2*z + 4

julia> h = derivative(f)
2*x*y + x + 1

julia> k = integral(g)
(1//3*z^2 + 2//3*z + 1//3)*w^3 + (1//2*z + 1//2)*w^2 + (-2*z + 4)*w

```

### Resultant and discriminant

```@docs
sylvester_matrix{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})
```

```@docs
resultant{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})
```

```@docs
resx{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})
```

```@docs
discriminant(a::PolyElem)
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> f = 3x*y^2 + (x + 1)*y + 3
3*x*y^2 + (x + 1)*y + 3

julia> g = 6(x + 1)*y + (x^3 + 2x + 2)
(6*x + 6)*y + x^3 + 2*x + 2

julia> S = sylvester_matrix(f, g)
[    3*x           x + 1               3]
[6*x + 6   x^3 + 2*x + 2               0]
[      0         6*x + 6   x^3 + 2*x + 2]

julia> h = resultant(f, g)
3*x^7 + 6*x^5 - 6*x^3 + 96*x^2 + 192*x + 96

julia> k = discriminant(f)
x^2 - 34*x + 1

```

### Newton representation

```@docs
monomial_to_newton!{T <: RingElem}(::Vector{T}, ::Vector{T})
```

```@docs
newton_to_monomial!{T <: RingElem}(::Vector{T}, ::Vector{T})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> f = 3x*y^2 + (x + 1)*y + 3
3*x*y^2 + (x + 1)*y + 3

julia> g = deepcopy(f)
3*x*y^2 + (x + 1)*y + 3

julia> roots = [R(1), R(2), R(3)]
3-element Vector{AbstractAlgebra.Generic.Poly{BigInt}}:
 1
 2
 3

julia> monomial_to_newton!(g.coeffs, roots)

julia> newton_to_monomial!(g.coeffs, roots)
```

### Roots

```@docs
roots(f::PolyElem)
roots(f::PolyElem, R::Field)
```

### Interpolation

```@docs
interpolate{T <: RingElem}(::PolyRing, ::Vector{T}, ::Vector{T})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> xs = [R(1), R(2), R(3), R(4)]
4-element Vector{AbstractAlgebra.Generic.Poly{BigInt}}:
 1
 2
 3
 4

julia> ys = [R(1), R(4), R(9), R(16)]
4-element Vector{AbstractAlgebra.Generic.Poly{BigInt}}:
 1
 4
 9
 16

julia> f = interpolate(S, xs, ys)
y^2

```

### Power sums

```@docs
polynomial_to_power_sums(::PolyElem{T}) where T <: RingElem
```

```@docs
power_sums_to_polynomial(::Vector{T}) where T <: RingElem
```

**Examples**

```@jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> f = x^4 - 2*x^3 + 10*x^2 + 7*x - 5
x^4 - 2*x^3 + 10*x^2 + 7*x - 5

julia> V = polynomial_to_power_sums(f)
4-element Vector{BigInt}:
   2
 -16
 -73
  20

julia> power_sums_to_polynomial(V)
x^4 - 2*x^3 + 10*x^2 + 7*x - 5
```

### Special functions

The following special functions can be computed for any polynomial ring.
Typically one uses the generator $x$ of a polynomial ring to get the respective
special polynomials expressed in terms of that generator.

```@docs
chebyshev_t(::Int, ::PolyElem)
```

```@docs
chebyshev_u(::Int, ::PolyElem)
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> f = chebyshev_t(20, y)
524288*y^20 - 2621440*y^18 + 5570560*y^16 - 6553600*y^14 + 4659200*y^12 - 2050048*y^10 + 549120*y^8 - 84480*y^6 + 6600*y^4 - 200*y^2 + 1

julia> g = chebyshev_u(15, y)
32768*y^15 - 114688*y^13 + 159744*y^11 - 112640*y^9 + 42240*y^7 - 8064*y^5 + 672*y^3 - 16*y

```

### Random generation

One may generate random polynomials with degrees in a given range. Additional
parameters are used to construct coefficients as elements of the coefficient
ring.

```julia
rand(R::PolyRing, deg_range::UnitRange{Int}, v...)
rand(R::PolyRing, deg::Int, v...)
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
f = rand(R, -1:3, -10:10)

S, y = PolynomialRing(GF(7), "y")
g = rand(S, 2:2)

U, z = PolynomialRing(R, "z")
h = rand(U, 3:3, -1:2, -10:10)
```


