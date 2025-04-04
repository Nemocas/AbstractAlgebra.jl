```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = AbstractAlgebra.doctestsetup()
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
`PolyRingElem{T}` and the polynomial ring types belong to the abstract
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

```@docs; canonical=false
polynomial_ring(R::NCRing, s::VarName; cached::Bool = true)
```

A shorthand version of this function is provided: given a base ring `R`, we abbreviate
the constructor as follows.

```julia
R[:x]
```

Here are some examples of creating polynomial rings and their associated
generators.

**Examples**

```jldoctest
julia> T, z = QQ[:z]
(Univariate polynomial ring in z over rationals, z)

julia> U, x = polynomial_ring(ZZ)
(Univariate polynomial ring in x over integers, x)
```

All of the examples here are generic polynomial rings, but specialised implementations
of polynomial rings provided by external modules will also usually provide a
`polynomial_ring` constructor to allow creation of their polynomial rings.

## Polynomial constructors

Once a polynomial ring is constructed, there are various ways to construct
polynomials in that ring.

The easiest way is simply using the generator returned by the `polynomial_ring`
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
polynomial(R::Ring, arr::Vector{T}, var::VarName=:x; cached::Bool=true)
```

Given an array of coefficients construct the polynomial with those coefficients
over the given ring and with the given variable.

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

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

julia> S, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over rationals, x)

julia> f = S(Rational{BigInt}[2, 3, 1])
x^2 + 3*x + 2

julia> g = S(BigInt[1, 0, 4])
4*x^2 + 1

julia> h = S([4, 7, 2, 9])
9*x^3 + 2*x^2 + 7*x + 4

julia> p = polynomial(ZZ, [1, 2, 3])
3*x^2 + 2*x + 1

julia> f = polynomial(ZZ, [1, 2, 3], :y)
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
similar(x::MyPoly{T}, R::Ring, var::VarName=var(parent(x))) where T <: RingElem
similar(x::MyPoly{T}, var::VarName=var(parent(x))) where T <: RingElem
zero(x::MyPoly{T}, R::Ring, var::VarName=var(parent(x))) where T <: RingElem
zero(x::MyPoly{T}, var::VarName=var(parent(x))) where T <: RingElem
```

Construct the zero polynomial with the given variable and coefficients in the
given ring, if specified, and in the coefficient ring of the given polynomial
otherwise.

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> f = 1 + 2x + 3x^2
3*x^2 + 2*x + 1

julia> g = similar(f)
0

julia> h = similar(f, QQ)
0

julia> k = similar(f, QQ, :y)
0
```

## Functions for types and parents of polynomial rings

```julia
base_ring(R::PolyRing)
base_ring(a::PolyRingElem)
```

Return the coefficient ring of the given polynomial ring or polynomial.

```julia
parent(a::NCRingElement)
```

Return the polynomial ring of the given polynomial..

```julia
characteristic(R::NCRing)
```

Return the characteristic of the given polynomial ring. If the characteristic
is not known, an exception is raised.

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

julia> U = base_ring(S)
Univariate polynomial ring in x over integers

julia> V = base_ring(y + 1)
Univariate polynomial ring in x over integers

julia> T = parent(y + 1)
Univariate polynomial ring in y over R
```

## Euclidean polynomial rings

For polynomials over a field, the [Euclidean Ring Interface](@ref) is implemented.

```julia
mod(f::PolyRingElem, g::PolyRingElem)
divrem(f::PolyRingElem, g::PolyRingElem)
div(f::PolyRingElem, g::PolyRingElem)
```

```julia
mulmod(f::PolyRingElem, g::PolyRingElem, m::PolyRingElem)
powermod(f::PolyRingElem, e::Int, m::PolyRingElem)
invmod(f::PolyRingElem, m::PolyRingElem)
```

```julia
divides(f::PolyRingElem, g::PolyRingElem)
remove(f::PolyRingElem, p::PolyRingElem)
valuation(f::PolyRingElem, p::PolyRingElem)
```

```julia
gcd(f::PolyRingElem, g::PolyRingElem)
lcm(f::PolyRingElem, g::PolyRingElem)
gcdx(f::PolyRingElem, g::PolyRingElem)
gcdinv(f::PolyRingElem, g::PolyRingElem)
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over rationals, x)

julia> S, = residue_ring(R, x^3 + 3x + 1);

julia> T, y = polynomial_ring(S, :y)
(Univariate polynomial ring in y over S, y)

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
iszero(a::PolyRingElem)
isone(a::PolyRingElem)
```

```julia
divexact(a::T, b::T) where T <: PolyRingElem
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
modulus{T <: ResElem}(::PolyRingElem{T})
```

```@docs
leading_coefficient(::PolyRingElem)
trailing_coefficient(::PolyRingElem)
constant_coefficient(::PolynomialElem)
```

```@docs
set_coefficient!(::PolynomialElem{T}, ::Int, c::T) where T <: RingElement
```

```@docs
tail(::PolynomialElem)
```

```@docs
gen(::PolyRingElem)
```

```@docs
is_gen(::PolyRingElem)
```

```@docs
is_monic(::PolyRingElem)
```

```@docs
is_square(::PolyRingElem)
```

```@docs
length(::PolynomialElem)
```

```@docs
degree(::PolynomialElem)
```

```@docs
is_monomial(::PolyRingElem)
```

```@docs
is_monomial_recursive(::PolyRingElem)
```

```@docs
is_term(::PolyRingElem)
```

```@docs
is_term_recursive(::PolyRingElem)
```

```@docs
is_constant(::PolynomialElem)
```

```@docs
is_separable(::PolyRingElem)
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

julia> T, z = polynomial_ring(QQ, :z)
(Univariate polynomial ring in z over rationals, z)

julia> U, = residue_ring(ZZ, 17);

julia> V, w = polynomial_ring(U, :w)
(Univariate polynomial ring in w over U, w)

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

julia> g = is_gen(w)
true

julia> divexact((2x + 1)*(x + 1), (x + 1))
2*x + 1

julia> m = is_unit(b)
true

julia> n = degree(d)
2

julia> r = modulus(w)
17

julia> is_term(2y^2)
true

julia> is_monomial(y^2)
true

julia> is_monomial_recursive(x*y^2)
true

julia> is_monomial(x*y^2)
false

julia> S, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

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
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

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
truncate(::PolyRingElem, ::Int)
```

```@docs
mullow{T <: RingElem}(::PolyRingElem{T}, ::PolyRingElem{T}, ::Int)
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

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
reverse(::PolyRingElem, ::Int)
reverse(::PolyRingElem)
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

julia> f = x*y^2 + (x + 1)*y + 3
x*y^2 + (x + 1)*y + 3

julia> g = reverse(f, 7)
3*y^6 + (x + 1)*y^5 + x*y^4

julia> h = reverse(f)
3*y^2 + (x + 1)*y + x

```

### Shifting

```@docs
shift_left(::PolyRingElem, ::Int)
```

```@docs
shift_right(::PolyRingElem, ::Int)
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

julia> f = x*y^2 + (x + 1)*y + 3
x*y^2 + (x + 1)*y + 3

julia> g = shift_left(f, 7)
x*y^9 + (x + 1)*y^8 + 3*y^7

julia> h = shift_right(f, 2)
x

```

### Inflation and deflation

```@docs
deflation(::PolyRingElem)
```

```@docs
inflate(::PolyRingElem, ::Int, ::Int)
inflate(::PolyRingElem, ::Int)
```

```@docs
deflate(::PolyRingElem, ::Int, ::Int)
deflate(::PolyRingElem, ::Int)
deflate(::PolyRingElem)
```

### Square root

```@docs
Base.sqrt(::PolyRingElem{T}; check::Bool) where T <: RingElement
```

**Examples**

```julia
R, x = polynomial_ring(ZZ, :x)
g = x^2+6*x+1
sqrt(g^2)
```

### Change of base ring

```@docs
change_base_ring(::Ring, ::PolyRingElem{T}) where T <: RingElement
change_coefficient_ring(::Ring, ::PolyRingElem{T}) where T <: RingElement
map_coefficients(::Any, ::PolyRingElem{<:RingElement})
```

**Examples**

```julia
R, x = polynomial_ring(ZZ, :x)
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
pseudorem{T <: RingElem}(::PolyRingElem{T}, ::PolyRingElem{T})
```

```@docs
pseudodivrem{T <: RingElem}(::PolyRingElem{T}, ::PolyRingElem{T})
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

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
content(::PolyRingElem)
```

```@docs
primpart(::PolyRingElem)
```

**Examples**

```
R, x = polynomial_ring(ZZ, :x)
S, y = polynomial_ring(R, :y)

k = x*y^2 + (x + 1)*y + 3

n = content(k)
p = primpart(k*(x^2 + 1))
```

### Evaluation, composition and substitution

```@docs
evaluate(::PolyRingElem, b)
```

```@docs
compose(::PolyRingElem, ::PolyRingElem)
```

```@docs
subst{T <: RingElem}(::PolyRingElem{T}, ::Any)
```

We also overload the functional notation so that the polynomial $f$ can be
evaluated at $a$ by writing $f(a)$.

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

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

julia> n = compose(f, g; inner = :second)
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
derivative(::PolyRingElem)
```

```@docs
integral{T <: Union{ResElem, FieldElem}}(::PolyRingElem{T})
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

julia> T, z = polynomial_ring(QQ, :z)
(Univariate polynomial ring in z over rationals, z)

julia> U, = residue_ring(T, z^3 + 3z + 1);

julia> V, w = polynomial_ring(U, :w)
(Univariate polynomial ring in w over U, w)

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
sylvester_matrix{T <: RingElem}(::PolyRingElem{T}, ::PolyRingElem{T})
```

```@docs
resultant{T <: RingElem}(::PolyRingElem{T}, ::PolyRingElem{T})
```

```@docs
resx{T <: RingElem}(::PolyRingElem{T}, ::PolyRingElem{T})
```

```@docs
discriminant(a::PolyRingElem)
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

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
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

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
roots(f::PolyRingElem)
roots(R::Field, f::PolyRingElem)
```

### Interpolation

```@docs
interpolate{T <: RingElem}(::PolyRing, ::Vector{T}, ::Vector{T})
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

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
polynomial_to_power_sums(::PolyRingElem{T}) where T <: RingElem
```

```@docs
power_sums_to_polynomial(::Vector{T}) where T <: RingElem
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

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
chebyshev_t(::Int, ::PolyRingElem)
```

```@docs
chebyshev_u(::Int, ::PolyRingElem)
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

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
rand(R::PolyRing, deg_range::AbstractUnitRange{Int}, v...)
rand(R::PolyRing, deg::Int, v...)
```

**Examples**

```julia
R, x = polynomial_ring(ZZ, :x)
f = rand(R, -1:3, -10:10)

S, y = polynomial_ring(GF(7), :y)
g = rand(S, 2:2)

U, z = polynomial_ring(R, :z)
h = rand(U, 3:3, -1:2, -10:10)
```

### Ring homomorphisms

```@docs
hom(::PolyRing, ::NCRing, ::Any, ::Any)
```
