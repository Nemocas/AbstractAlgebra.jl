```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Generic univariate polynomials

AbstractAlgebra.jl provides a module, implemented in `src/generic/Poly.jl` for generic
polynomials over any commutative ring belonging to the AbstractAlgebra abstract type
hierarchy.

As well as implementing the Univariate Polynomial interface, and relevant parts of the
Euclidean Ring interface for polynomials over a field, there are many additional
generic algorithms implemented for such polynomial rings. We describe this generic
functionality below.

All of the generic functionality is part of a submodule of AbstractAlgebra called
`Generic`. This is exported by default so that it is not necessary to qualify the
function names with the submodule name.

## Types and parent objects

Polynomials implemented using the AbstractAlgebra generics have type `Generic.Poly{T}`
where `T` is the type of elements of the coefficient ring. Internally they consist of
a Julia array of coefficients and some additional fields for length and a parent object,
etc. See the file `src/generic/GenericTypes.jl` for details.

Parent objects of such polynomials have type `Generic.PolyRing{T}`.

The string representation of the variable of the polynomial ring and the
base/coefficient ring $R$ is stored in the parent object. 

The polynomial element types belong to the abstract type `AbstractAlgebra.PolyElem{T}`
and the polynomial ring types belong to the abstract type `AbstractAlgebra.PolyRing{T}`. This enables one to write generic functions that can accept any AbstractAlgebra
polynomial type.

Note that both the generic polynomial ring type `Generic.PolyRing{T}` and the abstract
type it belongs to, `AbstractAlgebra.PolyRing{T}` are both called `PolyRing`. The 
former is a (parameterised) concrete type for a polynomial ring over a given base ring
whose elements have type `T`. The latter is an abstract type representing all
polynomial ring types in AbstractAlgebra.jl, whether generic or very specialised (e.g.
supplied by a C library).

## Polynomial ring constructors

In order to construct polynomials in AbstractAlgebra.jl, one must first construct the
polynomial ring itself. This is accomplished with the following constructor.

```julia
PolynomialRing(R::AbstractAlgebra.Ring, s::AbstractString; cached::Bool = true)
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

Here are some examples of creating polynomial rings and making use of the
resulting parent objects to coerce various elements into the polynomial ring.

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> T, z = QQ["z"]
(Univariate Polynomial Ring in z over Rationals, z)

julia> f = R()
0

julia> g = R(123)
123

julia> h = S(BigInt(1234))
1234

julia> k = S(x + 1)
(x+1)

julia> m = T(z + 1)
z+1//1

```

All of the examples here are generic polynomial rings, but specialised implementations
of polynomial rings provided by external modules will also usually provide a
`PolynomialRing` constructor to allow creation of their polynomial rings.

## Basic ring functionality

Once a polynomial ring is constructed, there are various ways to construct
polynomials in that ring.

The easiest way is simply using the generator returned by the `PolynomialRing`
constructor and build up the polynomial using basic arithmetic, as described in
the Ring interface. 

The Julia language also has special syntax for the construction of polynomials in terms
of a generator, e.g. we can write `2x` instead of `2*x`.

The polynomial rings in AbstractAlgebra.jl implement the full Ring interface. Of course
the entire Univariate Polynomial Ring interface is also implemented.

We give some examples of such functionality.

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> f = x^3 + 3x + 21
x^3+3*x+21

julia> g = (x + 1)*y^2 + 2x + 1
(x+1)*y^2+(2*x+1)

julia> h = zero(S)
0

julia> k = one(R)
1

julia> isone(k)
true

julia> iszero(f)
false

julia> n = length(g)
3

julia> U = base_ring(S)
Univariate Polynomial Ring in x over Integers

julia> V = base_ring(y + 1)
Univariate Polynomial Ring in x over Integers

julia> v = var(S)
:y

julia> T = parent(y + 1)
Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers

julia> g == deepcopy(g)
true

julia> t = divexact(2g, 2)
(x+1)*y^2+(2*x+1)

```

For polynomials over a field, the Euclidean Ring interface is implemented.

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> S = ResidueRing(R, x^3 + 3x + 1)
Residue ring of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1

julia> T, y = PolynomialRing(S, "y")
(Univariate Polynomial Ring in y over Residue ring of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1, y)

julia> f = (3*x^2 + x + 2)*y + x^2 + 1
(3//1*x^2+x+2//1)*y+(x^2+1//1)

julia> g = (5*x^2 + 2*x + 1)*y^2 + 2x*y + x + 1
(5//1*x^2+2//1*x+1//1)*y^2+(2//1*x)*y+(x+1//1)

julia> h = (3*x^3 + 2*x^2 + x + 7)*y^5 + 2x*y + 1
(2//1*x^2-8//1*x+4//1)*y^5+(2//1*x)*y+1//1

julia> invmod(f, g)
(707//3530*x^2+2151//1765*x+123//3530)*y+(-178//1765*x^2-551//3530*x+698//1765)

julia> mulmod(f, g, h)
(-30//1*x^2-43//1*x-9//1)*y^3+(-7//1*x^2-23//1*x-7//1)*y^2+(4//1*x^2-10//1*x-3//1)*y+(x^2-2//1*x)

julia> powmod(f, 3, h)
(69//1*x^2+243//1*x+79//1)*y^3+(78//1*x^2+180//1*x+63//1)*y^2+(27//1*x^2+42//1*x+18//1)*y+(3//1*x^2+3//1*x+2//1)

julia> h = mod(f, g)
(3//1*x^2+x+2//1)*y+(x^2+1//1)

julia> q, r = divrem(f, g)
(0//1, (3//1*x^2+x+2//1)*y+(x^2+1//1))

julia> d = gcd(f*h, g*h)
y+(1//11*x^2+6//11)

julia> k = gcdinv(f, h)
(y+(1//11*x^2+6//11), 0//1)

julia> m = lcm(f, h)
(-14//1*x^2-23//1*x-2//1)*y+(-4//1*x^2-5//1*x+1//1)

julia> flag, q = divides(g^2, g)
(true, (5//1*x^2+2//1*x+1//1)*y^2+(2//1*x)*y+(x+1//1))

julia> valuation(3g^3, g) == 3
true

julia> val, q = remove(5g^3, g)
(3, 5//1)

julia> r, s, t = gcdx(g, h)
(1//1, (311//3530*x^2-2419//3530*x+947//1765), (707//3530*x^2+2151//1765*x+123//3530)*y+(-178//1765*x^2-551//3530*x+698//1765))

```

Functions in the Euclidean Ring interface are supported over residue rings that are
not fields, except that if an impossible inverse is encountered during the computation
an error is thrown.

## Polynomial functionality provided by AbstractAlgebra.jl

The functionality listed below is automatically provided by AbstractAlgebra.jl for
any polynomial module that implements the full Univariate Polynomial Ring interface.
This includes AbstractAlgebra.jl's own generic polynomial rings.

But if a C library provides all the functionality documented in the Univariate
Polynomial Ring interface, then all the functions described here will also be 
automatically supplied by AbstractAlgebra.jl for that polynomial type.

Of course, modules are free to provide specific implementations of the functions
described here, that override the generic implementation.

### Basic functionality

```@docs
modulus{T <: ResElem}(::PolyElem{T})
```

```@docs
lead(::PolyElem)
trail(::PolyElem)
```

```@docs
gen(::PolyElem)
```

```@docs
isgen(::PolyElem)
```

```@docs
isunit(::PolyElem)
```

```@docs
degree(::Generic.PolynomialElem)
```

```@docs
ismonomial(::PolyElem)
```

```@docs
isterm(::PolyElem)
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

julia> a = zero(S)
0

julia> b = one(S)
1

julia> c = BigInt(1)//2*z^2 + BigInt(1)//3
1//2*z^2+1//3

julia> d = x*y^2 + (x + 1)*y + 3
(x)*y^2+(x+1)*y+3

julia> f = lead(d)
x

julia> y = gen(S)
y

julia> g = isgen(w)
true

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
(x)*y^2+(x+1)*y+3

julia> g = (x + 1)*y + (x^3 + 2x + 2)
(x+1)*y+(x^3+2*x+2)

julia> h = truncate(f, 1)
3

julia> k = mullow(f, g, 4)
(x^2+x)*y^3+(x^4+3*x^2+4*x+1)*y^2+(x^4+x^3+2*x^2+7*x+5)*y+(3*x^3+6*x+6)

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
(x)*y^2+(x+1)*y+3

julia> g = reverse(f, 7)
3*y^6+(x+1)*y^5+(x)*y^4

julia> h = reverse(f)
3*y^2+(x+1)*y+(x)

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
(x)*y^2+(x+1)*y+3

julia> g = shift_left(f, 7)
(x)*y^9+(x+1)*y^8+3*y^7

julia> h = shift_right(f, 2)
(x)

```

### Change of base ring

```@docs
change_base_ring(::PolyElem{T}, ::Any) where T <: RingElement
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
g = x^3+6*x + 1
change_base_ring(g, GF(2))
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
(x)*y^2+(x+1)*y+3

julia> g = (x + 1)*y + (x^3 + 2x + 2)
(x+1)*y+(x^3+2*x+2)

julia> h = pseudorem(f, g)
(x^7+3*x^5+2*x^4+x^3+5*x^2+4*x+1)

julia> q, r = pseudodivrem(f, g)
((x^2+x)*y+(-x^4-x^2+1), (x^7+3*x^5+2*x^4+x^3+5*x^2+4*x+1))

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
(x)*y^2+(x+1)*y+3

julia> g = (x + 1)*y + (x^3 + 2x + 2)
(x+1)*y+(x^3+2*x+2)

julia> M = R[x + 1 2x; x - 3 2x - 1]
[x+1 2*x]
[x-3 2*x-1]

julia> k = evaluate(f, 3)
12*x+6

julia> m = evaluate(f, x^2 + 2x + 1)
x^5+4*x^4+7*x^3+7*x^2+4*x+4

julia> n = compose(f, g)
(x^3+2*x^2+x)*y^2+(2*x^5+2*x^4+4*x^3+9*x^2+6*x+1)*y+(x^7+4*x^5+5*x^4+5*x^3+10*x^2+8*x+5)

julia> p = subst(f, M)
[3*x^3-3*x^2+3*x+4 6*x^3+2*x^2+2*x]
[3*x^3-8*x^2-2*x-3 6*x^3-8*x^2+2*x+2]

julia> q = f(M)
[3*x^3-3*x^2+3*x+4 6*x^3+2*x^2+2*x]
[3*x^3-8*x^2-2*x-3 6*x^3-8*x^2+2*x+2]

julia> r = f(23)
552*x+26

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
Residue ring of Univariate Polynomial Ring in z over Rationals modulo z^3+3//1*z+1//1

julia> V, w = PolynomialRing(U, "w")
(Univariate Polynomial Ring in w over Residue ring of Univariate Polynomial Ring in z over Rationals modulo z^3+3//1*z+1//1, w)

julia> f = x*y^2 + (x + 1)*y + 3
(x)*y^2+(x+1)*y+3

julia> g = (z^2 + 2z + 1)*w^2 + (z + 1)*w - 2z + 4
(z^2+2//1*z+1//1)*w^2+(z+1//1)*w+(-2//1*z+4//1)

julia> h = derivative(f)
(2*x)*y+(x+1)

julia> k = integral(g)
(1//3*z^2+2//3*z+1//3)*w^3+(1//2*z+1//2)*w^2+(-2//1*z+4//1)*w

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
(3*x)*y^2+(x+1)*y+3

julia> g = 6(x + 1)*y + (x^3 + 2x + 2)
(6*x+6)*y+(x^3+2*x+2)

julia> S = sylvester_matrix(f, g)
[3*x x+1 3]
[6*x+6 x^3+2*x+2 0]
[0 6*x+6 x^3+2*x+2]

julia> h = resultant(f, g)
3*x^7+6*x^5-6*x^3+96*x^2+192*x+96

julia> k = discriminant(f)
x^2-34*x+1

```

### Newton representation

```@docs
monomial_to_newton!{T <: RingElem}(::Array{T, 1}, ::Array{T, 1})
```

```@docs
newton_to_monomial!{T <: RingElem}(::Array{T, 1}, ::Array{T, 1})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> f = 3x*y^2 + (x + 1)*y + 3
(3*x)*y^2+(x+1)*y+3

julia> g = deepcopy(f)
(3*x)*y^2+(x+1)*y+3

julia> roots = [R(1), R(2), R(3)]
3-element Array{AbstractAlgebra.Generic.Poly{BigInt},1}:
 1
 2
 3

julia> monomial_to_newton!(g.coeffs, roots)

julia> newton_to_monomial!(g.coeffs, roots)

```

### Interpolation

```@docs
interpolate{T <: RingElem}(::PolyRing, ::Array{T, 1}, ::Array{T, 1})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integers, y)

julia> xs = [R(1), R(2), R(3), R(4)]
4-element Array{AbstractAlgebra.Generic.Poly{BigInt},1}:
 1
 2
 3
 4

julia> ys = [R(1), R(4), R(9), R(16)]
4-element Array{AbstractAlgebra.Generic.Poly{BigInt},1}:
 1
 4
 9
 16

julia> f = interpolate(S, xs, ys)
y^2

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
524288*y^20-2621440*y^18+5570560*y^16-6553600*y^14+4659200*y^12-2050048*y^10+549120*y^8-84480*y^6+6600*y^4-200*y^2+1

julia> g = chebyshev_u(15, y)
32768*y^15-114688*y^13+159744*y^11-112640*y^9+42240*y^7-8064*y^5+672*y^3-16*y

```
