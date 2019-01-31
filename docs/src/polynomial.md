```@meta
CurrentModule = AbstractAlgebra
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

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")
T, z = QQ["z"]

f = R()
g = R(123)
h = S(BigInt(1234))
k = S(x + 1)
m = T(z + 1)
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

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")

f = x^3 + 3x + 21
g = (x + 1)*y^2 + 2x + 1

h = zero(S)
k = one(R)
isone(k) == true
iszero(f) == false
n = length(g)
U = base_ring(S)
V = base_ring(y + 1)
v = var(S)
T = parent(y + 1)
g == deepcopy(g)
t = divexact(2g, 2)
```

For polynomials over a field, the Euclidean Ring interface is implemented.

**Examples**

```julia
R, x = PolynomialRing(QQ, "x")
S = ResidueRing(R, x^3 + 3x + 1)
T, y = PolynomialRing(S, "y")

f = (3*x^2 + x + 2)*y + x^2 + 1
g = (5*x^2 + 2*x + 1)*y^2 + 2x*y + x + 1
h = (3*x^3 + 2*x^2 + x + 7)*y^5 + 2x*y + 1

invmod(f, g)
mulmod(f, g, h)
powmod(f, 3, h)
h = mod(f, g)
q, r = divrem(f, g)
d = gcd(f*h, g*h)
k = gcdinv(f, h)
m = lcm(f, h)
flag, q = divides(g^2, g)
valuation(3g^3, g) == 3
val, q = remove(5g^3, g)
r, s, t = gcdx(g, h)
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
ismonomial(::PolyElem)
```

```@docs
isterm(::PolyElem)
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")
T, z = PolynomialRing(QQ, "z")
U = ResidueRing(ZZ, 17)
V, w = PolynomialRing(U, "w")

a = zero(S)
b = one(S)

c = BigInt(1)//2*z^2 + BigInt(1)//3
d = x*y^2 + (x + 1)*y + 3

f = lead(d)
y = gen(S)
g = isgen(w)
m = isunit(b)
n = degree(d)
r = modulus(w)
isterm(2y^2) == true
ismonomial(y^2) == true
```

### Truncation

```@docs
truncate(::PolyElem, ::Int)
```

```@docs
mullow{T <: RingElem}(::PolyElem{T}, ::PolyElem{T}, ::Int)
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")

f = x*y^2 + (x + 1)*y + 3
g = (x + 1)*y + (x^3 + 2x + 2)

h = truncate(f, 1)
k = mullow(f, g, 4)
```

### Reversal

```@docs
reverse(::PolyElem, ::Int)
reverse(::PolyElem)
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")

f = x*y^2 + (x + 1)*y + 3

g = reverse(f, 7)
h = reverse(f)
```

### Shifting

```@docs
shift_left(::PolyElem, ::Int)
```

```@docs
shift_right(::PolyElem, ::Int)
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")

f = x*y^2 + (x + 1)*y + 3

g = shift_left(f, 7)
h = shift_right(f, 2)
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

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")

f = x*y^2 + (x + 1)*y + 3
g = (x + 1)*y + (x^3 + 2x + 2)

h = pseudorem(f, g)
q, r = pseudodivrem(f, g)
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

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")
   
f = x*y^2 + (x + 1)*y + 3
g = (x + 1)*y + (x^3 + 2x + 2)
M = R[x + 1 2x; x - 3 2x - 1]

k = evaluate(f, 3)
m = evaluate(f, x^2 + 2x + 1)
n = compose(f, g)
p = subst(f, M)
q = f(M)
r = f(23)
```

### Derivative and integral

```@docs
derivative(::PolyElem)
```

```@docs
integral{T <: Union{ResElem, FieldElem}}(::PolyElem{T})
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")
T, z = PolynomialRing(QQ, "z")
U = ResidueRing(T, z^3 + 3z + 1)
V, w = PolynomialRing(U, "w")

f = x*y^2 + (x + 1)*y + 3
g = (z^2 + 2z + 1)*w^2 + (z + 1)*w - 2z + 4

h = derivative(f)
k = integral(g)   
```

### Resultant and discriminant

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

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")

f = 3x*y^2 + (x + 1)*y + 3
g = 6(x + 1)*y + (x^3 + 2x + 2)

h = resultant(f, g)
k = discriminant(f)
```

### Newton representation

```@docs
monomial_to_newton!{T <: RingElem}(::Array{T, 1}, ::Array{T, 1})
```

```@docs
newton_to_monomial!{T <: RingElem}(::Array{T, 1}, ::Array{T, 1})
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")

f = 3x*y^2 + (x + 1)*y + 3
g = deepcopy(f)
roots = [R(1), R(2), R(3)]

monomial_to_newton!(g.coeffs, roots)
newton_to_monomial!(g.coeffs, roots)
```

### Interpolation

```@docs
interpolate{T <: RingElem}(::PolyRing, ::Array{T, 1}, ::Array{T, 1})
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")

xs = [R(1), R(2), R(3), R(4)]
ys = [R(1), R(4), R(9), R(16)]

f = interpolate(S, xs, ys)
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

```julia
R, x = PolynomialRing(ZZ, "x")
S, y = PolynomialRing(R, "y")

f = chebyshev_t(20, y)
g = chebyshev_u(15, y)
```
