```@meta
CurrentModule = AbstractAlgebra
```

# Multvariate Polynomial Ring Interface

Multivariate polynomial rings are supported in AbstractAlgebra.jl, and in addition to
the standard Ring interface, numerous additional functions are provided.

Unlike other kinds of rings, even complex operations such as GCD depend heavily on the
multivariate representation. Therefore AbstractAlgebra.jl cannot provide much in the
way of additional functionality to external multivariate implementations.

This means that external libraries must be able to implement their multivariate
formats in whatever way they see fit. The required interface here should be implemented,
even if it is not optimal. But it can be extended, either by implementing one of the
optional interfaces, or by extending the required interface in some other way.

Naturally, any multivariate polynomial ring implementation provides the full Ring
interface, in order to be treated as a ring for the sake of AbstractAlgebra.jl.

Considerations which make it impossible for AbstractAlgebra.jl to provide generic
functionality on top of an arbitrary multivariate module include:

  * orderings (lexical, degree, weighted, block, arbitrary)
  * sparse or dense representation
  * distributed or recursive representation
  * packed or unpacked exponents
  * exponent bounds (and whether adaptive or not)
  * random access or iterators
  * whether monomials and polynomials have the same type
  * whether special cache aware data structures such as Geobuckets are used

## Types and parents

AbstractAlgebra.jl provides two abstract types for multivariate polynomial rings and
their elements:

  * `MPolyRing{T}` is the abstract type for multivariate polynomial ring parent types
  * `MPolyElem{T}` is the abstract type for multivariate polynomial types

We have that `MPolyRing{T} <: AbstractAlgebra.Ring` and 
`MPolyElem{T} <: AbstractAlgebra.RingElem`.

Note that both abstract types are parameterised. The type `T` should usually be the type
of elements of the coefficient ring of the polynomial ring. For example, in the case of
$\mathbb{Z}[x, y]$ the type `T` would be the type of an integer, e.g. `BigInt`.

Multivariate polynomial rings should be made unique on the system by caching parent
objects (unless an optional `cache` parameter is set to `false`). Multivariate
polynomial rings should at least be distinguished based on their base (coefficient)
ring and number of variables. But if they have the same base ring, symbols (for their
variables/generators) and ordering, they should certainly have the same parent object.

See `src/generic/GenericTypes.jl` for an example of how to implement such a cache (which
usually makes use of a dictionary).

## Required functionality for multivariate polynomials

In addition to the required functionality for the Ring interface, the Multivariate
Polynomial interface has the following required functions.

We suppose that `R` is a fictitious base ring (coefficient ring) and that `S` is a
multivariate polynomial ring over `R` (i.e. $S = R[x, y, \ldots]$) with parent object
`S` of type `MyMPolyRing{T}`. We also assume the polynomials in the ring have type
`MyMPoly{T}`, where `T` is the type of elements of the base (coefficient) ring.

Of course, in practice these types may not be parameterised, but we use parameterised
types here to make the interface clearer.

Note that the type `T` must (transitively) belong to the abstract type `RingElem`.

### Data type and parent object methods

```julia
vars(S::MyMPolyRing{T}) where T <: AbstractAlgebra.RingElem
```

Return an array of `Symbol`s representing the variables (generators) of the polynomial
ring. Note that these are `Symbol`s not `String`s, though their string values will
usually be used when printing polynomials.

```julia
nvars(f::MyMPolyRing{T}) where T <: AbstractAlgebra.RingElem
```

Return the number of variables of the polynomial ring.


```julia
gens(S::MyMPolyRing{T}) where T <: AbstractAlgebra.RingElem
```

Return an array of all the generators (variables) of the given polynomial ring (as
polynomials).

The first entry in the array will be the variable with most significance with respect
to the ordering.

```julia
ordering(S::MyMPolyRing{T})
```

Return the ordering of the given polynomial ring as a symbol. Supported values currently
include `:lex`, `:deglex` and `:degrevlex`.

**Examples**

```julia
S, (x, y) = PolynomialRing(QQ, ["x", "y"]; ordering=:deglex)

V = vars(S)
X = gens(S)
ord = ordering(S)
```

### Basic manipulation of rings and elements

```julia
length(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return the number of nonzero terms of the given polynomial. The length of the zero
polynomial is defined to be $0$. The return value should be of type `Int`.

```julia
isgen(x::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return `true` if $x$ is a generator of the polynomial ring.

```julia
max_degrees(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Returns a tuple `(B, b)` consisting of an array of `Int`s specifying the highest power
of each variable that appears in the given polynomial and `b` the largest of the values
in `B`.

```julia
isunit(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return `true` if $f$ is a unit in its parent polynomial ring.

```julia
isconstant(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return `true` if $f$ is a constant polynomial. The zero polynomial is considered
constant for the purposes of this function.

```julia
isterm(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return `true` if $f$ consists of a single term.

```julia
ismonomial(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return `true` if $f$ consists of a single term with coefficient $1$.

**Examples**

```julia
S, (x, y) = PolynomialRing(ZZ, ["x", "y"])

f = x^3*y + 3x*y^2 + 1

n = length(f)
isgen(y) == true
B, b = max_degrees(f)
nvars(f) == 2
isunit(f) == false
isconstant(f) == false
isterm(2x*y) == true
ismonomial(x*y) == false
```

### Exact division

For any ring that implements exact division, the following can be implemented.

```julia
divexact(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return the exact quotient of $f$ by $g$ if it exists, otherwise throw an error.

```julia
divides(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return a tuple `(flag, q)` where `flag` is `true` if $g$ divides $f$, in which case
$q$ will be the exact quotient, or `flag` is false and $q$ is set to zero.

```julia
remove(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Returns a tuple $(v, q)$ such that the highest power of $g$ that divides $f$ is $g^v$
and the cofactor is $q$.

```julia
valuation(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Returns $v$ such that the highest power of $g$ that divides $f$ is $g^v$.

**Examples**

```julia
R, (x, y) = PolynomialRing(ZZ, ["x", "y"])

f = 2x^2*y + 2x + y + 1
g = x^2*y^2 + 1

flag, q = divides(f*g, f)
d = divexact(f*g, f)
v, q = remove(f*g^3, g)
n = valuation(f*g^3, g)
```

### Euclidean division

Although multivariate polynomial rings are not in general Euclidean, it is possible to
define a quotient with remainder function that depends on the polynomial ordering in
the case that the quotient ring is a field or a Euclidean domain. In the case that
a polynomial $g$ divides a polynomial $f$, the result no longer depends on the ordering
and the remainder is zero, with the quotient agreeing with the exact quotient.

```julia
divrem(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return a tuple $(q, r)$ such that $f = qg + r$, where the coefficients of terms of
$r$ whose monomials are divisible by the leading monomial of $g$ are reduced modulo the
leading coefficient of $g$ (according to the Euclidean function on the coefficients).

Note that the result of this function depends on the ordering of the polynomial ring.

```julia
div(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

As per the `divrem` function, but returning the quotient only. Especially when the
quotient happens to be exact, this function can be exceedingly fast.

```julia
divrem(f::MyMPoly{T}, G::Array{MyMPoly{T}, 1}) where T <: AbstractAlgebra.RingElem
```

As per the `divrem` function above, except that each term of $r$ starting with the
most significant term, is reduced modulo the leading terms of each of the polynomials
in the array $G$ for which the leading monomial is a divisor.

A tuple $(Q, r)$ is returned from the function, where $Q$ is an array of polynomials
of the same length as $G$, and such that $f = r + \sum Q[i]G[i]$.

The result is again dependent on the ordering in general, but if the polynomials in $G$
are over a field and the reduced generators of a Groebner basis, then the result is
unique.

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

f = 2x^2*y + 2x + y + 1
g = x + y
h = y + 1

q = div(f, g)
q, r = divrem(f, g)
Q, r = divrem(f, [g, h])
```

### Evaluation

```julia
evaluate(f::MyMPoly{T}, A::Array{T, 1}) where T <: AbstractAlgebra.RingElem
```

Evaluate the polynomial $f$ at the values specified by the entries of the array $A$.


```julia
evaluate(f::MPoly{T}, A::Array{T, 1}) where T <: Integer
```

Evaluate the polynomial $f$ at the values specified by the entries of the array $A$.

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

f = 2x^2*y + 2x + y + 1

m = evaluate(f, Rational{BigInt}[2, 3])
n = evaluate(f, [2, 3])
```

### GCD

In cases where there is a meaningful Euclidean structure on the coefficient ring, it is
possible to compute the GCD of multivariate polynomials.

```julia
gcd(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return a greatest common divisor of $f$ and $g$.

**Examples**

```julia
R, (x, y) = PolynomialRing(ZZ, ["x", "y"])

f = 2x^2*y + 2x + y + 1
g = x^2*y^2 + 1

d = gcd(f*g^2, f^2*g)
```

## Interface for sparse distributed, random access multivariates

The following additional functions should be implemented by libraries that provide a
sparse distributed polynomial format, stored in a representation for which terms can
be accessed in constant time (e.g. where arrays are used to store coefficients and
exponent vectors).

### Sparse distributed, random access constructors

In addition to the standard constructors, the following constructor, taking arrays of
coefficients and exponent vectors, should be provided.

```julia
(S::MyMPolyRing{T})(A::Array{T, 1}, m::Array{UInt, 2}) where T <: AbstractAlgebra.RingEle
m
```

Create the polynomial in the given ring with nonzero coefficients specified by the
elements of `A` and corresponding exponent vectors given by the elements of `m`. For
efficiency reason, the exponents of term $i$ are given by the vector `m[:, i]` since
Julia uses column major two dimensional arrays.

For maximum compatibility with external libraries, the coefficient (and term) at index
$1$ correspond to the most significant term with respect to the polynomial ring
ordering.

Each exponent vector uses a separate word for each exponent field, the first of which
should be any degree or weight, and otherwise should be the exponent for the most
significant variable with respect to the ordering. The top bit of each word is reserved
to detect overflows.

If a full word is not used for exponents, a check should be done to ensure there are
no overflows before setting the exponents.

A library may also optionally provide an interface that makes use of `BigInt` (or any
other big integer type) for exponents instead of `UInt`.

**Examples**

```julia
S, (x, y) = PolynomialRing(QQ, ["x", "y"])

f = S(Rational{BigInt}[2, 3, 1], UInt[3 2 1; 0 1 0])
```

### Sparse distributed, random access basic manipulation

```julia
coeff(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElem
```

Return the coefficient of the $(n+1)$-th term of $f$. The first term should be the most
significant term with respect to the ordering.

```julia
exponent(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElem
```

Return an array of `Int`s giving the vector of exponents for the $n + 1$-th term of $f$.
The first entry of the array should correspond to the exponent of the most significant
variable with respect to the ordering.

```julia
exponent!(A::Array{Int, 1}, f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElem
```

As per `exponent`, but set the values in the array `A` rather than allocating an array
for this purpose. The array is also returned by the function after being mutated.

```julia
fit!(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElem
```

Ensure that the polynomial $f$ internally has space for $n$ nonzero terms. This function
must mutate the function in-place if it is mutable. It does not return the mutated
polynomial. Immutable types can still be supported by defining this function to do
nothing.

**Examples**

```julia
S, (x, y) = PolynomialRing(ZZ, ["x", "y"])

f = x^3*y + 3x*y^2 + 1

c = coeff(f, 1)
fit!(f, 8)
```

### Derivations

The following function allows to compute derivations of multivariate polynomials of type MPoly.

```@docs
derivative{T <: AbstractAlgebra.RingElem}(f::AbstractAlgebra.Generic.MPoly{T}, x::AbstractAlgebra.Generic.MPoly{T})
```

**Example**

```julia
R,(x,y) = AbstractAlgebra.PolynomialRing(ZZ,["x","y"])
f = x*y + x + y + 1
derivative(f,x)
derivative(f,y)
```

