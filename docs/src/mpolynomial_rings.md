```@meta
CurrentModule = AbstractAlgebra
```

# Multivariate Polynomial Ring Interface

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

Note that the type `T` must (transitively) belong to the abstract type `RingElem` or
more generally the union type `RingElement` which includes the Julia integer, rational
and floating point types.

### Constructors

To construct a multivariate polynomial ring, there is the following constructor.

```julia
PolynomialRing(R::Ring, s::Vector{AbstractString}; ordering=:lex, cached=true)
```

Returns a tuple, `S, vars` consisting of a polynomial ring $S$ and an array of
generators (variables) which print according to the strings in the supplied
vector $s$. The ordering can at present be `:lex`, `:deglex` or `:degrevlex`.
By default, the polynomial ring is cached, and creating a polynomial ring with
the same data will return the same ring object $S$. If this caching is not
desired, it can be switched off by setting `cached=false`.

Polynomials in a given ring can be constructed using the generators and basic
polynomial arithmetic. However, this is inefficient and the following build
context is provided for building polynomials term-by-term. It assumes the
polynomial data type is random access, and so the constructor functions must
be reimplemented for all other types of polynomials.

```julia
MPolyBuildCtx(R::MPolyRing)
```

Return a build context for creating polynomials in the given polynomial ring.

```julia
push_term!(M::MPolyBuildCtx, c::RingElem, v::Vector{Int})
```

Add the term with coefficient $c$ and exponent vector $v$ to the polynomial
under construction in the build context $M$.

```julia
finish(M::MPolyBuildCtx)
```

Finish construction of the polynomial, sort the terms, remove duplicate and
zero terms and return the created polynomial.

**Examples**

```julia
R, (x, y) = PolynomialRing(ZZ, ["x", "y"])

M = MPolyBuildCtx(R)
push_term!(C, ZZ(3), [1, 2])
push_term!(C, ZZ(2), [1, 1])
push_term!(C, ZZ(4), [0, 0])
p = finish(C)
```

### Data type and parent object methods

```julia
symbols(S::MyMPolyRing{T}) where T <: AbstractAlgebra.RingElem
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

Return an array of all the generators (variables) of the given polynomial ring
(as polynomials).

The first entry in the array will be the variable with most significance with
respect to the ordering.

```julia
gen(S::MyMPolyRing{T}, i::Int) where T <: AbstractAlgebra.RingElem
```

Return the $i$-th generator (variable) of the given polynomial ring (as a
polynomial).

```julia
ordering(S::MyMPolyRing{T})
```

Return the ordering of the given polynomial ring as a symbol. Supported values currently
include `:lex`, `:deglex` and `:degrevlex`.

**Examples**

```julia
S, (x, y) = PolynomialRing(QQ, ["x", "y"]; ordering=:deglex)

V = symbols(S)
X = gens(S)
ord = ordering(S)
```

### Basic manipulation of rings and elements

```julia
length(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return the number of nonzero terms of the given polynomial. The length of the zero
polynomial is defined to be $0$. The return value should be of type `Int`.

Return an array of the degrees of the polynomial $f$ in each of the variables.

```julia
total_degree(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return the total degree of the polynomial $f$, i.e. the highest sum of
exponents occuring in any term of $f$.

```julia
isgen(x::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return `true` if $x$ is a generator of the polynomial ring.

```julia
coeffs(p::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return an iterator for the coefficients of the polynomial $p$, starting
with the coefficient of the most significant term with respect to the
ordering. Generic code will provide this function automatically for
random access polynomials that implement the `coeff` function.

```julia
monomials(p::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return an iterator for the monomials of the polynomial $p$, starting with
the monomial of the most significant term with respect to the ordering.
Monomials in AbstractAlgebra are defined to have coefficient $1$. See the
function `terms` if you also require the coefficients, however note that
only monomials can be compared. Generic code will provide this function
automatically for random access polynomials that implement the `monomial`
function.

```julia
terms(p::MyMPoly{T}) where T <: AbstractAlgebra.RingElem
```

Return an iterator for the terms of the polynomial $p$, starting with
the most significant term with respect to the ordering. Terms in
AbstractAlgebra include the coefficient. Generic code will provide this
function automatically for random access polynomials that implement the
`term` function.

```julia
exponent_vectors(a::MyMPoly{T}) where T <: RingElement
```

Return an iterator for the exponent vectors for each of the terms of the
polynomial starting with the most significant term with respect to the
ordering. Each exponent vector is an array of `Int`s, one for each
variable, in the order given when the polynomial ring was created.
Generic code will provide this function automatically for random access
polynomials that implement the `exponent_vector` function.

**Examples**

```julia
S, (x, y) = PolynomialRing(ZZ, ["x", "y"])

f = x^3*y + 3x*y^2 + 1

n = length(f)
isgen(y) == true
B, b = max_degrees(f)
nvars(f) == 2
C = collect(coeffs(f))
M = collect(monomials(f))
T = collect(terms(f))
V = collect(exponent_vectors(f))
d = total_degree(f)
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

### Ad hoc exact division

For any ring that implements exact division, the following can be implemented.

```julia
divexact(f::MyMPoly{T}, c::Integer) where T <: AbstractAlgebra.RingElem
divexact(f::MyMPoly{T}, c::Rational) where T <: AbstractAlgebra.RingElem
divexact(f::MyMPoly{T}, c::T) where T <: AbstractAlgebra.RingElem
```

Divide the polynomial exactly by the constant $c$.

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

f = 3x^2*y^2 + 2x + 1

f1 = divexact(f, 5)
f2 = divexact(f, QQ(2, 3))
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

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

f = 2x^2*y + 2x + y + 1
g = x + y

q = div(f, g)
q, r = divrem(f, g)
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
(S::MyMPolyRing{T})(A::Vector{T}, m::Vector{Vector{Int}}) where T <: AbstractAlgebra.RingElem
```

Create the polynomial in the given ring with nonzero coefficients specified by
the elements of $A$ and corresponding exponent vectors given by the elements of
$m$.

There is no assumption about coefficients being nonzero or terms being in order
or unique. Zero terms are removed by the function, duplicate terms are combined
(added) and the terms are sorted so that they are in the correct order.

Each exponent vector uses a separate integer for each exponent field, the first
of which should be the exponent for the most significant variable with respect
to the ordering. All exponents must be non-negative.

A library may also optionally provide an interface that makes use of `BigInt`
(or any other big integer type) for exponents instead of `Int`.

**Examples**

```julia
S, (x, y) = PolynomialRing(QQ, ["x", "y"])

f = S(Rational{BigInt}[2, 3, 1], [[3, 2], [1, 0], [0, 1]])
```

### Sparse distributed, random access basic manipulation

```julia
coeff(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElem
```

Return the coefficient of the $n$-th term of $f$. The first term should be the most
significant term with respect to the ordering.

```julia
coeff(a::MyMPoly{T}, exps::Vector{Int}) where T <: RingElement
```

Return the coefficient of the term with the given exponent vector, or zero
if there is no such term.

```julia
monomial(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElem
monomial!(m::MyMPoly{T}, f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElem
```

Return the $n$-th monomial of $f$ or set $m$ to the $n$-th monomial of $f$,
respectively. The first monomial should be the most significant term with
respect to the ordering. Monomials have coefficient $1$ in AbstractAlgebra.
See the function `term` if you also require the coefficient, however, note
that only monomials can be compared.

```julia
term(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElem
```

Return the $n$-th term of $f$. The first term should be the one whose
monomial is most significant with respect to the ordering.

```julia
exponent(f::MyMPoly{T}, i::Int, j::Int) where T <: AbstractAlgebra.RingElem
```

Return the exponent of the $j$-th variable in the $n$-th term of the polynomial
$f$. The first term is the one with whose monomial is most significant with
respect to the ordering.

```julia
exponent_vector(a::MyMPoly{T}, i::Int) where T <: RingElement
```

Return a vector of exponents, corresponding to the exponent vector of the
i-th term of the polynomial. Term numbering begins at $1$ and the exponents
are given in the order of the variables for the ring, as supplied when the
ring was created.

```julia
setcoeff!(a::MyMPoly, exps::Vector{Int}, c::S) where S <: RingElement
```

Set the coefficient of the term with the given exponent vector to the given
value $c$. If no such term exists (and $c \neq 0$), one will be inserted. This
function takes $O(\log n)$ operations if a term with the given exponent already
exists and $c \neq 0$, or if the term is inserted at the end of the polynomial.
Otherwise it can take $O(n)$ operations in the worst case. This function must
return the modified polynomial.

**Examples**

```julia
S, (x, y) = PolynomialRing(ZZ, ["x", "y"])

f = x^3*y + 3x*y^2 + 1

c1 = coeff(f, 1)
c2 = coeff(f, x^3*y)
m = monomial(f, 2)
e1 = exponent(f, 1, 1)
v1 = exponent_vector(f, 1)
t1 = term(f, 1)
setcoeff!(f, [3, 1], 12)
```

### Unsafe functions

The following functions must be provided, but are considered unsafe, as they
may leave the polynomials in an inconsistent state and they mutate their
inputs. As usual, such functions should only be applied on polynomials that
have no references elsewhere in the system and are mainly intended to be used
in carefully written library code, rather than by users.

Users should instead build polynomials using the constructors described above.

```julia
fit!(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElem
```

Ensure that the polynomial $f$ internally has space for $n$ nonzero terms. This
function must mutate the function in-place if it is mutable. It does not return
the mutated polynomial. Immutable types can still be supported by defining this
function to do nothing.

```julia
setcoeff!(a::MyMPoly{T}, i::Int, c::T) where T <: RingElement
setcoeff!(a::MyMPoly{T}, i::Int, c::U) where {T <: RingElement, U <: Integer}
```

Set the $i$-th coefficient of the polynomial $a$ to $c$. No check is performed
on the index $i$ or for $c = 0$. It may be necessary to call
`combine_like_terms` after calls to this function, to remove zero terms. The
function must return the modified polynomial.

```julia
combine_like_terms!(a::MyMPoly{T}) where T <: RingElement
```

Remove zero terms and combine any adjacent terms with the same exponent
vector (by adding them). It is assumed that all the exponent vectors are
already in the correct order with respect to the ordering. The function
must return the resulting polynomial.

```julia
set_exponent_vector!(a::MyMPoly{T}, i::Int, exps::Vector{Int}) where T <: RingElement 
```

Set the $i$-th exponent vector to the given exponent vector. No check is
performed on the index $i$, which is assumed to be valid (or that the
polynomial has enough space allocated). No sorting of exponents is performed
by this function. To sort the terms after setting any number of exponents
with this function, run the `sort_terms!` function. The function must return
the modified polynomial.

```julia
sort_terms!(a::MyMPoly{T}) where {T <: RingElement}
```

Sort the terms of the given polynomial according to the polynomial ring
ordering. Zero terms and duplicate exponents are ignored. To deal with those
call `combine_like_terms`. The sorted polynomial must be returned by the
function.

## Optional functionality for multivariate polynomials

The following functions can optionally be implemented for multivariate
polynomial types.

### Reduction by an ideal

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

Q, r = divrem(f, [g, h])
```

### Evaluation

```julia
evaluate(a::MyMPoly{T}, A::Vector{T}) where T <: RingElem
```

Evaluate the polynomial at the given values in the coefficient ring of the
polynomial. The result should be an element of the coefficient ring.

```julia
evaluate(f::MyMPoly{T}, A::Vector{U}) where {T <: AbstractAlgebra.RingElem, U <: Intege
r}
```

Evaluate the polynomial $f$ at the values specified by the entries of the array $A$.

```julia
(a::MyMPoly{T})(vals::Union{NCRingElem, RingElement}...) where T <: RingElement
```

Evaluate the polynomial at the given arguments. This provides functional
notation for polynomial evaluation, i.e. $f(a, b, c)$. It must be defined
for each supported polynomial type (Julia does not allow functional
notation to be defined for an abstract type).

The code for this function in MPoly.jl can be used when implementing this
as it provides the most general possible evaluation, which is much more
general than the case of evaluation at elements of the same ring.

The evaluation should succeed for any set of values for which a
multiplication is defined with the product of a coefficient and all the
values before it.

Note that the values may be in non-commutative rings. Products are
performed in the order of the variables in the polynomial ring that the
polynomial belongs to, preceded by a multiplication by the coefficient
on the left.

**Examples**

```julia
R, (x, y) = AbstractAlgebra.PolynomialRing(ZZ, ["x", "y"])
S = MatrixAlgebra(ZZ, 2)

f = x*y + x + y + 1

evaluate(f, [ZZ(1), ZZ(2)])
evaluate(f, [2, 3])
f(1, 2)
f(ZZ(1), ZZ(2))
f(x - y, x + y)

M1 = S([1 2; 3 4])
M2 = S([2 4; 1 -1])
M3 = S([1 -1; 1 1])

f(M1, M2)
f(M1, ZZ(2))
```

### Derivations

The following function allows to compute derivations of multivariate
polynomials of type MPoly.

```julia
derivative(f::MyMPoly{T}, j::Int) where T <: AbstractAlgebra.RingElem
```

Compute the derivative of $f$ with respect to the $j$-th variable of the
polynomial ring.

**Examples**

```julia
R, (x, y) = AbstractAlgebra.PolynomialRing(ZZ, ["x", "y"])

f = x*y + x + y + 1

derivative(f, 1)
derivative(f, 2)
```
