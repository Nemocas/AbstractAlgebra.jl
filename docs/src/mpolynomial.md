```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Sparse distributed multivariate polynomials

AbstractAlgebra.jl provides a module, implemented in `src/MPoly.jl` for
sparse distributed multivariate polynomials over any commutative ring belonging to the
AbstractAlgebra abstract type hierarchy.

## Generic sparse distributed multivariable polynomial types

AbstractAlgebra provides a generic multivariate polynomial type `Generic.MPoly{T}`
where `T` is the type of elements of the coefficient ring.

The polynomials are implemented using a Julia array of coefficients and a 2-dimensional
Julia array of `UInt`s for the exponent vectors. Note that exponent $n$ is represented
by the $n$-th column of the exponent array, not the $n$-th row. This is because Julia
uses a column major representation. See the file `src/generic/GenericTypes.jl` for
details.

The top bit of each `UInt` is reserved for overflow detection.

Parent objects of such polynomials have type `Generic.MPolyRing{T}`.

The string representation of the variables of the polynomial ring and the
base/coefficient ring $R$ and the ordering are stored in the parent object.

## Abstract types

The polynomial element types belong to the abstract type `MPolyRingElem{T}`
and the polynomial ring types belong to the abstract type `MPolyRing{T}`.

!!! note

    Note that both the generic polynomial ring type `Generic.MPolyRing{T}` and the abstract
    type it belongs to, `MPolyRing{T}` are both called `MPolyRing`. The
    former is a (parameterised) concrete type for a polynomial ring over a given base ring
    whose elements have type `T`. The latter is an abstract type representing all
    multivariate polynomial ring types in AbstractAlgebra.jl, whether generic or very
    specialised (e.g. supplied by a C library).

## Polynomial ring constructors

In order to construct multivariate polynomials in AbstractAlgebra.jl, one must first
construct the polynomial ring itself. This is accomplished with the following
constructors.

```@docs
polynomial_ring(::Ring, ::Vector{Symbol})
polynomial_ring(::Ring, ::Vararg)
polynomial_ring(::Ring, ::Int)
@polynomial_ring
```

Like for univariate polynomials, a shorthand constructor is
provided when the number of generators is greater than `1`: given a base ring
`R`, we abbreviate the constructor as follows:

```julia
R["x", "y", ...]
```

Here are some examples of creating multivariate polynomial rings and making use
of the resulting parent objects to coerce various elements into the polynomial
ring.

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"]; ordering=:deglex)
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> T, (z, t) = QQ["z", "t"]
(Multivariate polynomial ring in 2 variables over rationals, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}[z, t])

julia> f = R()
0

julia> g = R(123)
123

julia> h = R(BigInt(1234))
1234

julia> k = R(x + 1)
x + 1

julia> m = R(x + y + 1)
x + y + 1

julia> derivative(k, 1)
1

julia> derivative(k, 2)
0

julia> R, x = polynomial_ring(ZZ, 10); R
Multivariate polynomial ring in 10 variables x1, x2, x3, x4, ..., x10
  over integers

```

## Polynomial constructors

Multivariate polynomials can be constructed from the generators in the
usual way using arithmetic operations.

Also, all of the standard ring element constructors may be used to construct
multivariate polynomials.

```julia
(R::MPolyRing{T})() where T <: RingElement
(R::MPolyRing{T})(c::Integer) where T <: RingElement
(R::MPolyRing{T})(a::elem_type(R)) where T <: RingElement
(R::MPolyRing{T})(a::T) where T <: RingElement
```

For more efficient construction of multivariate polynomial, one can use the
`MPoly` build context, where terms (coefficient followed by an exponent vector)
are pushed onto a context one at a time and then the polynomial constructed
from those terms in one go using the `finish` function.

```@docs
MPolyBuildCtx(R::MPolyRing)
push_term!(M::MPolyBuildCtx, c::RingElem, v::Vector{Int})
finish(M::MPolyBuildCtx)
```

Note that the `finish` function resets the build context so that it can be
used to construct multiple polynomials..

When a multivariate polynomial type has a representation that allows constant
time access (e.g. it is represented internally by arrays), the following
additional constructor is available. It takes and array of coefficients and
and array of exponent vectors.

```julia
(S::MPolyRing{T})(A::Vector{T}, m::Vector{Vector{Int}}) where T <: RingElem
```

Create the polynomial in the given ring with nonzero coefficients specified by
the elements of $A$ and corresponding exponent vectors given by the elements of
$m$.

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> C = MPolyBuildCtx(R)
Builder for an element of Multivariate polynomial ring in 2 variables over integers

julia> push_term!(C, ZZ(3), [1, 2]);


julia> push_term!(C, ZZ(2), [1, 1]);


julia> push_term!(C, ZZ(4), [0, 0]);


julia> f = finish(C)
3*x*y^2 + 2*x*y + 4

julia> push_term!(C, ZZ(4), [1, 1]);


julia> f = finish(C)
4*x*y

julia> S, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over rationals, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}[x, y])

julia> f = S(Rational{BigInt}[2, 3, 1], [[3, 2], [1, 0], [0, 1]])
2*x^3*y^2 + 3*x + y
```

## Functions for types and parents of multivariate polynomial rings

```julia
base_ring(R::MPolyRing)
base_ring(a::MPolyRingElem)
```

Return the coefficient ring of the given polynomial ring or polynomial,
respectively.


```julia
parent(a::MPolyRingElem)
```

Return the polynomial ring of the given polynomial.

```julia
characteristic(R::MPolyRing)
```

Return the characteristic of the given polynomial ring. If the characteristic
is not known, an exception is raised.

## Polynomial functions

### Basic manipulation

All the standard ring functions are available, including the following.

```julia
zero(R::MPolyRing)
one(R::MPolyRing)
iszero(a::MPolyRingElem)
isone(a::MPolyRingElem)
```

```
divexact(a::T, b::T) where T <: MPolyRingElem
```

All basic functions from the Multivariate Polynomial interface are provided.

```julia
symbols(S::MPolyRing)
number_of_variables(f::MPolyRing)
gens(S::MPolyRing)
gen(S::MPolyRing, i::Int)
```

```julia
ordering(S::MPolyRing{T})
```

Note that the currently supported orderings are `:lex`, `:deglex` and
`:degrevlex`.

```julia
length(f::MPolyRingElem)
degrees(f::MPolyRingElem)
total_degree(f::MPolyRingElem)
```

```julia
is_gen(x::MPolyRingElem)
```

```julia
divexact(f::T, g::T) where T <: MPolyRingElem
```

For multivariate polynomial types that allow constant time access to
coefficients, the following are also available, allowing access to the given
coefficient, monomial or term. Terms are numbered from the most significant
first.

```julia
coeff(f::MPolyRingElem, n::Int)
coeff(a::MPolyRingElem, exps::Vector{Int})
```

Access a coefficient by term number or exponent vector.

```julia
monomial(f::MPolyRingElem, n::Int)
monomial!(m::T, f::T, n::Int) where T <: MPolyRingElem
```

The second version writes the result into a preexisting polynomial
object to save an allocation.

```julia
term(f::MPolyRingElem, n::Int)
```

```julia
exponent(f::MyMPolyRingElem, i::Int, j::Int)
```

Return the exponent of the $j$-th variable in the $i$-th term of the polynomial
$f$.

```julia
exponent_vector(a::MPolyRingElem, i::Int)
```

```julia
setcoeff!(a::MPolyRingElem{T}, exps::Vector{Int}, c::T) where T <: RingElement
```

Although multivariate polynomial rings are not usually Euclidean, the following
functions from the Euclidean interface are often provided.

```julia
divides(f::T, g::T) where T <: MPolyRingElem
remove(f::T, g::T) where T <: MPolyRingElem
valuation(f::T, g::T) where T <: MPolyRingElem
```

```julia
divrem(f::T, g::T) where T <: MPolyRingElem
div(f::T, g::T) where T <: MPolyRingElem
```

Compute a tuple $(q, r)$ such that $f = qg + r$, where the coefficients of terms of
$r$ whose monomials are divisible by the leading monomial of $g$ are reduced modulo the
leading coefficient of $g$ (according to the Euclidean function on the coefficients).
The `divrem` version returns both quotient and remainder whilst the `div` version only
returns the quotient.

Note that the result of these functions depend on the ordering of the polynomial ring.

```julia
gcd(f::T, g::T) where T <: MPolyRingElem
```

The following functionality is also provided for all multivariate polynomials.

```@docs
is_univariate(::MPolyRing{T}) where T <: RingElement
```

```@docs
vars(p::MPolyRingElem{T}) where T <: RingElement
```

```@docs
var_index(::MPolyRingElem{T}) where T <: RingElement
```

```@docs
degree(::MPolyRingElem{T}, ::Int) where T <: RingElement
```

```@docs
degree(::MPolyRingElem{T}, ::MPolyRingElem{T}) where T <: RingElement
```

```@docs
degrees(::MPolyRingElem{T}) where T <: RingElement
```

```@docs
is_constant(::MPolyRingElem{T}) where T <: RingElement
```

```@docs
is_term(::MPolyRingElem{T}) where T <: RingElement
```

```@docs
is_monomial(::MPolyRingElem{T}) where T <: RingElement
```

```@docs
is_univariate(::MPolyRingElem{T}) where T <: RingElement
```

```@docs
coeff(::MPolyRingElem{T}, ::MPolyRingElem{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = x^2 + 2x + 1
x^2 + 2*x + 1

julia> V = vars(f)
1-element Vector{AbstractAlgebra.Generic.MPoly{BigInt}}:
 x

julia> var_index(y) == 2
true

julia> degree(f, x) == 2
true

julia> degree(f, 2) == 0
true

julia> d = degrees(f)
2-element Vector{Int64}:
 2
 0

julia> is_constant(R(1))
true

julia> is_term(2x)
true

julia> is_monomial(y)
true

julia> is_unit(R(1))
true

julia> S, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = x^3*y + 3x*y^2 + 1
x^3*y + 3*x*y^2 + 1

julia> c1 = coeff(f, 1)
1

julia> c2 = coeff(f, x^3*y)
1

julia> m = monomial(f, 2)
x*y^2

julia> e1 = exponent(f, 1, 1)
3

julia> v1 = exponent_vector(f, 1)
2-element Vector{Int64}:
 3
 1

julia> t1 = term(f, 1)
x^3*y

julia> setcoeff!(f, [3, 1], 12)
12*x^3*y + 3*x*y^2 + 1

julia> S, (x, y) = polynomial_ring(QQ, ["x", "y"]; ordering=:deglex)
(Multivariate polynomial ring in 2 variables over rationals, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}[x, y])

julia> V = symbols(S)
2-element Vector{Symbol}:
 :x
 :y

julia> X = gens(S)
2-element Vector{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}:
 x
 y

julia> ord = ordering(S)
:deglex

julia> S, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = x^3*y + 3x*y^2 + 1
x^3*y + 3*x*y^2 + 1

julia> n = length(f)
3

julia> is_gen(y)
true

julia> number_of_variables(S) == 2
true

julia> d = total_degree(f)
4

julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = 2x^2*y + 2x + y + 1
2*x^2*y + 2*x + y + 1

julia> g = x^2*y^2 + 1
x^2*y^2 + 1

julia> flag, q = divides(f*g, f)
(true, x^2*y^2 + 1)

julia> d = divexact(f*g, f)
x^2*y^2 + 1

julia> v, q = remove(f*g^3, g)
(3, 2*x^2*y + 2*x + y + 1)

julia> n = valuation(f*g^3, g)
3

julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over rationals, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}[x, y])

julia> f = 3x^2*y^2 + 2x + 1
3*x^2*y^2 + 2*x + 1

julia> f1 = divexact(f, 5)
3//5*x^2*y^2 + 2//5*x + 1//5

julia> f2 = divexact(f, QQ(2, 3))
9//2*x^2*y^2 + 3*x + 3//2
```

### Square root

Over rings for which an exact square root is available, it is possible to take
the square root of a polynomial or test whether it is a square.

```julia
sqrt(f::MPolyRingElem, check::bool=true)
is_square(::MPolyRingElem)
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = -4*x^5*y^4 + 5*x^5*y^3 + 4*x^4 - x^3*y^4
-4*x^5*y^4 + 5*x^5*y^3 + 4*x^4 - x^3*y^4

julia> sqrt(f^2)
4*x^5*y^4 - 5*x^5*y^3 - 4*x^4 + x^3*y^4

julia> is_square(f)
false
```

### Iterators

The following iterators are provided for multivariate polynomials.

```julia
coefficients(p::MPoly)
monomials(p::MPoly)
terms(p::MPoly)
exponent_vectors(a::MPoly)
```

**Examples**

```jldoctest
julia> S, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = x^3*y + 3x*y^2 + 1
x^3*y + 3*x*y^2 + 1

julia> C = collect(coefficients(f))
3-element Vector{BigInt}:
 1
 3
 1

julia> M = collect(monomials(f))
3-element Vector{AbstractAlgebra.Generic.MPoly{BigInt}}:
 x^3*y
 x*y^2
 1

julia> T = collect(terms(f))
3-element Vector{AbstractAlgebra.Generic.MPoly{BigInt}}:
 x^3*y
 3*x*y^2
 1

julia> V = collect(exponent_vectors(f))
3-element Vector{Vector{Int64}}:
 [3, 1]
 [1, 2]
 [0, 0]
```

### Changing base (coefficient) rings

In order to substitute the variables of a polynomial $f$ over a ring $T$ by
elements in a $T$-algebra $S$, you first have to change the base ring of $f$
using the following function, where $g$ is a function representing the
structure homomorphism of the $T$-algebra $S$.

```@docs
change_base_ring(::Ring, p::MPolyRingElem{T}) where {T <: RingElement}
change_coefficient_ring(::Ring, p::MPolyRingElem{T}) where {T <: RingElement}
map_coefficients(::Any, p::MPolyRingElem)
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> fz = x^2*y^2 + x + 1
x^2*y^2 + x + 1

julia> fq = change_base_ring(QQ, fz)
x^2*y^2 + x + 1

julia> fq = change_coefficient_ring(QQ, fz)
x^2*y^2 + x + 1

```

In case a specific parent ring is constructed, it can also be passed to the function.

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> S,  = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over rationals, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}[x, y])

julia> fz = x^5 + y^3 + 1
x^5 + y^3 + 1

julia> fq = change_base_ring(QQ, fz, parent=S)
x^5 + y^3 + 1
```

### Multivariate coefficients

In order to return the "coefficient" (as a multivariate polynomial in the same
ring), of a given monomial (in which some of the variables may not appear and
others may be required to appear to exponent zero), we can use the following
function.

```@docs
coeff(a::MPolyRingElem{T}, vars::Vector{Int}, exps::Vector{Int}) where T <: RingElement
coeff(a::T, vars::Vector{T}, exps::Vector{Int}) where T <: MPolyRingElem
```

**Examples**

```jldoctest
julia> R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y, z])

julia> f = x^4*y^2*z^2 - 2x^4*y*z^2 + 4x^4*z^2 + 2x^2*y^2 + x + 1
x^4*y^2*z^2 - 2*x^4*y*z^2 + 4*x^4*z^2 + 2*x^2*y^2 + x + 1

julia> coeff(f, [1, 3], [4, 2]) == coeff(f, [x, z], [4, 2])
true

```

### Inflation/deflation

```@docs
deflation(f::MPolyRingElem{T}) where T <: RingElement
```

```@docs
deflate(f::MPolyRingElem{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: RingElement
deflate(f::MPolyRingElem{T}, defl::Vector{Int}) where T <: RingElement
deflate(f::MPolyRingElem{T}) where T <: RingElement
deflate(f::MPolyRingElem, vars::Vector{Int}, shift::Vector{Int}, defl::Vector{Int})
deflate(f::T, vars::Vector{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: MPolyRingElem
```

```@docs
inflate(f::MPolyRingElem{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: RingElement
inflate(f::MPolyRingElem{T}, defl::Vector{Int}) where T <: RingElement
inflate(f::MPolyRingElem, vars::Vector{Int}, shift::Vector{Int}, defl::Vector{Int})
inflate(f::T, vars::Vector{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: MPolyRingElem
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = x^7*y^8 + 3*x^4*y^8 - x^4*y^2 + 5x*y^5 - x*y^2
x^7*y^8 + 3*x^4*y^8 - x^4*y^2 + 5*x*y^5 - x*y^2

julia> def, shift = deflation(f)
([1, 2], [3, 3])

julia> f1 = deflate(f, def, shift)
x^2*y^2 + 3*x*y^2 - x + 5*y - 1

julia> f2 = inflate(f1, def, shift)
x^7*y^8 + 3*x^4*y^8 - x^4*y^2 + 5*x*y^5 - x*y^2

julia> f2 == f
true

julia> g = (x+y+1)^2
x^2 + 2*x*y + 2*x + y^2 + 2*y + 1

julia> g0 = coeff(g, [y], [0])
x^2 + 2*x + 1

julia> g1 = deflate(g - g0, [y], [1], [1])
2*x + y + 2

julia> g == g0 + y * g1
true

```

### Conversions

```@docs
to_univariate(R::PolyRing{T}, p::MPolyRingElem{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> S, z = polynomial_ring(ZZ, "z")
(Univariate polynomial ring in z over integers, z)

julia> f = 2x^5 + 3x^4 - 2x^2 - 1
2*x^5 + 3*x^4 - 2*x^2 - 1

julia> g = to_univariate(S, f)
2*z^5 + 3*z^4 - 2*z^2 - 1

```

### Evaluation

The following function allows evaluation of a polynomial at all its variables.
The result is always in the ring that a product of a coefficient and one of the
values belongs to, i.e. if all the values are in the coefficient ring, the
result of the evaluation will be too.

```@docs
evaluate(::MPolyRingElem{T}, ::Vector{U}) where {T <: RingElement, U <: RingElement}
```

The following functions allow evaluation of a polynomial at some of its
variables. Note that the result will be a product of values and an element
of the polynomial ring, i.e. even if all the values are in the coefficient
ring and all variables are given values, the result will be a constant
polynomial, not a coefficient.

```@docs
evaluate(::MPolyRingElem{T}, ::Vector{Int}, ::Vector{U}) where {T <: RingElement, U <: RingElement}
```

```@docs
evaluate(::S, ::Vector{S}, ::Vector{U}) where {S <: MPolyRingElem{T}, U <: RingElement} where T <: RingElement
```

The following function allows evaluation of a polynomial at values in a
not necessarily commutative ring, e.g. elements of a matrix algebra.

```@docs
evaluate(::MPolyRingElem{T}, ::Vector{U}) where {T <: RingElement, U <: NCRingElem}
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = 2x^2*y^2 + 3x + y + 1
2*x^2*y^2 + 3*x + y + 1

julia> evaluate(f, BigInt[1, 2])
14

julia> evaluate(f, [QQ(1), QQ(2)])
14//1

julia> evaluate(f, [1, 2])
14

julia> f(1, 2) == 14
true

julia> evaluate(f, [x + y, 2y - x])
2*x^4 - 4*x^3*y - 6*x^2*y^2 + 8*x*y^3 + 2*x + 8*y^4 + 5*y + 1

julia> f(x + y, 2y - x)
2*x^4 - 4*x^3*y - 6*x^2*y^2 + 8*x*y^3 + 2*x + 8*y^4 + 5*y + 1

julia> R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y, z])

julia> f = x^2*y^2 + 2x*z + 3y*z + z + 1
x^2*y^2 + 2*x*z + 3*y*z + z + 1

julia> evaluate(f, [1, 3], [3, 4])
9*y^2 + 12*y + 29

julia> evaluate(f, [x, z], [3, 4])
9*y^2 + 12*y + 29

julia> evaluate(f, [1, 2], [x + z, x - z])
x^4 - 2*x^2*z^2 + 5*x*z + z^4 - z^2 + z + 1

julia> S = MatrixAlgebra(ZZ, 2)
Matrix algebra of degree 2
  over integers

julia> M1 = S([1 2; 3 4])
[1   2]
[3   4]

julia> M2 = S([2 3; 1 -1])
[2    3]
[1   -1]

julia> M3 = S([-1 1; 1 1])
[-1   1]
[ 1   1]

julia> evaluate(f, [M1, M2, M3])
[ 64    83]
[124   149]
```

### Leading and constant coefficients, leading monomials and leading terms

The leading and trailing coefficient, constant coefficient, leading monomial
and leading term of a polynomial p are returned by the following functions:

```@docs
leading_coefficient(::MPolyRingElem{T}) where T <: RingElement
trailing_coefficient(p::MPolyRingElem{T}) where T <: RingElement
leading_monomial(::MPolyRingElem{T}) where T <: RingElement
leading_term(::MPolyRingElem{T}) where T <: RingElement
constant_coefficient(::MPolyRingElem{T}) where T <: RingElement
tail(::MPolyRingElem{T}) where T <: RingElement
```

**Examples**


```julia
using AbstractAlgebra
R,(x,y) = polynomial_ring(ZZ, ["x", "y"], ordering=:deglex)
p = 2*x*y + 3*y^3 + 1
leading_term(p)
leading_monomial(p)
leading_coefficient(p)
leading_term(p) == leading_coefficient(p) * leading_monomial(p)
constant_coefficient(p)
tail(p)
```

### Least common multiple, greatest common divisor

The greatest common divisor of two polynomials a and b is returned by

```@docs
gcd(a::Generic.MPoly{T}, b::Generic.MPoly{T}) where {T <: RingElement}
```

Note that this functionality is currently only provided for AbstractAlgebra
generic polynomials. It is not automatically provided for all multivariate
rings that implement the multivariate interface.

However, if such a gcd is provided, the least common multiple of two
polynomials a and b is returned by

```@docs
lcm(a::MPolyRingElem{T}, b::MPolyRingElem{T}) where {T <: RingElement}
```

**Examples**

```jldoctest
julia> using AbstractAlgebra

julia> R,(x,y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> a = x*y + 2*y
x*y + 2*y

julia> b = x^3*y + y
x^3*y + y

julia> gcd(a,b)
y

julia> lcm(a,b)
x^4*y + 2*x^3*y + x*y + 2*y

julia> lcm(a,b) == a * b // gcd(a,b)
true

```

### Derivations

```@docs
derivative(::MPolyRingElem{T}, ::MPolyRingElem{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> R, (x, y) = AbstractAlgebra.polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = x*y + x + y + 1
x*y + x + y + 1

julia> derivative(f, x)
y + 1

julia> derivative(f, y)
x + 1

julia> derivative(f, 1)
y + 1

julia> derivative(f, 2)
x + 1
```

### Homogeneous polynomials

It is possible to test whether a polynomial is homogeneous with respect to the standard grading using the function

```@docs
is_homogeneous(x::MPolyRingElem{T}) where T <: RingElement
```

## Random generation

Random multivariate polynomials in a given ring can be constructed by passing
a range of degrees for the variables and a range on the number of terms.
Additional parameters are used to generate the coefficients of the polynomial.

Note that zero coefficients may currently be generated, leading to less than
the requested number of terms.

```julia
rand(R::MPolyRing, exp_range::AbstractUnitRange{Int}, term_range::AbstractUnitRange{Int}, v...)
```

**Examples**

```jldoctest; setup = :(import Random; Random.seed!(42))
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = rand(R, -1:2, 3:5, -10:10)
4*x^4*y^4

julia> S, (s, t) = polynomial_ring(GF(7), ["x", "y"])
(Multivariate polynomial ring in 2 variables over finite field F_7, AbstractAlgebra.Generic.MPoly{AbstractAlgebra.GFElem{Int64}}[x, y])

julia> g = rand(S, -1:2, 3:5)
4*x^3*y^4
```
