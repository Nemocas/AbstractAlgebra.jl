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

The polynomial element types belong to the abstract type `MPolyElem{T}`
and the polynomial ring types belong to the abstract type `MPolyRing{T}`.

Note that both the generic polynomial ring type `Generic.MPolyRing{T}` and the abstract
type it belongs to, `MPolyRing{T}` are both called `MPolyRing`. The
former is a (parameterised) concrete type for a polynomial ring over a given base ring
whose elements have type `T`. The latter is an abstract type representing all
multivariate polynomial ring types in AbstractAlgebra.jl, whether generic or very
specialised (e.g. supplied by a C library).

## Polynomial ring constructors

In order to construct multivariate polynomials in AbstractAlgebra.jl, one must first
construct the polynomial ring itself. This is accomplished with one of the following
constructors.

```julia
PolynomialRing(R::Ring, S::Array{String, 1}; cached::Bool = true, ordering::Symbol=:lex)
PolynomialRing(R::Ring, n::Int, s::String="x"; cached::Bool = false, ordering::Symbol = :lex)
```

The first constructor, given a base ring `R` and an array `S` of strings
specifying how the generators (variables) should be printed, will return a
tuple `S, (x, ...)` representing the new polynomial ring $S = R[x, \ldots]$ and
a tuple of the generators $(x, ...)$ of the ring.

The second constructor given a string `s` and a number of variables `n` will
do the same as the first constructor except that the variables will be
automatically numbered. For example if `s` is the string `x` and `n = 3` then
the variables will print as `x1`, `x2`, `x3`.

By default the parent object `S` will depend only on `R` and  `(x, ...)` and
will be cached. Setting the optional argument `cached` to `false` will prevent
the parent object `S` from being cached.

The optional named argument `ordering` can be used to specify an ordering. The
currently supported options are `:lex`, `:deglex` and `:degrevlex`.

Like for univariate polynomials, a shorthand version of this function is
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
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"]; ordering=:deglex)
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> T, (z, t) = QQ["z", "t"]
(Multivariate Polynomial Ring in z, t over Rationals, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}[z, t])

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

```

## Polynomial constructors

Multivariate polynomials can be constructed from the generators in the
usual way using arithmetic operations.

Also, all of the standard ring element constructors may be used to construct
multivariate polynomials.

```julia
(R::MPolyRing)() # constructs zero
(R::MPolyRing)(c::Integer)
(R::MPolyRing)(c::elem_type(R))
(R::MPolyRing{T})(a::T) where T <: RingElement
```

For more efficient construction of multivariate polynomial, one can use the
`MPoly` build context, where terms (coefficient followed by an exponent vector)
are pushed onto a context one at a time and then the polynomial constructed
from those terms in one go using the `finish` function.

```julia
MPolyBuildCtx(R::MPolyRing)
push_term!(M::MPolyBuildCtx, c::RingElem, v::Vector{Int})
finish(M::MPolyBuildCtx)
```

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
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> C = MPolyBuildCtx(R)
Builder for a polynomial in Multivariate Polynomial Ring in x, y over Integers

julia> push_term!(C, ZZ(3), [1, 2]);


julia> push_term!(C, ZZ(2), [1, 1]);


julia> push_term!(C, ZZ(4), [0, 0]); 


julia> f = finish(C)
3*x*y^2 + 2*x*y + 4

julia> S, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rationals, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}[x, y])

julia> f = S(Rational{BigInt}[2, 3, 1], [[3, 2], [1, 0], [0, 1]])
2*x^3*y^2 + 3*x + y
```

## Functions for types and parents of multivariate polynomial rings

```julia
base_ring(R::MPolyRing)
base_ring(a::MPolyElem)
```

Return the coefficient ring of the given polynomial ring or polynomial,
respectively.


```julia
parent(a::MPolyElem)
```

Return the polynomial ring of the given polynomial.

```julia
characteristic(R::MPolyRing)
```

Return the characteristic of the given polynomial ring.

## Polynomial functions

### Basic manipulation

All the standard ring functions are available, including the following.

```julia
zero(R::MPolyRing)
one(R::MPolyRing)
iszero(a::MPolyElem)
isone(a::MPolyElem)
```

```
divexact(a::T, b::T) where T <: MPolyElem
```

All basic functions from the Multivariate Polynomial interface are provided.

```julia
symbols(S::MPolyRing)
nvars(f::MPolyRing)
gens(S::MPolyRing)
gen(S::MPolyRing, i::Int)
```

```julia
ordering(S::MPolyRing{T})
```

Note that the currently supported orderings are `:lex`, `:deglex` and
`:degrevlex`.

```julia
length(f::MPolyElem)
degrees(f::MPolyElem)
total_degree(f::MPolyElem)
```

```julia
isgen(x::MPolyElem)
```

```julia
divexact(f::T, g::T) where T <: MPolyElem
```

For multivariate polynomial types that allow constant time access to
coefficients, the following are also available, allowing access to the given
coefficient, monomial or term. Terms are numbered from the most significant
first.

```julia
coeff(f::MPolyElem, n::Int)
coeff(a::MPolyElem, exps::Vector{Int})
```

Access a coefficient by term number or exponent vector.

```julia
monomial(f::MPolyElem, n::Int)
monomial!(m::T, f::T, n::Int) where T <: MPolyElem
```

The second version writes the result into a preexisting polynomial
object to save an allocation.

```julia
term(f::MPolyElem, n::Int)
```

```julia
exponent(f::MyMPolyElem}, i::Int, j::Int)
```

Return the exponent of the $j$-th variable in the $i$-th term of the polynomial
$f$.

```julia
exponent_vector(a::MPolyElem, i::Int)
```

```julia
setcoeff!(a::MPolyElem{T}, exps::Vector{Int}, c::T) where T <: RingElement
```

Although multivariate polynomial rings are not usually Euclidean, the following
functions from the Euclidean interface are often provided.

```julia
divides(f::T, g::T) where T <: MPolyElem
remove(f::T, g::T) where T <: MPolyElem
valuation(f::T, g::T) where T <: MPolyElem
```

```julia
divrem(f::T, g::T) where T <: MPolyElem
div(f::T, g::T) where T <: MPolyElem
```

Compute a tuple $(q, r)$ such that $f = qg + r$, where the coefficients of terms of
$r$ whose monomials are divisible by the leading monomial of $g$ are reduced modulo the
leading coefficient of $g$ (according to the Euclidean function on the coefficients).
The `divrem` version returns both quotient and remainder whilst the `div` version only
returns the quotient.

Note that the result of these functions depend on the ordering of the polynomial ring.

```julia
gcd(f::T, g::T) where T <: MPolyElem
```

The following functionality is also provided for all multivariate polynomials.

```@docs
isunivariate(::MPolyRing{T}) where T <: RingElement
```

```@docs
vars(p::MPolyElem{T}) where T <: RingElement
```

```@docs
var_index(::MPolyElem{T}) where T <: RingElement
```

```@docs
degree(::MPolyElem{T}, ::Int) where T <: RingElement
```

```@docs
degree(::MPolyElem{T}, ::MPolyElem{T}) where T <: RingElement
```

```@docs
degrees(::MPolyElem{T}) where T <: RingElement
```

```@docs
isconstant(::MPolyElem{T}) where T <: RingElement
```

```@docs
isterm(::MPolyElem{T}) where T <: RingElement
```

```@docs
ismonomial(::MPolyElem{T}) where T <: RingElement
```

```@docs
isunivariate(::MPolyElem{T}) where T <: RingElement
```

```@docs
coeff(::MPolyElem{T}, ::MPolyElem{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = x^2 + 2x + 1
x^2 + 2*x + 1

julia> V = vars(f)
1-element Array{AbstractAlgebra.Generic.MPoly{BigInt},1}:
 x

julia> var_index(y) == 2
true

julia> degree(f, x) == 2
true

julia> degree(f, 2) == 0
true

julia> d = degrees(f)
2-element Array{Int64,1}:
 2
 0

julia> isconstant(R(1))
true

julia> isterm(2x)
true

julia> ismonomial(y)
true

julia> isunit(R(1))
true

julia> c = coeff(f, x^2)
1

```

### Square root

Over rings for which an exact square root is available, it is possible to take
the square root of a polynomial or test whether it is a square.

```julia
sqrt(f::MPolyElem, check::bool=true)
issquare(::MPolyElem)
```

### Iterators

The following iterators are provided for multivariate polynomials.

```julia
coefficients(p::MPoly)
monomials(p::MPoly)
terms(p::MPoly)
exponent_vectors(a::MPoly)
```

### Changing base (coefficient) rings

In order to substitute the variables of a polynomial $f$ over a ring $T$ by
elements in a $T$-algebra $S$, you first have to change the base ring of $f$
using the following function, where $g$ is a function representing the
structure homomorphism of the $T$-algebra $S$.

```@docs
change_base_ring(::Ring, p::MPolyElem{T}) where {T <: RingElement}
change_coefficient_ring(::Ring, p::MPolyElem{T}) where {T <: RingElement}
map_coefficients(::Any, p::MPolyElem)
```

**Examples**

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

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
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> S,  = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rationals, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}[x, y])

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
coeff(a::MPolyElem{T}, vars::Vector{Int}, exps::Vector{Int}) where T <: RingElement
coeff(a::T, vars::Vector{T}, exps::Vector{Int}) where T <: MPolyElem
```

**Examples**

```jldoctest
julia> R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y, z])

julia> f = x^4*y^2*z^2 - 2x^4*y*z^2 + 4x^4*z^2 + 2x^2*y^2 + x + 1
x^4*y^2*z^2 - 2*x^4*y*z^2 + 4*x^4*z^2 + 2*x^2*y^2 + x + 1

julia> coeff(f, [1, 3], [4, 2]) == coeff(f, [x, z], [4, 2])
true

```

### Inflation/deflation

```@docs
deflation(f::MPolyElem{T}) where T <: RingElement
```

```@docs
deflate(f::MPolyElem{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: RingElement
deflate(f::MPolyElem{T}, defl::Vector{Int}) where T <: RingElement
deflate(f::MPolyElem{T}) where T <: RingElement
```

```@docs
inflate(f::MPolyElem{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: RingElement
inflate(f::MPolyElem{T}, defl::Vector{Int}) where T <: RingElement
```

**Examples**

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

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

```

### Conversions

```@docs
to_univariate(R::PolyRing{T}, p::MPolyElem{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> S, z = PolynomialRing(ZZ, "z")
(Univariate Polynomial Ring in z over Integers, z)

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
evaluate(::MPolyElem{T}, ::Vector{U}) where {T <: RingElement, U <: RingElement}
```

The following functions allow evaluation of a polynomial at some of its
variables. Note that the result will be a product of values and an element
of the polynomial ring, i.e. even if all the values are in the coefficient
ring and all variables are given values, the result will be a constant
polynomial, not a coefficient.

```@docs
evaluate(::MPolyElem{T}, ::Vector{Int}, ::Vector{U}) where {T <: RingElement, U <: RingElement}
```

```@docs
evaluate(::S, ::Vector{S}, ::Vector{U}) where {S <: MPolyElem{T}, U <: RingElement} where T <: RingElement
```

The following function allows evaluation of a polynomial at values in a
not necessarily commutative ring, e.g. elements of a matrix algebra.

```@docs
evaluate(::MPolyElem{T}, ::Vector{U}) where {T <: RingElement, U <: NCRingElem}
```

**Examples**

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

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

julia> R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y, z])

julia> f = x^2*y^2 + 2x*z + 3y*z + z + 1
x^2*y^2 + 2*x*z + 3*y*z + z + 1

julia> evaluate(f, [1, 3], [3, 4])
9*y^2 + 12*y + 29

julia> evaluate(f, [x, z], [3, 4])
9*y^2 + 12*y + 29

julia> evaluate(f, [1, 2], [x + z, x - z])
x^4 - 2*x^2*z^2 + 5*x*z + z^4 - z^2 + z + 1

julia> S = MatrixAlgebra(ZZ, 2)
Matrix Algebra of degree 2 over Integers

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
leading_coefficient(::MPolyElem{T}) where T <: RingElement
trailing_coefficient(p::MPolyElem{T}) where T <: RingElement
leading_monomial(::MPolyElem{T}) where T <: RingElement
leading_term(::MPolyElem{T}) where T <: RingElement
constant_coefficient(::MPolyElem{T}) where T <: RingElement
tail(::MPolyElem{T}) where T <: RingElement
```

**Examples**


```julia
using AbstractAlgebra
R,(x,y) = PolynomialRing(ZZ, ["x", "y"], ordering=:deglex)
p = 2*x*y + 3*y^3 + 1
leading_term(p)
leading_monomial(p)
leading_coefficient(p)
leading_term(p) == leading_coefficient(p) * leading_monomial(p)
constant_coefficient(p)
tail(p)
```

### Least common multiple, greatest common divisor

The greated common divisor of two polynomials a and b is returned by

```@docs
gcd(a::Generic.MPoly{T}, b::Generic.MPoly{T}) where {T <: RingElement}
```

Note that this functionality is currently only provided for AbstractAlgebra
generic polynomials. It is not automatically provided for all multivariate
rings that implement the multivariate interface.

However, if such a gcd is provided, the least common multiple of two
polynomials a and b is returned by

```@docs
lcm(a::MPolyElem{T}, b::MPolyElem{T}) where {T <: RingElement}
```

**Examples**

```jldoctest
julia> using AbstractAlgebra

julia> R,(x,y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

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
derivative(::MPolyElem{T}, ::MPolyElem{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> R, (x, y) = AbstractAlgebra.PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = x*y + x + y + 1
x*y + x + y + 1

julia> derivative(f, x)
y + 1

julia> derivative(f, y)
x + 1

```

### Homogeneous polynomials

It is possible to test whether a polynomial is homogeneous with respect to the standard grading using the function

```@docs
ishomogeneous(x::MPolyElem{T}) where T <: RingElement
```

## Random generation

Random multivariate polynomials in a given ring can be constructed by passing
a range of degrees for the variables and a range on the number of terms.
Additional parameters are used to generate the coefficients of the polynomial.

Note that zero coefficients may currently be generated, leading to less than
the requested number of terms.

```julia
rand(R::MPolyRing, exp_range::UnitRange{Int}, term_range::UnitRange{Int}, v...)
```

** Examples **

```@repl
R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
f = rand(R, -1:2, 3:5, -10:10)

S, (s, t) = PolynomialRing(GF(7), ["x", "y"])
g = rand(S, -1:2, 3:5)
```
