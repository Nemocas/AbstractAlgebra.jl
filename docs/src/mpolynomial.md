```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Generic sparse distributed multivariate polynomials

AbstractAlgebra.jl provides a module, implemented in `src/generic/MPoly.jl` for generic
sparse distributed multivariate polynomials over any commutative ring belonging to the
AbstractAlgebra abstract type hierarchy.

This modules implements the Multivariate Polynomial interface, including the sparse
distributed, random access part of the interface.

All of the generic functionality is part of a submodule of AbstractAlgebra called
`Generic`. This is exported by default so that it is not necessary to qualify the
function names with the submodule name.

Multivariates are implemented in this module using a Julia array of coefficients and
a 2-dimensional Julia array of `UInt`s for the exponent vectors. Note that exponent $n$
is represented by the $n$-th column of the exponent array, not the $n$-th row. This is
because Julia uses a column major representation.

## Types and parent objects

Multivariate polynomials implemented in AbstractAlgebra.jl have type `Generic.MPoly{T}`
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

The polynomial element types belong to the abstract type `AbstractAlgebra.MPolyElem{T}`
and the polynomial ring types belong to the abstract type `AbstractAlgebra.MPolyRing{T}`.

Note that both the generic polynomial ring type `Generic.MPolyRing{T}` and the abstract
type it belongs to, `AbstractAlgebra.MPolyRing{T}` are both called `MPolyRing`. The
former is a (parameterised) concrete type for a polynomial ring over a given base ring
whose elements have type `T`. The latter is an abstract type representing all
multivariate polynomial ring types in AbstractAlgebra.jl, whether generic or very
specialised (e.g. supplied by a C library).

## Polynomial ring constructors

In order to construct multivariate polynomials in AbstractAlgebra.jl, one must first
construct the polynomial ring itself. This is accomplished with the following
constructor.

```julia
PolynomialRing(R::AbstractAlgebra.Ring, S::Array{String, 1}; cached::Bool = true, ordering::Symbol=:lex)
```

Given a base ring `R` and and array `S` of strings specifying how the generators
(variables) should be printed, return a tuple `S, (x, ...)` representing the new
polynomial ring $S = R[x, \ldots]$ and a tuple of the generators $(x, ...)$ of the ring.
By default the parent object `S` will depend only on `R` and  `(x, ...)` and will be
cached. Setting the optional argument `cached` to `false` will prevent the parent object
`S` from being cached.

The optional named argument `ordering` can be used to specify an ordering. The currently
supported options are `:lex`, `:deglex` and `:degrevlex`.

Here are some examples of creating multivariate polynomial rings and making use of the
resulting parent objects to coerce various elements into the polynomial ring.

**Examples**

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"]; ordering=:deglex)
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = R()
0

julia> g = R(123)
123

julia> h = R(BigInt(1234))
1234

julia> k = R(x + 1)
x+1

julia> m = R(x + y + 1)
x+y+1

julia> derivative(k, 1)
1

julia> derivative(k, 2)
0

```

All of the examples here are generic polynomial rings, but specialised implementations
of polynomial rings provided by external modules will also usually provide a
`PolynomialRing` constructor to allow creation of their polynomial rings.

## Polynomial functionality provided by AbstractAlgebra.jl

### Basic manipulation

```@docs
vars(p::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```

```@docs
var_index(::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```

```@docs
degree(::AbstractAlgebra.MPolyElem{T}, ::Int) where T <: RingElement
```

```@docs
degree(::AbstractAlgebra.MPolyElem{T}, ::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```

```@docs
degrees(::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```

```@docs
isconstant(::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```

```@docs
isterm(::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```

```@docs
ismonomial(::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```

```@docs
coeff(::AbstractAlgebra.MPolyElem{T}, ::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = x^2 + 2x + 1
x^2+2*x+1

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

### Changing base (coefficient) rings

In order to substitute the variables of a polynomial $f$ over a ring $T$ by
elements in a $T$-algebra $S$, you first have to change the base ring of $f$
using the following function, where $g$ is a function representing the
structure homomorphism of the $T$-algebra $S$.

```@docs
change_base_ring(::Ring, p::AbstractAlgebra.MPolyElem{T}) where {T <: RingElement}
map_coeffs(::Any, p::AbstractAlgebra.MPolyElem)
```

**Examples**

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> fz = x^2*y^2 + x + 1
x^2*y^2+x+1

julia> fq = change_base_ring(QQ, fz)
x^2*y^2+x+1//1

```

Incase a specific parent ring is constructed, it can also be passed to the function.

```

**Examples**

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> S,  = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rationals, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}[x, y])

julia> fz = x^5 + y^3
x^5+y^3+1

julia> fq = change_base_ring(QQ, fz, S)
x^5+y^3+1//1


```

### Multivariate coefficients

In order to return the "coefficient" (as a multivariate polynomial in the same
ring), of a given monomial (in which some of the variables may not appear and
others may be required to appear to exponent zero), we can use the following
function.

```@docs
coeff(a::AbstractAlgebra.MPolyElem{T}, vars::Vector{Int}, exps::Vector{Int}) where T <: RingElement
coeff(a::T, vars::Vector{T}, exps::Vector{Int}) where T <: AbstractAlgebra.MPolyElem
```

**Examples**

```jldoctest
julia> R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y, z])

julia> f = x^4*y^2*z^2 - 2x^4*y*z^2 + 4x^4*z^2 + 2x^2*y^2 + x + 1
x^4*y^2*z^2-2*x^4*y*z^2+4*x^4*z^2+2*x^2*y^2+x+1

julia> coeff(f, [1, 3], [4, 2]) == coeff(f, [x, z], [4, 2])
true

```

### Inflation/deflation

```@docs
deflation(f::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```

```@docs
deflate(f::AbstractAlgebra.MPolyElem{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: RingElement
```

```@docs
inflate(f::AbstractAlgebra.MPolyElem{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: RingElement
```

**Examples**

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = x^7*y^8 + 3*x^4*y^8 - x^4*y^2 + 5x*y^5 - x*y^2
x^7*y^8+3*x^4*y^8-x^4*y^2+5*x*y^5-x*y^2

julia> def, shift = deflation(f)
([1, 2], [3, 3])

julia> f1 = deflate(f, def, shift)
x^2*y^2+3*x*y^2-x+5*y-1

julia> f2 = inflate(f1, def, shift)
x^7*y^8+3*x^4*y^8-x^4*y^2+5*x*y^5-x*y^2

julia> f2 == f
true

```

### Conversions

```@docs
to_univariate(R::AbstractAlgebra.PolyRing{T}, p::AbstractAlgebra.MPolyElem{T}) where T <: AbstractAlgebra.RingElement
```

**Examples**

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> S, z = PolynomialRing(ZZ, "z")
(Univariate Polynomial Ring in z over Integers, z)

julia> f = 2x^5 + 3x^4 - 2x^2 - 1
2*x^5+3*x^4-2*x^2-1

julia> g = to_univariate(S, f)
2*z^5+3*z^4-2*z^2-1

```

### Evaluation

The following function allows evaluation of a polynomial at all its variables.
The result is always in the ring that a product of a coefficient and one of the
values belongs to, i.e. if all the values are in the coefficient ring, the
result of the evaluation will be too.

```@docs
evaluate(::AbstractAlgebra.MPolyElem{T}, ::Vector{U}) where {T <: RingElement, U <: RingElement}
```

The following functions allow evaluation of a polynomial at some of its
variables. Note that the result will be a product of values and an element
of the polynomial ring, i.e. even if all the values are in the coefficient
ring and all variables are given values, the result will be a constant
polynomial, not a coefficient.

```@docs
evaluate(::AbstractAlgebra.MPolyElem{T}, ::Vector{Int}, ::Vector{U}) where {T <: RingElement, U <: RingElement}
```

```@docs
evaluate(::S, ::Vector{S}, ::Vector{U}) where {S <: AbstractAlgebra.MPolyElem{T}, U <: RingElement} where T <: RingElement
```

The following function allows evaluation of a polynomial at values in a
not necessarily commutative ring, e.g. elements of a matrix algebra.

```@docs
evaluate(::AbstractAlgebra.MPolyElem{T}, ::Vector{U}) where {T <: RingElement, U <: NCRingElem}
```

**Examples**

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = 2x^2*y^2 + 3x + y + 1
2*x^2*y^2+3*x+y+1

julia> evaluate(f, BigInt[1, 2])
14

julia> evaluate(f, [QQ(1), QQ(2)])
14//1

julia> evaluate(f, [1, 2])
14

julia> f(1, 2) == 14
true

julia> evaluate(f, [x + y, 2y - x])
2*x^4-4*x^3*y-6*x^2*y^2+8*x*y^3+2*x+8*y^4+5*y+1

julia> f(x + y, 2y - x)
2*x^4-4*x^3*y-6*x^2*y^2+8*x*y^3+2*x+8*y^4+5*y+1

julia> R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y, z])

julia> f = x^2*y^2 + 2x*z + 3y*z + z + 1
x^2*y^2+2*x*z+3*y*z+z+1

julia> evaluate(f, [1, 3], [3, 4])
9*y^2+12*y+29

julia> evaluate(f, [x, z], [3, 4])
9*y^2+12*y+29

julia> evaluate(f, [1, 2], [x + z, x - z])
x^4-2*x^2*z^2+5*x*z+z^4-z^2+z+1

julia> evaluate(f, [2, 4, 6], QQ)
167//1

julia> evaluate(f, [x, z], [2, 4], QQ)
4//1*y^2+12//1*y+21//1

julia> S = MatrixAlgebra(ZZ, 2)
Matrix Algebra of degree 2 over Integers

julia> M1 = S([1 2; 3 4])
[1  2]
[3  4]

julia> M2 = S([2 3; 1 -1])
[2   3]
[1  -1]

julia> M3 = S([-1 1; 1 1])
[-1  1]
[ 1  1]

julia> evaluate(f, [M1, M2, M3])
[ 64   83]
[124  149]
```

### Leading coefficients, leading monomials and leading terms

The leading coefficient, leading monomial and leading term of a polynomial p are returned by the following functions:

```@docs
lc(p::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```
```@docs
lm(p::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```
```@docs
lt(p::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```

**Examples**


```julia
using AbstractAlgebra
R,(x,y) = PolynomialRing(ZZ, ["x", "y"], ordering=:deglex)
p = 2*x*y + 3*y^3
lt(p)
lm(p)
lc(p)
lt(p) == lc(p) * lm(p)
```

### Least common multiple, greatest common divisor

The greated common divisor of two polynomials a and b is returned by

```@docs
gcd(a::AbstractAlgebra.Generic.MPoly{T}, b::AbstractAlgebra.Generic.MPoly{T}) where {T <: RingElement}
```

Note that this functionality is currently only provided for AbstractAlgebra
generic polynomials. It is not automatically provided for all multivariate
rings that implement the multivariate interface.

However, if such a gcd is provided, the least common multiple of two
polynomials a and b is returned by

```@docs
lcm(a::AbstractAlgebra.MPolyElem{T}, b::AbstractAlgebra.MPolyElem{T}) where {T <: RingElement}
```

**Examples**

```jldoctest
julia> using AbstractAlgebra

julia> R,(x,y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> a = x*y + 2*y
x*y+2*y

julia> b = x^3*y + y
x^3*y+y

julia> gcd(a,b)
y

julia> lcm(a,b)
x^4*y+2*x^3*y+x*y+2*y

julia> lcm(a,b) == a * b // gcd(a,b)
true

```

### Derivations

```@docs
derivative(::AbstractAlgebra.MPolyElem{T}, ::AbstractAlgebra.MPolyElem{T}) where T <: AbstractAlgebra.RingElement
```

**Examples**

```jldoctest
julia> R, (x, y) = AbstractAlgebra.PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y])

julia> f = x*y + x + y + 1
x*y+x+y+1

julia> derivative(f, x)
y+1

julia> derivative(f, y)
x+1

```

### Homogeneous polynomials

It is possible to test whether a polynomial is homogeneous with respect to the standard grading using the function

```@docs
ishomogeneous(x::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```
