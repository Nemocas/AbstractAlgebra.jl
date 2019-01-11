```@meta
CurrentModule = AbstractAlgebra
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

The string representation of the variables of the polynomial ring, the
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
supported options are `:lex`, `:deglex` and `:degrevlex	.

Here are some examples of creating multivariate polynomial rings and making use of the
resulting parent objects to coerce various elements into the polynomial ring.

**Examples**

```julia
R, (x, y) = PolynomialRing(ZZ, ["x", "y"]; ordering=:deglex)

f = R()
g = R(123)
h = R(BigInt(1234))
k = R(x + 1)
m = R(x + y + 1)
derivative(k,x)
derivative(k,y)
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
var_index(p::AbstractAlgebra.MPolyElem{T}) where T <: RingElement
```

**Examples**

```julia
R, (x, y) = PolynomialRing(ZZ, ["x", "y"])

f = x^2 + 2x + 1

V = vars(f)
var_index(y) == 2
```

### Changing base (coefficient) rings

In order to substitute the variables of a polynomial $f$ over a ring $T$ by
elements in a $T$-algebra $S$, you first have to change the base ring of $f$
using the following function, where $g$ is a function representing the
structure homomorphism of the $T$-algebra $S$.

```@docs
change_base_ring(p::AbstractAlgebra.MPolyElem{T}, g) where {T <: RingElement}
```

Note that $g$ can also be a Nemo parent, e.g. `QQ`.

**Examples**

```julia
R, (x, y) = PolynomialRing(ZZ, ["x", "y"])

fz = x^2*y^2 + x + 1

fq = change_base_ring(fz, QQ)
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

```julia
R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])

f = x^4*y^2*z^2 - 2x^4*y*z^2 + 4x^4*z^2 + 2x^2*y^2 + x + 1

coeff(f, [1, 3], [4, 2]) == coeff(f, [x, z], [4, 2])
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

```julia
R, (x, y) = PolynomialRing(ZZ, ["x", "y"])

f = x^7*y^8 + 3*x^4*y^8 - x^4*y^2 + 5x*y^5 - x*y^2

def, shift = deflation(f)
f1 = deflate(f, def, shift)
f2 = inflate(f1, def, shift)
f2 == f
```

### Conversions

```@docs
to_univariate(R::AbstractAlgebra.PolyRing{T}, p::AbstractAlgebra.MPolyElem{T}) where T <: AbstractAlgebra.RingElement
```

**Examples**

```julia
R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
S, z = PolynomialRing(ZZ, "z")

f = 2x^5 + 3x^4 - 2x^2 - 1

g = to_univariate(S, f)
```

