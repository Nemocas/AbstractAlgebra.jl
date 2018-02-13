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
R, (x, y) = PolynomialRing(JuliaZZ, ["x", "y"]; ordering=:deglex)

f = R()
g = R(123)
h = R(BigInt(1234))
k = R(x + 1)
m = R(x + y + 1)
```

All of the examples here are generic polynomial rings, but specialised implementations
of polynomial rings provided by external modules will also usually provide a
`PolynomialRing` constructor to allow creation of their polynomial rings.

