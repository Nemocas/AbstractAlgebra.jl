```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Universal polynomial

AbstractAlgebra.jl provides a module, implemented in `src/generic/UnivPoly.jl` for
a universal polynomial ring. This is very similar to the multivariate polynomial
rings, except that variables can be added to the ring at any time.

To compensate for the fact that the number of variables may change, many of the
functions relax their restrictions on exponent vectors. For example, if one
creates a polynomial when the ring only has two variables, each exponent vector
would consist of two integers. Later, when the ring has more variable, these
exponent vectors will still be accepted. The exponent vectors are simply padded
out to the full number of variables behind the scenes.

## Generic sparse distributed universal multivariable polynomial types

AbstractAlgebra provides a generic universal polynomial type `Generic.UnivPoly{T, U}`
where `T` is the type of elements of the coefficient ring and `U` is the type of
the elements of the underlying multivariate polynomial ring. Essentially, `U` can
be any type belonging to `MPolyElem{T}`.

Parent objects of such polynomials have type `Generic.UnivPolyRing{T, U}`.

## Abstract types

AbstractAlgebra also provides abstract types for universal polynomials and their
rings. These are `UnivPolyElem{T, U}` and `UnivPolyRing{T, U}` respectively.
These in turn belong to `Ring`.

## Polynomial ring constructors

In order to construct universal polynomials in AbstractAlgebra.jl, one must first
construct the universal polynomial ring itself. This is unique given a base ring.

The universal polynomial ring over a given base ring `R` is constructed with
one of the following constructor functions.

```julia
UniversalPolynomialRing(R::Ring; cached::Bool = true, ordering::Symbol=:lex)
```

Given a base ring `R` and an array `S` of strings, return an object representing
the universal polynomial ring $S = R[\ldots]$ with no variables in it initially.

**Examples**

```jldoctest
julia> S = UniversalPolynomialRing(ZZ)
Universal Polynomial Ring over Integers
```

## Adding variables

There are two ways to add variables to a universal polynomial ring `S`.

```julia
gen(S::UnivPolyRing, var::T) where T <: Union{Char, Symbol, String}
gens(S::UnivPolyRing, vars::Vector{T}) where T <: Union{Char, Symbol, String}
```

**Examples**

```jldoctest
julia> S = UniversalPolynomialRing(ZZ)
Universal Polynomial Ring over Integers

julia> x = gen(S, "x")
x

julia> y, z = gens(S, ["y", "z"])
(y, z)
```

## Universal polynomial functionality

The universal polynomial ring behaves exactly like a multivariate polynomial
ring with the few differences noted above.

The only functionality not implemented is the ability to do `divrem` by an
ideal of polynomials.

The universal polynomial ring is very useful for doing symbolic manipulation.
However, it is important to understand that AbstractAlgebra is not a symbolic
system and the performance of the universal polynomial ring will closely match
that of a multivariate polynomial ring with the same number of variables.

The disadvantage of this approach to symbolic manipulation is that some
manipulations that would be offered by a symbolic system are not available,
as variables are not identified by their names alone in AbstractAlgebra, as
would be the case symbolically, but by objects.

The most powerful symbolic tools we offer are the generalised evaluation
functions, the multivariate coefficient functionality, the ability to
change coefficient ring and to map coefficients according to a supplied
function and the ability to convert a multivariate which happens to have
just one variable into a dense univariate polynomial.

Further facilities may be added in future to ease symbolic manipulations.
