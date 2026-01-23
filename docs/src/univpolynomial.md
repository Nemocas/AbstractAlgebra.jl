```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
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
be any type belonging to `MPolyRingElem{T}`.

Parent objects of such polynomials have type `Generic.UniversalPolyRing{T, U}`.

## Abstract types

AbstractAlgebra also provides abstract types for universal polynomials and their
rings. These are `UniversalPolyRingElem{T, U}` and `UniversalPolyRing{T, U}` respectively.
These in turn belong to `Ring`.

## Polynomial ring constructors

In order to construct universal polynomials in AbstractAlgebra.jl, one must first
construct the universal polynomial ring itself. This is unique given a base ring.

The universal polynomial ring over a given base ring `R` is constructed with
one of the following constructor functions.

```@docs
universal_polynomial_ring
```

## Adding variables

There are two ways to add variables to a universal polynomial ring `S`.

```julia
gen(S::UniversalPolyRing, var::VarName)
gens(S::UniversalPolyRing, vars::Vector{VarName})
```

**Examples**

```jldoctest
julia> S = universal_polynomial_ring(ZZ)
Universal polynomial ring over Integers

julia> x = gen(S, :x)
x

julia> number_of_generators(S)
1

julia> y, z = gens(S, [:y, :z])
2-element Vector{UniversalRingElem{AbstractAlgebra.Generic.MPoly{BigInt}}}:
 y
 z

julia> number_of_generators(S)
3
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
