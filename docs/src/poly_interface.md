```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Univariate Polynomial Ring Interface

Univariate polynomial rings are supported in AbstractAlgebra, and in addition to the
standard Ring interface, numerous additional functions are required to be present for
univariate polynomial rings.

Univariate polynomial rings can be built over both commutative and noncommutative
rings.

Univariate polynomial rings over a field are also Euclidean and therefore such rings
must implement the Euclidean interface.

Since a sparse distributed multivariate format can generally also handle sparse
univariate polynomials, the univariate polynomial interface is designed around the
assumption that they are dense. This is not a requirement, but it may be easier to
use the multivariate interface for sparse univariate types.

## Types and parents

AbstractAlgebra provides two abstract types for polynomial rings and their elements
over a commutative ring:

  * `PolyRing{T}` is the abstract type for univariate polynomial ring parent types
  * `PolyRingElem{T}` is the abstract type for univariate polynomial types

Similarly there are two abstract types for polynomial rings and their elements
over a noncommutative ring:

  * `NCPolyRing{T}` is the abstract type for univariate polynomial ring parent types
  * `NCPolyRingElem{T}` is the abstract type for univariate polynomial types

We have that `PolyRing{T} <: Ring` and
`PolyRingElem{T} <: RingElem`. Similarly we have that
`NCPolyRing{T} <: NCRing` and 
`NCPolyRingElem{T} <: NCRingElem`.

Note that the abstract types are parameterised. The type `T` should usually be the type
of elements of the coefficient ring of the polynomial ring. For example, in the case of
$\mathbb{Z}[x]$ the type `T` would be the type of an integer, e.g. `BigInt`.

If the parent object for such a ring has type `MyZX` and polynomials in that ring have
type `MyZXPoly` then one would have:

  * `MyZX <: PolyRing{BigInt}`
  * `MyZXPoly <: PolyRingElem{BigInt}`

Polynomial rings should be made unique on the system by caching parent objects (unless
an optional `cache` parameter is set to `false`). Polynomial rings should at least be
distinguished based on their base (coefficient) ring. But if they have the same base
ring and symbol (for their variable/generator), they should certainly have the same
parent object.

See `src/generic/GenericTypes.jl` for an example of how to implement such a cache (which
usually makes use of a dictionary).

## Required functionality for univariate polynomials

In addition to the required functionality for the Ring/NCRing interface (and in the case
of polynomials over a field, the Euclidean Ring interface), the Polynomial Ring interface
has the following required functions.

We suppose that `R` is a fictitious base ring (coefficient ring) and that `S` is a
univariate polynomial ring over `R` (i.e. $S = R[x]$) with parent object `S` of type
`MyPolyRing{T}`. We also assume the polynomials in the ring have type `MyPoly{T}`, where
`T` is the type of elements of the base (coefficient) ring.

Of course, in practice these types may not be parameterised, but we use parameterised
types here to make the interface clearer.

Note that the type `T` must (transitively) belong to the abstract type `RingElem` or
`NCRingElem`.

We describe the functionality below for polynomials over commutative rings, i.e. with
element type belonging to `RingElem`, however similar constructors should be available
for element types belonging to `NCRingElem` instead, if the coefficient ring is
noncommutative.

### Constructors

In addition to the standard constructors, the following constructors, taking an array of
coefficients, must be available.

```julia
(S::MyPolyRing{T})(A::Vector{T}) where T <: RingElem
(S::MyPolyRing{T})(A::Vector{U}) where T <: RingElem, U <: RingElem
(S::MyPolyRing{T})(A::Vector{U}) where T <: RingElem, U <: Integer
```

Create the polynomial in the given ring whose degree $i$ coefficient is given by `A[1 + i]`.
The elements of the array are assumed to be able to be coerced into the base ring `R`.
If the argument is an empty vector, the zero polynomial shall be returned.

It may be desirable to have a additional version of the function that accepts an array
of Julia `Int` values  if this can be done more efficiently.

It is also possible to create polynomials directly without first creating the
corresponding polynomial ring.

```julia
polynomial(R::Ring, arr::Vector{T}, var::VarName=:x; cached::Bool=true)
```

Given an array of coefficients construct the polynomial with those coefficients
over the given ring and with the given variable.

!!! note

    If `cached` is set to `false` then the parent ring of the created
    polynomial is not cached. However, this means that subsequent polynomials
    created in the same way will not be compatible. Instead, one should use the
    parent object of the first polynomial to create subsequent polynomials instead
    of calling this function repeatedly with `cached=false`.

### Data type and parent object methods

```julia
var(S::MyPolyRing{T}) where T <: RingElem
```

Return a `Symbol` representing the variable (generator) of the polynomial ring. Note
that this is a `Symbol` not a `String`, though its string value will usually be used
when printing polynomials.

```julia
symbols(S::MyPolyRing{T}) where T <: RingElem
```

Return the array `[s]` where `s` is a `Symbol` representing the variable of the given
polynomial ring. This is provided for uniformity with the multivariate interface, where
there is more than one variable and hence an array of symbols.

```julia
dense_poly_type(::Type{T}) where T <: RingElement
```

Return the type of a polynomial whose coefficients have the given type. In our
example `MyPoly{T}`.

This function is defined for generic polynomials and only needs to be defined for
custom polynomial rings, e.g. ones defined by a C implementation.

```@docs
polynomial_ring_only(R::Ring, s::Symbol; cached::Bool=true)
```

The default implementation figures out the appropriate polynomial ring type via
`dense_poly_type` and calls its constructor with `R, s, cached` as arguments.
In our example, this would be

```julia
MyPolyRing{T}(R, s, cached)
```

Accordingly, `polynomial_ring_only` only needs to be defined, if such a
constructor does not exist or other behaviour is wanted.

### Basic manipulation of rings and elements

```julia
length(f::MyPoly{T}) where T <: RingElem
```

Return the length of the given polynomial. The length of the zero polynomial is defined
to be $0$, otherwise the length is the degree plus $1$. The return value should be of
type `Int`.

```julia
set_length!(f::MyPoly{T}, n::Int) where T <: RingElem
```

This function must zero any coefficients beyond the requested length $n$ and then set
the length of the polynomial to $n$. This function does not need to normalise the
polynomial and is not useful to the user, but is used extensively by the AbstractAlgebra
generic functionality.

This function returns the resulting polynomial.

```julia
coeff(f::MyPoly{T}, n::Int) where T <: RingElem
```

Return the coefficient of the polynomial `f` of degree `n`. If `n` is larger than the
degree of the polynomial, it should return zero in the coefficient ring. 

```julia
setcoeff!(f::MyPoly{T}, n::Int, a::T) where T <: RingElem
```

Set the degree $n$ coefficient of $f$ to $a$. This mutates the polynomial in-place
if possible and returns the mutated polynomial (so that immutable types can also be
supported). The function must not assume that the polynomial already has space for
$n + 1$ coefficients. The polynomial must be resized if this is not the case.

Note that this function is not required to normalise the polynomial and is not
necessarily useful to the user, but is used extensively by the generic functionality in
AbstractAlgebra.jl. It is for setting raw coefficients in the representation.

```julia
normalise(f::MyPoly{T}, n::Int) where T <: RingElem
```

Given a polynomial whose length is currently $n$, including any leading zero
coefficients, return the length of the normalised polynomial (either zero or the length
of the polynomial with nonzero leading coefficient). Note that the function does not
actually perform the normalisation.

```julia
fit!(f::MyPoly{T}, n::Int) where T <: RingElem
```

Ensure that the polynomial $f$ internally has space for $n$ coefficients. This function
must mutate the function in-place if it is mutable. It does not return the mutated
polynomial. Immutable types can still be supported by defining this function to do
nothing.

Some interfaces for C polynomial types automatically manage the internal allocation of
polynomials in every function that can be called on them. Explicit adjustment by
the generic code in AbstractAlgebra.jl is not required. In such cases, this function
can also be defined to do nothing.

## Optional functionality for polynomial rings

Sometimes parts of the Euclidean Ring interface can and should be implemented for
polynomials over a ring that is not necessarily a field.

When divisibility testing can be implemented for a polynomial ring over a field, it 
should be possible to implement the following functions from the Euclidean Ring
interface:

  * divides
  * remove
  * valuation

When the given polynomial ring is a GCD domain, with an effective GCD algorithm, it
may be possible to implement the following functions:

  * gcd
  * lcm

Polynomial rings can optionally implement any part of the generic univariate polynomial
functionality provided by AbstractAlgebra.jl, using the same interface. 

Obviously additional functionality can also be added to that provided by
AbstractAlgebra.jl on an ad hoc basis.

### Similar

The `similar` function is available for all univariate polynomial types, but
new polynomial rings can define a specialised version of it if required.

```julia
similar(x::MyPoly{T}, R::Ring=base_ring(x), var::VarName=var(parent(x))) where T <: RingElem
```

Construct the zero polynomial with the given variable and coefficients in the
given ring, if specified, and with the defaults shown if not.

Custom polynomial rings may choose which polynomial type is best-suited to
return for any given arguments. If they don't specialise the function the
default polynomial type returned is a `Generic.Poly`.

