```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Generic Laurent polynomials

Laurent polynomials are similar to polynomials but can have terms of negative degrees,
and form a ring denoted by $R[x, x^{-1}]$ where `R`is the coefficient ring.
AbstractAlgebra.jl provides a generic implementation of Laurent polynomials,
built in terms of regular polynomials.

## Types and parent objects

Two abstract types `LaurentPolyElem{T}` and `LaurentPolynomialRing{T}`
are defined to represent Laurent polynomials and rings thereof, parameterized
on a base ring `T`.

The type `LaurentPolyWrap{T, ...} <: LaurentPolyElem{T}` implements generic Laurent
polynomials by wrapping regular polynomials: a Laurent polynomial `l` wraps a
polynomial `p` and an integer `n` such that $l = x^{-n} * p$.

The corresponding parent type is `LaurentPolyWrapRing{T, ...} <: LaurentPolynomialRing{T}`.

## Laurent polynomials ring constructor

In order to instantiate Laurent polynomials, one must first construct the parent ring:

```@docs
LaurentPolynomialRing
```

## Basic functionality

Laurent polynomials implement the ring interface (except for `divexact`), and some methods
from the polynomial interface, for example:

```jldoctest
julia> R, x = LaurentPolynomialRing(ZZ, "x")
(Univariate Laurent Polynomial Ring in x over Integers, x)

julia> var(R)
:x

julia> symbols(R)
1-element Array{Symbol,1}:
 :x

julia> nvars(R)
1

julia> f = x^-2 + 2x
2*x + x^-2

julia> coeff.(f, -2:2)
5-element Array{BigInt,1}:
 1
 0
 0
 2
 0

julia> setcoeff!(f, 3, big(5))
5*x^3 + 2*x + x^-2

julia> isgen(f)
false

julia> shift_left(f,2)
5*x^5 + 2*x^3 + 1

julia> map_coeffs(x->2x, f)
10*x^3 + 4*x + 2*x^-2

julia> change_base_ring(RealField, f)
5.0*x^3 + 2.0*x + x^-2

julia> lead(f), trail(f)
(5, 1)
```
