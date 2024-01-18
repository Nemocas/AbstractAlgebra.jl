```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Generic Laurent polynomials

Laurent polynomials are similar to polynomials but can have terms of negative degrees,
and form a ring denoted by $R[x, x^{-1}]$ where `R` is the coefficient ring.

## Generic Laurent polynomial types

AbstractAlgebra.jl provides a generic implementation of Laurent polynomials,
built in terms of regular polynomials in the file `src/generic/LaurentPoly.jl`.

The type `LaurentPolyWrap{T, ...} <: LaurentPolyRingElem{T}` implements generic Laurent
polynomials by wrapping regular polynomials: a Laurent polynomial `l` wraps a
polynomial `p` and an integer `n` such that $l = x^{-n} * p$.

The corresponding parent type is `LaurentPolyWrapRing{T, ...} <: LaurentPolyRing{T}`.

## Abstract types

Two abstract types `LaurentPolyRingElem{T}` and `LaurentPolyRing{T}`
are defined to represent Laurent polynomials and rings thereof, parameterized
on a base ring `T`.

## Laurent polynomials ring constructor

In order to instantiate Laurent polynomials, one must first construct the parent ring:

```@docs
laurent_polynomial_ring
```

## Basic functionality

Laurent polynomials implement the ring interface, and some methods
from the polynomial interface, for example:

```jldoctest
julia> R, x = laurent_polynomial_ring(ZZ, "x")
(Univariate Laurent polynomial ring in x over integers, x)

julia> var(R)
:x

julia> symbols(R)
1-element Vector{Symbol}:
 :x

julia> number_of_variables(R)
1

julia> f = x^-2 + 2x
2*x + x^-2

julia> coeff.(f, -2:2)
5-element Vector{BigInt}:
 1
 0
 0
 2
 0

julia> set_coefficient!(f, 3, ZZ(5))
5*x^3 + 2*x + x^-2

julia> is_gen(f)
false

julia> shift_left(f,2)
5*x^5 + 2*x^3 + 1

julia> map_coefficients(x->2x, f)
10*x^3 + 4*x + 2*x^-2

julia> change_base_ring(RealField, f)
5.0*x^3 + 2.0*x + x^-2

julia> leading_coefficient(f), trailing_coefficient(f)
(5, 1)
```
