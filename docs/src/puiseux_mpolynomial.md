```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Multivariate Puiseux polynomials

A Puiseux polynomial over a field $K$ is a finite sum of terms
$c x_1^{e_1} \cdots x_n^{e_n}$ with $c \in K$ and rational exponents $e_i$,
which may be negative. Puiseux polynomials form a ring containing the Laurent
polynomial ring $K[x_1^{\pm 1}, \dots, x_n^{\pm 1}]$.

## Generic Puiseux polynomial types

AbstractAlgebra.jl provides a generic implementation of multivariate Puiseux
polynomials in the file `src/generic/PuiseuxMPoly.jl`, built on top of the
generic multivariate Laurent polynomials.

The element type is `PuiseuxMPolyRingElem{T} <: RingElem` and the parent type
is `PuiseuxMPolyRing{T} <: Ring`, parameterized on the coefficient type `T`.

A Puiseux polynomial $f$ is stored as a pair $(g, d)$ of a Laurent polynomial
$g$ and a positive integer $d$, its *scale*, with

```math
f(x_1, \dots, x_n) = g(x_1^{1/d}, \dots, x_n^{1/d}).
```

The representation is normalized: $d$ is coprime to the greatest common divisor
of all exponents of $g$, so $d$ is the least common denominator of the exponents
of $f$. Every Puiseux polynomial has a unique normalized representation, which
makes comparing elements cheap. Arithmetic rescales both operands to a common
scale, operates on the underlying Laurent polynomials and normalizes the result.

## Constructor

```@docs
puiseux_polynomial_ring(::Field, ::Vector{String})
```

Elements can be created from the generators by the usual arithmetic
operations, and constants can be coerced into the ring by calling it:

```jldoctest
julia> R, (x, y) = puiseux_polynomial_ring(QQ, ["x", "y"]);

julia> R(3)
3

julia> 3*x^(1//2)*y^(-2//3) - 1
3*x^(1//2)*y^(-2//3) - 1
```

## Basic functionality

Puiseux polynomial rings implement the Ring interface. In addition, the
following functions are available.

```julia
base_ring(R::PuiseuxMPolyRing)
coefficient_ring(R::PuiseuxMPolyRing)
symbols(R::PuiseuxMPolyRing)
number_of_variables(R::PuiseuxMPolyRing)
gens(R::PuiseuxMPolyRing)
is_univariate(R::PuiseuxMPolyRing)
```

Here `base_ring` returns the underlying Laurent polynomial ring, and
`coefficient_ring` returns the field $K$.

```julia
length(f::PuiseuxMPolyRingElem)
coefficients(f::PuiseuxMPolyRingElem)
monomials(f::PuiseuxMPolyRingElem)
exponent_vectors(f::PuiseuxMPolyRingElem)
leading_coefficient(f::PuiseuxMPolyRingElem)
is_gen(f::PuiseuxMPolyRingElem)
is_monomial(f::PuiseuxMPolyRingElem)
is_term(f::PuiseuxMPolyRingElem)
is_unit(f::PuiseuxMPolyRingElem)
canonical_unit(f::PuiseuxMPolyRingElem)
```

The exponent vectors have entries of type `Rational{Int}`. The units are the
nonzero terms, i.e. the nonzero constant multiples of monomials. As for
Laurent polynomials, the canonical unit of `f` is the canonical unit of its
leading coefficient times $x_1^{m_1} \cdots x_n^{m_n}$, where $m_i$ is the
smallest exponent of $x_i$ in `f`.

```jldoctest
julia> R, (x, y) = puiseux_polynomial_ring(QQ, ["x", "y"]);

julia> f = x^(1//2) + 2*y^(1//3);

julia> collect(coefficients(f))
2-element Vector{Rational{BigInt}}:
 1
 2

julia> exponent_vectors(f)
2-element Vector{Vector{Rational{Int64}}}:
 [1//2, 0]
 [0, 1//3]

julia> is_unit(f), is_unit(2*x^(1//2))
(false, true)

julia> canonical_unit(3*x*y + 2*x^(1//2))
3*x^(1//2)
```

The normalized representation can be inspected with the following functions.

```@docs
poly(::PuiseuxMPolyRingElem)
scale(::PuiseuxMPolyRingElem)
```

## Powering

```julia
^(f::PuiseuxMPolyRingElem, n::Integer)
^(f::PuiseuxMPolyRingElem, q::Rational)
```

Any Puiseux polynomial can be raised to a nonnegative integer power. A negative
integer power requires `f` to be a monomial times a constant, and a non-integer
rational power requires `f` to be a monomial with coefficient one.

```jldoctest
julia> R, (x, y) = puiseux_polynomial_ring(QQ, ["x", "y"]);

julia> (x^(1//2) + y)^2
x + 2*x^(1//2)*y + y^2

julia> (x*y^2)^(1//3)
x^(1//3)*y^(2//3)

julia> (x + y)^(1//2)
ERROR: ArgumentError: only monomials can be exponentiated to non-integer powers
[...]
```

## Divisibility

```julia
divexact(f::PuiseuxMPolyRingElem, g::PuiseuxMPolyRingElem; check::Bool = true)
divexact(f::PuiseuxMPolyRingElem, a::RingElement; check::Bool = true)
divides(f::PuiseuxMPolyRingElem, g::PuiseuxMPolyRingElem)
```

`divexact` returns the quotient of `f` by a nonzero Puiseux polynomial or
constant `g`, provided the division is exact. `divides(f, g)` returns a tuple
`(flag, q)` where `flag` is `true` if $f = qg$ for some Puiseux polynomial $q$;
otherwise `q` is zero.

```jldoctest
julia> R, (x, y) = puiseux_polynomial_ring(QQ, ["x", "y"]);

julia> f = x^(1//2) + 2*y^(1//3);

julia> divexact(f^2, f)
x^(1//2) + 2*y^(1//3)

julia> divexact(f, x^(1//6))
x^(1//3) + 2*x^(-1//6)*y^(1//3)

julia> divides(f^2, f)
(true, x^(1//2) + 2*y^(1//3))

julia> divides(x, f)
(false, 0)
```

## Changing the coefficient ring

```julia
change_base_ring(K::Ring, f::PuiseuxMPolyRingElem; cached::Bool = true, parent::PuiseuxMPolyRing)
map_coefficients(g, f::PuiseuxMPolyRingElem; cached::Bool = true, parent::PuiseuxMPolyRing)
```

Return the Puiseux polynomial obtained from `f` by coercing each coefficient
into `K`, respectively by applying `g` to each coefficient. Unless `parent` is
given, the result lies in the Puiseux polynomial ring over `K`, respectively
over the parent of the images of `g`, in the same variables as `parent(f)`.

```jldoctest
julia> R, (x, y) = puiseux_polynomial_ring(QQ, ["x", "y"]);

julia> map_coefficients(c -> 2*c, x^(1//2) + 3*y^(1//3))
2*x^(1//2) + 6*y^(1//3)
```

## Valuation

```@docs
valuation(::PuiseuxMPolyRingElem)
```
