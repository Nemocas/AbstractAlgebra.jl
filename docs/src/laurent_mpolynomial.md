```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Sparse distributed multivariate Laurent polynomials

Every element of the multivariate Laurent polynomial ring
$R[x_1, x_1^{-1}, \dots, x_n, x_n^{-1}]$ can be presented as a sum of products
of powers of the $x_i$ where the power can be *any* integer. Therefore, the
interface for sparse multivarate polynomials carries over with the additional
feature that exponents can be negative.

## Generic multivariate Laurent polynomial types

AbstractAlgebra.jl provides a generic implementation of multivariate Laurent
polynomials, built in terms of regular multivariate polynomials, in the file
`src/generic/LaurentMPoly.jl`.

The type `LaurentMPolyWrap{T, ...} <: LaurentMPolyRingElem{T}` implements generic
multivariate Laurent polynomials by wrapping regular polynomials:
a Laurent polynomial `l` wraps a polynomial `p` and a vector of integers $n_i$
such that $l = \prod_i x_i^{n_i} * p$. The representation is said to be
normalized when each $n_i$ is as large as possible (or zero when `l` is zero),
but the representation of a given element is not required to be normalized
internally.

The corresponding parent type is `LaurentMPolyWrapRing{T, ...} <: LaurentMPolyRing{T}`.

## Abstract types

Two abstract types `LaurentMPolyRingElem{T}` and `LaurentMPolyRing{T}`
are defined to represent Laurent polynomials and rings thereof, parameterized
on a base ring `T`.

## Multivate Laurent polynomial operations

Since, from the point of view of the interface, Laurent polynomials are simply
regular polynomials with possibly negative exponents, the following functions
from the polynomial interface are completely analogous. As with regular
polynomials, an implementation must provide access to the elements as a sum of
individual terms *in some order*. This order currently cannot be specified in
the constructor.

```julia
laurent_polynomial_ring(R::Ring, S::Vector{<:VarName}; cached::Bool = true)
laurent_polynomial_ring(R::Ring, n::Int, s::VarName; cached::Bool = false)
```

```julia
(S::LaurentMPolyRing{T})(A::Vector{T}, m::Vector{Vector{Int}})
```

```julia
MPolyBuildCtx(R::LaurentMPolyRing)
push_term!(M::LaurentMPolyBuildCtx, c::RingElem, v::Vector{Int})
finish(M::LaurentMPolyBuildCtx)
```

```julia
symbols(S::LaurentMPolyRing)
number_of_variables(f::LaurentMPolyRing)
gens(S::LaurentMPolyRing)
gen(S::LaurentMPolyRing, i::Int)
is_gen(x::LaurentMPolyRingElem)
var_index(p::LaurentMPolyRingElem)
length(f::LaurentMPolyRingElem)
```

```julia
coefficients(p::LaurentMPolyRingElem)
monomials(p::LaurentMPolyRingElem)
terms(p::LaurentMPolyRingElem)
exponent_vectors(p::LaurentMPolyRingElem)
leading_coefficient(p::LaurentMPolyRingElem)
leading_monomial(p::LaurentMPolyRingElem)
leading_term(p::LaurentMPolyRingElem)
leading_exponent_vector(p::LaurentMPolyRingElem)
```

```julia
change_base_ring(::Ring, p::LaurentMPolyRingElem)
change_coefficient_ring(::Ring, p::LaurentMPolyRingElem)
map_coefficients(::Any, p::LaurentMPolyRingElem)
```

```julia
evaluate(p::LaurentMPolyRingElem, ::Vector)
```

```julia
derivative(p::LaurentMPolyRingElem, x::LaurentMPolyRingElem)
derivative(p::LaurentMPolyRingElem, i::Int)
```

```julia
rand(R::LaurentMPolyRingElem, length_range::AbstractUnitRange{Int}, exp_range::AbstractUnitRange{Int}, v...)
```

The choice of canonical unit for Laurent polynomials includes the product
$\prod_i x_i^{n_i}$ from the normalized representation. In particular,
this means that the output of `gcd` will not have any negative exponents.

```jldoctest
julia> R, (x, y) = laurent_polynomial_ring(ZZ, ["x", "y"]);

julia> canonical_unit(2*x^-5 - 3*x + 4*y^-4 + 5*y^2)
-x^-5*y^-4

julia> gcd(x^-3 - y^3, x^-2 - y^2)
x*y - 1
```

