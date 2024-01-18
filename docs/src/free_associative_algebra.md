```@meta
CurrentModule = AbstractAlgebra.Generic
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Free algebras

AbstractAlgebra.jl provides a module, implemented in `src/FreeAssAlgebra.jl` for
free associative algebras over any commutative ring belonging to the
AbstractAlgebra abstract type hierarchy.

## Generic free algebra types

AbstractAlgebra provides a generic type `Generic.FreeAssAlgElem{T}`
where `T` is the type of elements of the coefficient ring. The elements are
implemented using a Julia array of coefficients and a vector of
vectors of `Int`s for the monomial words. Parent objects of such elements have
type `Generic.FreeAssAlgebra{T}`.

The element types belong to the abstract type `NCRingElem`,
and the algebra types belong to the abstract type `NCRing`.

The following basic functions are implemented.
```julia
base_ring(R::FreeAssAlgebra)
base_ring(a::FreeAssAlgElem)
parent(a::FreeAssAlgElem)
characteristic(R::FreeAssAlgebra)
```

## Free algebra constructors

```julia
free_associative_algebra(R::Ring, s::AbstractVector{<:VarName}; cached::Bool = true)
free_associative_algebra(R::Ring, n::Int, s::VarName; cached::Bool = false)
```

The first constructor, given a base ring `R` and an array `s` of variables,
will return a tuple `S, (x, ...)` representing the new algebra
$S = R \left<x, \ldots \right>$ and a tuple of generators $(x, ...)$.

The second constructor given a string `s` and a number of variables `n` will
do the same as the first constructor except that the variables will be
automatically numbered as, `s1`, `s2`, ..., `sn`.

By default the parent object `S` will depend only on `R` and  `(x, ...)` and
will be cached. Setting the optional argument `cached` to `false` will prevent
the parent object `S` from being cached.

**Examples**

```jldoctest
julia> R, (x, y) = free_associative_algebra(ZZ, ["x", "y"])
(Free associative algebra on 2 indeterminates over integers, AbstractAlgebra.Generic.FreeAssAlgElem{BigInt}[x, y])

julia> (x + y + 1)^2
x^2 + x*y + y*x + y^2 + 2*x + 2*y + 1


julia> (x*y*x*x)^4
x*y*x^3*y*x^3*y*x^3*y*x^2
```

## Free algebra element constructors

Elements of a free algebra can be constructed from the generators in the
usual way using arithmetic operations. Also, all of the standard ring element
constructors may be used. Finally, the `MPolyBuildCtx` is overloaded to work
with coefficients and monomial words and not exponent vectors.

**Examples**

```jldoctest
julia> R, (x, y, z) = free_associative_algebra(ZZ, ["x", "y", "z"])
(Free associative algebra on 3 indeterminates over integers, AbstractAlgebra.Generic.FreeAssAlgElem{BigInt}[x, y, z])

julia> B = MPolyBuildCtx(R)
Builder for an element of Free associative algebra on 3 indeterminates over integers

julia> push_term!(B, ZZ(1), [1,2,3,1]); push_term!(B, ZZ(2), [3,3,1]); finish(B)
x*y*z*x + 2*z^2*x

julia> push_term!(B, ZZ(3), [3,3,3]); push_term!(B, ZZ(4), Int[]); finish(B)
3*z^3 + 4

julia> [gen(R, 2), R(9)]
2-element Vector{AbstractAlgebra.Generic.FreeAssAlgElem{BigInt}}:
 y
 9
```

## Element functions

### Basic manipulation

The standard ring functions are available. The following functions from the
multivariate polynomial interface are provided.

```julia
symbols(S::FreeAssAlgebra)
number_of_variables(f::FreeAssAlgebra)
gens(S::FreeAssAlgebra)
gen(S::FreeAssAlgebra, i::Int)
is_gen(x::FreeAssAlgElem)
total_degree(a::FreeAssAlgElem)
length(f::FreeAssAlgElem)
```

As with multivariate polynomials, an implementation must provide access to
the elements as a sum of individual terms *in some order*. The `length`
function provides the number of such terms, and the following functions
provide the first such term.


```julia
leading_coefficient(a::FreeAssAlgElem)
leading_monomial(a::FreeAssAlgElem)
leading_term(a::FreeAssAlgElem)
leading_exponent_word(a::FreeAssAlgElem)
```

For types that allow constant time access to coefficients, the following are
also available, allowing access to the given coefficient, monomial or term.
Terms are numbered from the most significant first.

```julia
coeff(f::FreeAssAlgElem, n::Int)
monomial(f::FreeAssAlgElem, n::Int)
term(f::FreeAssAlgElem, n::Int)
```

In contrast with the interface for multivariable polynomials, the function
`exponent_vector` is replaced by `exponent_word`

```@docs
exponent_word(a::Generic.FreeAssAlgElem{T}, i::Int) where T <: RingElement
```

**Examples**

```jldoctest
julia> R, (x, y, z) = free_associative_algebra(ZZ, ["x", "y", "z"])
(Free associative algebra on 3 indeterminates over integers, AbstractAlgebra.Generic.FreeAssAlgElem{BigInt}[x, y, z])

julia> map(total_degree, (R(0), R(1), -x^2*y^2*z^2*x + z*y))
(-1, 0, 7)

julia> leading_term(-x^2*y^2*z^2*x + z*y)
-x^2*y^2*z^2*x

julia> leading_monomial(-x^2*y^2*z^2*x + z*y)
x^2*y^2*z^2*x

julia> leading_coefficient(-x^2*y^2*z^2*x + z*y)
-1

julia> exponent_word(-x^2*y^2*z^2*x + z*y, 1)
7-element Vector{Int64}:
 1
 1
 2
 2
 3
 3
 1
```

```@docs
evaluate(a::AbstractAlgebra.FreeAssAlgElem{T}, vals::Vector{U}) where {T <: RingElement, U <: NCRingElem}
```

### Iterators

The following iterators are provided for elements of a free associative algebra,
with `exponent_words` providing the analogous functionality that `exponent_vectors`
provides for multivariate polynomials.

```julia
terms(p::FreeAssAlgElem)
coefficients(p::FreeAssAlgElem)
monomials(p::FreeAssAlgElem)
```

```@docs
exponent_words(a::FreeAssAlgElem{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> R, (a, b, c) = free_associative_algebra(ZZ, ["a", "b", "c"])
(Free associative algebra on 3 indeterminates over integers, AbstractAlgebra.Generic.FreeAssAlgElem{BigInt}[a, b, c])

julia> collect(terms(3*b*a*c - b + c + 2))
4-element Vector{Any}:
 3*b*a*c
 -b
 c
 2

julia> collect(coefficients(3*b*a*c - b + c + 2))
4-element Vector{Any}:
  3
 -1
  1
  2

julia> collect(monomials(3*b*a*c - b + c + 2))
4-element Vector{Any}:
 b*a*c
 b
 c
 1

julia> collect(exponent_words(3*b*a*c - b + c + 2))
4-element Vector{Vector{Int64}}:
 [2, 1, 3]
 [2]
 [3]
 []
```

### Groebner bases

The function `groebner_basis` provides the computation of a Groebner basis of an ideal, given a set of 
generators of that ideal.
Since such a Groebner basis is not necessarily finite, one can additionally pass a `reduction_bound`
to the function, to only compute a partial Groebner basis.

**Examples**

```jldoctest; setup = :(using AbstractAlgebra)
julia> R, (x, y, u, v, t, s) = free_associative_algebra(GF(2), ["x", "y", "u", "v", "t", "s"])
(Free associative algebra on 6 indeterminates over finite field F_2, AbstractAlgebra.Generic.FreeAssAlgElem{AbstractAlgebra.GFElem{Int64}}[x, y, u, v, t, s])

julia> g = Generic.groebner_basis([u*(x*y)^3 + u*(x*y)^2 + u + v, (y*x)^3*t + (y*x)^2*t + t + s])
5-element Vector{AbstractAlgebra.Generic.FreeAssAlgElem{AbstractAlgebra.GFElem{Int64}}}:
 u*x*y*x*y*x*y + u*x*y*x*y + u + v
 y*x*y*x*y*x*t + y*x*y*x*t + t + s
 u*x*s + v*x*t
 u*x*y*x*s + v*x*y*x*t
 u*x*y*x*y*x*s + v*x*y*x*y*x*t
```

In order to check whether a given element of the algebra is in the ideal generated by a Groebner 
basis `g`, one can compute its normal form.
```jldoctest; setup = :(using AbstractAlgebra)
julia> R, (x, y, u, v, t, s) = free_associative_algebra(GF(2), ["x", "y", "u", "v", "t", "s"]);

julia> g = Generic.groebner_basis([u*(x*y)^3 + u*(x*y)^2 + u + v, (y*x)^3*t + (y*x)^2*t + t + s]);

julia> normal_form(u*(x*y)^3*s*t + u*(x*y)^2*s*t +u*s*t + v*s*t, g)
0
 ```
