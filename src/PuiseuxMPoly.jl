#################################################################################
#
# Constructors
#
#################################################################################

@doc raw"""
    puiseux_polynomial_ring(K::Field, variableName::Vector{String})

Return a tuple `(R, x)` consisting of the ring `R` of Puiseux polynomials over
`K` in the given variables and the vector `x` of its generators.

The elements of `R` are finite sums of terms $c x_1^{e_1} \cdots x_n^{e_n}$ with
$c \in K$ and rational, possibly negative, exponents $e_i$.

# Examples

```jldoctest
julia> R, (x, y) = puiseux_polynomial_ring(QQ, ["x", "y"]);

julia> R
Puiseux polynomial ring in 2 variables x, y
  over rationals

julia> f = x^(1//2) + 2*y^(-1//3)
x^(1//2) + 2*y^(-1//3)

julia> f^2
x + 4*x^(1//2)*y^(-1//3) + 4*y^(-2//3)
```
"""
function puiseux_polynomial_ring(K::Field, variableSymbols::Vector{Symbol})
    base_ring, _ = laurent_polynomial_ring(K, variableSymbols)
    Kt = Generic.PuiseuxMPolyRing(base_ring)
    return Kt, gens(Kt)
end

@varnames_interface Generic.puiseux_polynomial_ring(R::Ring, s)
