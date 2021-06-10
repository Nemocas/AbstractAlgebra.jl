###############################################################################
#
#   LaurentPoly.jl : Laurent polynomials over rings
#
###############################################################################

export LaurentPolynomialRing

###############################################################################
#
#   LaurentPolynomialRing constructor
#
###############################################################################

@doc doc"""
    LaurentPolynomialRing(R::Ring, s::Union{AbstractString, Char, Symbol})

Given a base ring `R` and string `s` specifying how the generator (variable)
should be printed, return a tuple `S, x` representing the new Laurent polynomial
ring $S = R[x, 1/x]$ and the generator $x$ of the ring.

## Examples
```julia
julia> R, x = LaurentPolynomialRing(ZZ, "x")
(Univariate Laurent Polynomial Ring in x over Integers, x)

julia> 2x^-3 + x^2
x^2 + 2*x^-3

julia> rand(R, -3:3, -9:9)
-3*x^2 - 8*x + 4 + 3*x^-1 - 6*x^-2 + 9*x^-3
```
"""
LaurentPolynomialRing(R::Ring, s::Union{AbstractString, Char, Symbol})

function LaurentPolynomialRing(R::Ring, s::Symbol)
   return Generic.LaurentPolynomialRing(R, s)
end

function LaurentPolynomialRing(R::Ring, s::Char)
   return LaurentPolynomialRing(R, Symbol(s))
end

function LaurentPolynomialRing(R::Ring, s::AbstractString)
   return LaurentPolynomialRing(R, Symbol(s))
end

