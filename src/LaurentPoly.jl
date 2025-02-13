###############################################################################
#
#   LaurentPoly.jl : Laurent polynomials over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

characteristic(R::LaurentPolyRing) = characteristic(base_ring(R))
is_known(::typeof(characteristic), R::LaurentPolyRing) = is_known(characteristic, base_ring(R))

is_finite(R::LaurentPolyRing) = is_trivial(R)

###############################################################################
#
#   laurent_polynomial_ring constructor
#
###############################################################################

@doc raw"""
    laurent_polynomial_ring(R::Ring, s::VarName)

Given a base ring `R` and string `s` specifying how the generator (variable)
should be printed, return a tuple `S, x` representing the new Laurent polynomial
ring $S = R[x, 1/x]$ and the generator $x$ of the ring.

# Examples
```julia
julia> R, x = laurent_polynomial_ring(ZZ, :x)
(Univariate Laurent Polynomial Ring in x over Integers, x)

julia> 2x^-3 + x^2
x^2 + 2*x^-3

julia> rand(R, -3:3, -9:9)
-3*x^2 - 8*x + 4 + 3*x^-1 - 6*x^-2 + 9*x^-3
```
"""
function laurent_polynomial_ring(R::Ring, s::VarName; cached::Bool = true)
   @req !is_trivial(R) "Zero rings are currently not supported as coefficient ring."
   return Generic.laurent_polynomial_ring(R, Symbol(s); cached)
end
