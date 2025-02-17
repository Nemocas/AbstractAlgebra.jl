###############################################################################
#
#   UnivPoly.jl: Generic universal polynomial ring (variables can be added)
#
###############################################################################

function content(a::UniversalPolyRingElem)
   return content(data(a))
end

###############################################################################
#
#   Iterators
#
###############################################################################

function coefficients(a::UniversalPolyRingElem)
   return Generic.UnivPolyCoeffs(a)
end

function exponent_vectors(a::UniversalPolyRingElem)
   return Generic.UnivPolyExponentVectors(a)
end

function monomials(a::UniversalPolyRingElem)
   return Generic.UnivPolyMonomials(a)
end

function terms(a::UniversalPolyRingElem)
   return Generic.UnivPolyTerms(a)
end

###############################################################################
#
#   universal_polynomial_ring constructors
#
###############################################################################

@doc raw"""
    universal_polynomial_ring(R::Ring, varnames::Vector{Symbol}; cached::Bool=true, internal_ordering::Symbol=:lex)
    universal_polynomial_ring(R::Ring; cached::Bool=true, internal_ordering::Symbol=:lex)

Given a coefficient ring `R` and variable names, say `varnames = [:x1, :x2, ...]`, return
a tuple `S, [x1, x2, ...]` of the universal polynomial ring `S = R[x1, x2, \dots]` and its generators `x1, x2, \dots`.

If `varnames` is omitted, return an object representing
the universal polynomial ring `S = R[\ldots]` with no variables in it initially.

# Examples

```jldoctest
julia> S, (x,y) = universal_polynomial_ring(ZZ, [:x,:y])
(Universal Polynomial Ring over Integers, AbstractAlgebra.Generic.UnivPoly{BigInt}[x, y])

julia> z = gen(S, :z)
z

julia> x*y - z
x*y - z

julia> S = universal_polynomial_ring(ZZ)
Universal Polynomial Ring over Integers

julia> x = gen(S, :x)
x

julia> y, z = gens(S, [:y, :z])
(y, z)

julia> x*y - z
x*y - z
```
"""
function universal_polynomial_ring(R::Ring, varnames::Vector{Symbol}; cached::Bool=true, internal_ordering::Symbol=:lex)
   @req !is_trivial(R) "Zero rings are currently not supported as coefficient ring."
   T = elem_type(R)
   U = mpoly_type(R)

   S = Generic.UniversalPolyRing{T}(R, varnames, internal_ordering, cached)
   return (S, gens(S))
end

function universal_polynomial_ring(R::Ring; cached::Bool=true, internal_ordering::Symbol=:lex)
   return universal_polynomial_ring(R, Symbol[]; internal_ordering, cached)[1]
end

@varnames_interface universal_polynomial_ring(R::Ring, s)
