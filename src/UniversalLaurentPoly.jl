###############################################################################
#
#   UniversalLaurentPoly.jl : Universal laurent polynomial ring (variables can be added)
#
###############################################################################

###############################################################################
#
#   String I/O
#
###############################################################################

universal_ring_name(R::UniversalRing{<:LaurentMPolyRingElem}) = "laurent polynomial ring"

###############################################################################
#
#   constructors
#
###############################################################################

@doc raw"""
    universal_laurent_polynomial_ring(R::Ring, varnames::Vector{Symbol}; cached::Bool=true)
    universal_laurent_polynomial_ring(R::Ring; cached::Bool=true)

Given a coefficient ring `R` and variable names, say `varnames = [:x1, :x2, ...]`, return
a tuple `S, [x1, x2, ...]` of the universal laurent polynomial ring `S = R[x1, x2, \dots]` and its generators `x1, x2, \dots`.

If `varnames` is omitted, return an object representing
the universal laurent polynomial ring `S = R[\ldots]` with no variables in it initially.

# Examples

```jldoctest
julia> S, (x,y) = universal_laurent_polynomial_ring(ZZ, [:x,:y])
(Universal laurent polynomial ring over Integers, UniversalRingElem{AbstractAlgebra.Generic.LaurentMPolyWrap{BigInt, AbstractAlgebra.Generic.MPoly{BigInt}, AbstractAlgebra.Generic.LaurentMPolyWrapRing{BigInt, AbstractAlgebra.Generic.MPolyRing{BigInt}}}, BigInt}[x, y])

julia> z = gen(S, :z)
z

julia> x*y - z
x*y - z

julia> S = universal_laurent_polynomial_ring(ZZ)
Universal laurent polynomial ring over Integers

julia> x = gen(S, :x)
x

julia> y, z = gens(S, [:y, :z])
2-element Vector{UniversalRingElem{AbstractAlgebra.Generic.LaurentMPolyWrap{BigInt, AbstractAlgebra.Generic.MPoly{BigInt}, AbstractAlgebra.Generic.LaurentMPolyWrapRing{BigInt, AbstractAlgebra.Generic.MPolyRing{BigInt}}}, BigInt}}:
 y
 z

julia> x^(-1)*y - z
-z + x^-1*y
```
"""
function universal_laurent_polynomial_ring(R::Ring, varnames::Vector{Symbol}; cached::Bool=true)
   @req !is_trivial(R) "Zero rings are currently not supported as coefficient ring."
   S = get_cached!(UniversalLaurentPolyRingID, R, cached) do
      P, _ = laurent_polynomial_ring(R, varnames; cached = false)
      T = elem_type(P)
      U = elem_type(coefficient_ring(P))
      UniversalRing{T, U}(P)
   end
   return (S, gens(S, varnames))
end

const UniversalLaurentPolyRingID = CacheDictType{Ring, Ring}()

function universal_laurent_polynomial_ring(R::Ring; cached::Bool=true)
   return universal_laurent_polynomial_ring(R, Symbol[]; cached)[1]
end

@varnames_interface universal_laurent_polynomial_ring(R::Ring, s)
