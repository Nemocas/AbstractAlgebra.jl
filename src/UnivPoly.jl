###############################################################################
#
#   UnivPoly.jl: Generic universal polynomial ring (variables can be added)
#
###############################################################################

function content(a::UniversalRingElem{<:MPolyRingElem})
   return content(data(a))
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function (a::UniversalPolyRingElem)(;kwargs...)
   ss = symbols(parent(a))
   vars = Int[]
   vals = RingElement[]
   for (var, val) in kwargs
     vari = findfirst(isequal(var), ss)
     vari === nothing && continue
     push!(vars, vari)
     push!(vals, val)
   end
   return evaluate(a, vars, vals)
end

###############################################################################
#
#   Iterators
#
###############################################################################

function coefficients(a::UniversalRingElem{<:MPolyRingElem})
   return Generic.UnivPolyCoeffs(a)
end

function exponent_vectors(a::UniversalRingElem{<:MPolyRingElem})
   return Generic.UnivPolyExponentVectors(a)
end

function monomials(a::UniversalRingElem{<:MPolyRingElem})
   return Generic.UnivPolyMonomials(a)
end

function terms(a::UniversalRingElem{<:MPolyRingElem})
   return Generic.UnivPolyTerms(a)
end

###############################################################################
#
#  Factorization
#
###############################################################################

function _wrap_factorization(f::Fac{<:MPolyRingElem}, S::UniversalRing{<:MPolyRingElem})
   res = Fac{elem_type(S)}()
   res.unit = Generic.UnivPoly(f.unit, S)
   for (fact, expo) in f
      mulpow!(res, Generic.UnivPoly(fact, S), expo)
   end
   return res
end

factor_squarefree(f::UniversalRingElem{<:MPolyRingElem}) = _wrap_factorization(factor_squarefree(data(f)), parent(f))

factor(f::UniversalRingElem{<:MPolyRingElem}) = _wrap_factorization(factor(data(f)), parent(f))

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
(Universal Polynomial Ring over Integers, AbstractAlgebra.Generic.UniversalRingElem{AbstractAlgebra.Generic.MPoly{BigInt}}[x, y])

julia> z = gen(S, :z)
z

julia> x*y - z
x*y - z

julia> S = universal_polynomial_ring(ZZ)
Universal Polynomial Ring over Integers

julia> x = gen(S, :x)
x

julia> y, z = gens(S, [:y, :z])
2-element Vector{AbstractAlgebra.Generic.UniversalRingElem{AbstractAlgebra.Generic.MPoly{BigInt}}}:
 y
 z

julia> x*y - z
x*y - z
```
"""
function universal_polynomial_ring(R::Ring, varnames::Vector{Symbol}; cached::Bool=true, internal_ordering::Symbol=:lex)
   @req !is_trivial(R) "Zero rings are currently not supported as coefficient ring."
   P = poly_ring(R, varnames; internal_ordering)
   T = elem_type(P)

   S = Generic.UniversalRing{T}(P; cached)
   return (S, gens(S))
end

function universal_polynomial_ring(R::Ring; cached::Bool=true, internal_ordering::Symbol=:lex)
   return universal_polynomial_ring(R, Symbol[]; internal_ordering, cached)[1]
end

@varnames_interface universal_polynomial_ring(R::Ring, s)
