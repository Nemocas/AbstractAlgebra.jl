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
#   Data type and parent object methods
#
###############################################################################

@doc raw"""
    universal_poly_type(::Type{T}) where T<:RingElement
    universal_poly_type(::T) where T<:RingElement
    universal_poly_type(::Type{S}) where S<:Ring
    universal_poly_type(::S) where S<:Ring

The type of universal polynomials with coefficients of type `T` respectively `elem_type(S)`.
Falls back to `Generic.UnivPoly{T}`.

See also [`universal_poly_ring_type`](@ref), [`mpoly_type`](@ref) and [`mpoly_ring_type`](@ref).
"""
universal_poly_type(::Type{T}) where T<:RingElement = Generic.UnivPoly{T}
universal_poly_type(::Type{S}) where S<:Ring = universal_poly_type(elem_type(S))
universal_poly_type(x) = universal_poly_type(typeof(x)) # to stop this method from eternally recursing on itself, we better add ...
universal_poly_type(::Type{T}) where T = throw(ArgumentError("Type `$T` must be subtype of `RingElement`."))
universal_poly_type(T::Type{Union{}}) = throw(MethodError(universal_poly_type, (T,)))

@doc raw"""
    universal_poly_ring_type(::Type{T}) where T<:RingElement
    universal_poly_ring_type(::T) where T<:RingElement
    universal_poly_ring_type(::Type{S}) where S<:Ring
    universal_poly_ring_type(::S) where S<:Ring

The type of universal polynomial rings with coefficients of type `T`
respectively `elem_type(S)`. Implemented via [`universal_poly_type`](@ref).

See also [`mpoly_type`](@ref) and [`mpoly_ring_type`](@ref).
"""
universal_poly_ring_type(x) = parent_type(universal_poly_type(x))

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
#  Factorization
#
###############################################################################

function _wrap_factorization(f::Fac{<:MPolyRingElem}, S::UniversalPolyRing)
   res = Fac{elem_type(S)}()
   res.unit = Generic.UnivPoly(f.unit, S)
   for (fact, expo) in f
      mulpow!(res, Generic.UnivPoly(fact, S), expo)
   end
   return res
end

factor_squarefree(f::UniversalPolyRingElem) = _wrap_factorization(factor_squarefree(data(f)), parent(f))

factor(f::UniversalPolyRingElem) = _wrap_factorization(factor(data(f)), parent(f))

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
2-element Vector{AbstractAlgebra.Generic.UnivPoly{BigInt}}:
 y
 z

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
