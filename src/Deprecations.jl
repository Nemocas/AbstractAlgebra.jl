# Deprecated in 0.28.*

# we can't use @deprecate because this is also a typename
function ResidueRing(R::Ring, a::RingElement; cached::Bool = true)
    Base.depwarn("'ResidueRing(R::Ring, a::RingElement; cached::Bool = true)' is deprecated, use "*
    "'residue_ring(R::Ring, a::RingElement; cached::Bool = true)' instead.", :ResidueRing)
    return residue_ring(R, a; cached)
end

# we can't use @deprecate because this is also a typename
function ResidueField(R::Ring, a::RingElement; cached::Bool = true)
    Base.depwarn("'ResidueField(R::Ring, a::RingElement; cached::Bool = true)' is deprecated, use "*
    "'residue_field(R::Ring, a::RingElement; cached::Bool = true)' instead.", :ResidueField)
    return residue_field(R, a; cached)
end

# Deprecated in 0.30.*

@deprecate factor(f::FracElem, R::Ring) factor(R, f)

@deprecate factor(f::PolyRingElem, R::Field) factor(R, f)

@deprecate roots(f::PolyRingElem, R::Field) roots(R, f)

# deprecated in 0.32.x
@deprecate mat(f::Map(FPModuleHomomorphism)) matrix(f::Map(FPModuleHomomorphism))

# deprecated in 0.34.x
@deprecate power_series_ring(R::AbstractAlgebra.Ring, weights::Vector{Int},
    prec::Int, s::Vector{Symbol}; kw...) power_series_ring(R, prec, s; weights, kw...)
