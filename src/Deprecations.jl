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
