#################################################################################
#
# Constructors
#
#################################################################################

function puiseux_polynomial_ring(K::Field, variableSymbols::Vector{Symbol})
    base_ring, _ = laurent_polynomial_ring(K, variableSymbols)
    Kt = Generic.PuiseuxMPolyRing(base_ring)
    return Kt, gens(Kt)
end

@varnames_interface Generic.puiseux_polynomial_ring(R::Ring, s)
