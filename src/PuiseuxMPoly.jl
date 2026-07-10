#################################################################################
#
# Constructors
#
#################################################################################

function puiseux_polynomial_ring(K::Field, variableName::Vector{String})
    #@req !isempty(variableName) "list of variables must not be empty"
    base_ring, _ = laurent_polynomial_ring(K, variableName)
    Kt = Generic.PuiseuxMPolyRing(base_ring)
    return Kt, gens(Kt)
end
