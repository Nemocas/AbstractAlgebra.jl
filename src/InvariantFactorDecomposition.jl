###############################################################################
#
#   SNFModule.jl : Invariant factor decomposition of modules
#
###############################################################################

###############################################################################
#
#   SNFModule constructor
#
###############################################################################

@doc raw"""
    invariant_factors(m::FPModule{T}) where T <: RingElement

Return a vector of the invariant factors of the module $M$.
"""
function invariant_factors(m::FPModule{T}) where T <: RingElement
   return Generic.invariant_factors(m)
end

