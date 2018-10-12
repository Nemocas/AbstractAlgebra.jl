###############################################################################
#
#   NCRings.jl : Generic not necessarily commutative rings
#
###############################################################################

elem_type(::T) where {T <: NCRing} = elem_type(T)

include("Rings.jl")

