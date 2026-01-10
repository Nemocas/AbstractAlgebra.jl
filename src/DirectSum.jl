###############################################################################
#
#   DirectSumModule.jl : Direct sums of modules
#
###############################################################################

###############################################################################
#
#   DirectSum constructor
#
###############################################################################

@doc raw"""
    direct_sum(m::Vector{<:FPModule{T}}) where T <: RingElement
    direct_sum(vals::FPModule{T}...) where T <: RingElement

Return a tuple $M, f, g$ consisting of $M$ the direct sum of the modules `m`
(supplied as a vector of modules), a vector $f$ of the injections
of the $m[i]$ into $M$ and a vector $g$ of the projections from
$M$ onto the $m[i]$.
"""
direct_sum(m::Vector{<:FPModule{<:RingElement}}) = Generic.direct_sum(m)
direct_sum(m::FPModule{<:RingElement}...) = Generic.direct_sum([m...])
