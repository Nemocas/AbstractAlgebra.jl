###############################################################################
#
#   DirectSumModule.jl : Direct sums of modules
#
###############################################################################

export DirectSumModule, DirectSumModuleElem, summands

###############################################################################
#
#   DirectSum constructor
#
###############################################################################

@doc Markdown.doc"""
    DirectSum(m::Vector{<:FPModule{T}}) where T <: RingElement

Return a tuple $M, f, g$ consisting of $M$ the direct sum of the modules `m`
(supplied as a vector of modules), a vector $f$ of the injections
of the $m[i]$ into $M$ and a vector $g$ of the projections from
$M$ onto the $m[i]$.
"""
function DirectSum(m::Vector{<:FPModule{T}}) where T <: RingElement
   return Generic.DirectSum(m)
end

function DirectSum(vals::FPModule{T}...) where T <: RingElement
   return DirectSum([vals...])
end

function direct_sum(m::Vector{<:Module{T}}) where T <: RingElement
   Generic.DirectSum(m)
end

function direct_sum(m::Module{T}...) where T <: RingElement
   Generic.DirectSum(m...)
end
