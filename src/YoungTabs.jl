export Partition, AllParts, YoungTableau, SkewDiagram
export rowlength, collength, hooklength, dim, isrimhook, leglength

##############################################################################
#
#   Partition constructor
#
##############################################################################

function Partition(part::AbstractVector{T}, check::Bool=true) where T
   Generic.Partition(part, check)
end
 
function AllParts(n::T) where T
   Generic.AllParts(n)
end

##############################################################################
#
#   SkewDiagram constructor
#
##############################################################################

function SkewDiagram(lambda::Vector{T}, mu::Vector{T}) where T
   Generic.SkewDiagram(lambda, mu)
end

# Also see src/AbstractAlgebra.jl for an additional constructor

##############################################################################
#
#   YoungTableau constructor
#
##############################################################################
 
function YoungTableau(p::Vector{Int})
   Generic.YoungTableau(p)
end

# Also see src/AbstractAlgebra.jl for additional constructors