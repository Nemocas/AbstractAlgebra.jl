###############################################################################
#
#   FreeModule.jl : Free modules over rings
#
###############################################################################

export FreeModule, FreeModuleElem

###############################################################################
#
#   FreeModule constructor
#
###############################################################################

@doc Markdown.doc"""
    FreeModule(R::NCRing, rank::Int; cached::Bool = true)

Return the free module over the ring $R$ with the given rank.
"""
function FreeModule(R::NCRing, rank::Int; cached::Bool = true)
   return Generic.FreeModule(R, rank, cached=cached)
end

function free_module(R::NCRing, rank::Int; cached::Bool = true)
   return Generic.FreeModule(R, rank; cached=cached)
end

