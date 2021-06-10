###############################################################################
#
#   FreeModule.jl : Free modules over rings
#
###############################################################################

export FreeModule, VectorSpace, free_module, vector_space

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

###############################################################################
#
#   VectorSpace constructor
#
###############################################################################

@doc Markdown.doc"""
    VectorSpace(R::Field, dim::Int; cached::Bool = true)

Return the vector space over the field $R$ with the given dimension.
"""
function VectorSpace(R::Field, dim::Int; cached::Bool = true)
   Generic.FreeModule(R, dim; cached=cached)
end

function vector_space(R::Field, dim::Int; cached::Bool = true)
   Generic.FreeModule(R, dim; cached=cached)
end
