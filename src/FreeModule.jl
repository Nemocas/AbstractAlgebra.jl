###############################################################################
#
#   FreeModule.jl : Free modules over rings
#
###############################################################################

###############################################################################
#
#   FreeModule constructor
#
###############################################################################

@doc raw"""
    free_module(R::NCRing, rank::Int; cached::Bool = true)

Return the free module over the ring $R$ with the given rank.
"""
function free_module(R::NCRing, rank::Int; cached::Bool = true)
  return Generic.FreeModule(R, rank; cached=cached)
end

###############################################################################
#
#   VectorSpace constructor
#
###############################################################################

@doc raw"""
    vector_space(R::Field, dim::Int; cached::Bool = true)

Return the vector space over the field $R$ with the given dimension.
"""
function vector_space(R::Field, dim::Int; cached::Bool = true)
  Generic.FreeModule(R, dim; cached=cached)
end
