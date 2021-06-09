################################################################################
#
#  MapCache.jl : Caching of maps
#
################################################################################

export cached, disable_cache!, enable_cache!, set_limit!

function cached(M::Map; limit::Int = 100, enabled::Bool = true)
   return Generic.MapCache(M, limit, enabled)
end
