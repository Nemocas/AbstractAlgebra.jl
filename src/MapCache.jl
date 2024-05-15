################################################################################
#
#  MapCache.jl : Caching of maps
#
################################################################################

function cached(M::Map; limit::Int = 100, enabled::Bool = true)
  return Generic.MapCache(M, limit, enabled)
end
