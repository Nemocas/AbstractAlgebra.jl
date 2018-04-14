################################################################################
#
#  MapCache.jl : Caching of maps
#
################################################################################

export cached, set_limit!

function cached(M::AbstractAlgebra.Map; limit::Int = 100, enable::Bool = true)
   return MapCache(M, limit, enable)
end

domain(M::MapCache) = domain(M.map)
codomain(M::MapCache) = codomain(M.map)

function set_limit!(M::MapCache, limit::Int)
   limit < 0 && error("Limit must be nonnegative")
   M.limit = limit
end

function (M::MapCache{D, C, T, De, Ce})(a::De) where {D, C, T, De, Ce}
   if isdefined(M, :image_cache)
      if haskey(M.image_cache, a)
         return M.image_cache[a]::Ce
      else
         b = M.map(a)
         if M.limit > 0
            M.image_cache[a] = b
            M.limit -= 1
         end
         return b
      end
   else
      return M.map(a)
   end
end

show(io::IO, M::MapCache) = show(M.map)

