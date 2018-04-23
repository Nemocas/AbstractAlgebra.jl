################################################################################
#
#  MapCache.jl : Caching of maps
#
################################################################################

export cached, disable_cache!, enable_cache!, set_limit!

get_field(M::MapCache, f) = getfield(M.map, f) # pass accessors through
set_field!(M::MapCache, f) = setfield(M.map, f) # pass accessors through

function cached(M::AbstractAlgebra.Map; limit::Int = 100, enabled::Bool = true)
   return MapCache(M, limit, enabled)
end

domain(M::MapCache{D, C}) where {D, C} = domain(M.map)::D
codomain(M::MapCache{D, C}) where {D, C} = codomain(M.map)::C

function set_limit!(M::MapCache, limit::Int)
   limit < 0 && error("Limit must be nonnegative")
   M.limit = limit
end

function enable_cache!(M::MapCache{D, C, S, T, De, Ce}) where {D, C, S, T, De, Ce}
   if !isdefined(M, :image_cache)
      M.image_cache = Dict{De, Ce}()
   end
   M.enabled = true
   nothing
end

function disable_cache!(M::MapCache)
   M.enabled = false
   nothing
end

function (M::MapCache{D, C, S, T, De, Ce})(a::De) where {D, C, S, T, De, Ce}
   if M.enabled
      if haskey(M.image_cache, a)
         return M.image_cache[a]::Ce
      else
         b = M.map(a)
         if M.limit > 0
            M.image_cache[a] = b
            M.limit -= 1
         end
         return b::Ce
      end
   else
      return M.map(a)::Ce
   end
end

show(io::IO, M::MapCache) = show(M.map)

