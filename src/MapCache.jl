################################################################################
#
#  MapCache.jl : Caching of maps
#
################################################################################

@doc raw"""
    cached(M::Map; limit::Int=100, enabled::Bool=true)

Return a cached version of the map `M`, having the same supertype as `M`, that
memoises up to `limit` values of `M` in a dictionary. Setting `enabled` to
`false` creates the map with its cache switched off; it can be switched on
later with [`enable_cache!`](@ref).

# Examples

```jldoctest
julia> f = map_from_func(x -> x + 1, ZZ, ZZ);

julia> g = cached(f);

julia> f(ZZ(1)) == g(ZZ(1))
true
```
"""
function cached(M::Map; limit::Int = 100, enabled::Bool = true)
   return Generic.MapCache(M, limit, enabled)
end
