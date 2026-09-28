################################################################################
#
#  MapCache.jl : Caching of maps
#
################################################################################

domain(M::MapCache{D, C}) where {D, C} = domain(M.map)::D
codomain(M::MapCache{D, C}) where {D, C} = codomain(M.map)::C
image_fn(M::MapCache{D, C}) where {D, C} = image_fn(M.map)

@doc raw"""
    set_limit!(M::Generic.MapCache, limit::Int)

Set the number of further values that may be stored in the cache of `M` to
`limit`. Setting it to `0` prevents any further values from being cached.

# Examples

```jldoctest
julia> f = cached(map_from_func(x -> x + 1, ZZ, ZZ));

julia> set_limit!(f, 200)
200
```
"""
function set_limit!(M::MapCache, limit::Int)
   limit < 0 && error("Limit must be non-negative")
   M.limit = limit
end

@doc raw"""
    enable_cache!(M::Generic.MapCache)

Switch the cache of `M` on. Values stored in the cache while it was disabled
are kept.

# Examples

```jldoctest
julia> f = cached(map_from_func(x -> x + 1, ZZ, ZZ); enabled=false);

julia> enable_cache!(f)

julia> f(ZZ(1))
2
```
"""
function enable_cache!(M::MapCache{D, C, S, T, De, Ce}) where {D, C, S, T, De, Ce}
   if !isdefined(M, :image_cache)
      M.image_cache = Dict{De, Ce}()
   end
   M.enabled = true
   nothing
end

@doc raw"""
    disable_cache!(M::Generic.MapCache)

Switch the cache of `M` off, keeping the values it already stores. See
[`enable_cache!`](@ref).

# Examples

```jldoctest
julia> f = cached(map_from_func(x -> x + 1, ZZ, ZZ));

julia> disable_cache!(f)

julia> f(ZZ(1))
2
```
"""
function disable_cache!(M::MapCache)
   M.enabled = false
   nothing
end

(M::MapCache{D, C, S, T, De, Ce})(a::De) where {D, C, S, T, De, Ce} = image(M, a)

function image(M::MapCache{D, C, S, T, De, Ce}, a::De) where {D, C, S, T, De, Ce}
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

