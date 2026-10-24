################################################################################
#
#  Caches
#
################################################################################

const _cache_clearers = Any[]

@doc raw"""
    AbstractAlgebra.register_cache_clearer(clear)

Register `clear` to be called without arguments by
[`AbstractAlgebra.clear_cache`](@ref). Returns `clear`.

Other packages must register in their `__init__`: changes to AbstractAlgebra's
registry made while another package precompiles are not saved.
"""
function register_cache_clearer(clear)
  push!(_cache_clearers, clear)
  return clear
end

@doc raw"""
    AbstractAlgebra.register_cache(cache::Union{AbstractDict, AbstractArray})

Register `cache` to be emptied with `empty!` by
[`AbstractAlgebra.clear_cache`](@ref). Returns `cache`, so a cache can be
registered where it is defined:

    const _my_cache = AbstractAlgebra.register_cache(Dict{Int, BigInt}())

The caveat on `__init__` of [`AbstractAlgebra.register_cache_clearer`](@ref)
applies.
"""
function register_cache(cache::Union{AbstractDict, AbstractArray})
  register_cache_clearer(() -> empty!(cache))
  return cache
end

@doc raw"""
    AbstractAlgebra.clear_cache()

Empty all caches registered with [`AbstractAlgebra.register_cache`](@ref) or
[`AbstractAlgebra.register_cache_clearer`](@ref).
"""
function clear_cache()
  for clear in _cache_clearers
    clear()
  end
  return nothing
end
