# Cached maps

All basic map (i.e. those not built up from other maps) in AbstractAlgebra can be
cached.

A cache is a dictionary that can be switched on and off at run time that keeps a cache
of previous evaluations of the map. This can be useful if the map is extremely difficult
to evaluate, e.g. a discrete logarithm map. Rather than evaluate the map afresh each
time, the map first looks up the dictionary of previous known values of the map.

To facilitate caching of maps, the Generic module provides a type `Generic.MapCache`,
which can be used to wrap any existing map object with a dictionary.

Importantly, the supertype of the resulting `MapCache` object is identical to that of
the map being cached. This means that any functions that would accept the original
map will also accept the cached version.

Note that caching of maps only works for maps that correctly abstract access to their
fields using accessor functions, as described in the map interface.

## Cached map constructors

To construct a cached map from an existing map object, we have the following function:

```julia
cached(M::Map; enabled=true, limit=100)
```

Return a cached map with the same supertype as $M$, caching up to `limit` values of the
map `M` in a dictionary, assuming that the cache is enabled.

Caches can be disabled by setting the value of the parameter `enabled` to `false`. This
allows for the user to quickly go through code and completely disable caches of maps that
were previously enabled, for testing purposes, etc.

Caches can also be turned on and off at run time (see below).

**Examples**

```julia
f = map_from_func(ZZ, ZZ, x -> x + 1)
g = cached(f)

f(ZZ(1)) == g(ZZ(1))
```

## Functionality for cached maps

The following functions are provided for cached maps.

```julia
enable_cache!(M::MapCache)
disable_cache!(M::MapCache)
```

Temporarily enable or disable the cache for the given map. The values stored in the cache
are not lost when it is disabled.

```julia
set_limit!(M::MapCache, limit::Int)
```

Set the limit on the number of values that can be cached in the dictionary, to the given
value. Setting the value to 0 will effectively disable further caching for this map.

**Examples**

```julia
f = cached(map_from_func(ZZ, ZZ, x -> x + 1))

a = f(ZZ(1))
disable_cache!(f)
b = f(ZZ(1))
enable_cache!(f)
c = f(ZZ(1))

set_limit!(f, 200)
d = f(ZZ(1))
```


