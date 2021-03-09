# Weak value dictionaries
#
# This file is a modification of base/weakkeydict.jl from the Julia project.

import Base: isempty, setindex!, getkey, length, iterate, empty, delete!, pop!, get, sizehint!, copy
"""
    WeakValueDict([itr])

`WeakValueDict()` constructs a hash table where the values are weak
references to objects which may be garbage collected even when
used as values in a hash table.
"""
mutable struct WeakValueDict{K,V} <: AbstractDict{K,V}
    ht::Dict{K,WeakRef}
    lock::ReentrantLock
    finalizer::Function
    dirty::Bool

    # Constructors mirror Dict's
    function WeakValueDict{K,V}() where V where K
        t = new(Dict{K,Any}(), ReentrantLock(), identity, false)
        t.finalizer = v -> t.dirty = true
        return t
    end
end
function WeakValueDict{K,V}(kv) where V where K
    h = WeakValueDict{K,V}()
    for (k,v) in kv

        h[k] = v
    end
    return h
end
WeakValueDict{K,V}(p::Pair) where V where K = setindex!(WeakValueDict{K,V}(), p.second, p.first)

function WeakValueDict{K,V}(ps::Pair...) where V where K
    h = WeakValueDict{K,V}()
    sizehint!(h, length(ps))
    for p in ps
        h[p.first] = p.second
    end
    return h
end
WeakValueDict() = WeakValueDict{Any,Any}()

WeakValueDict(kv::Tuple{}) = WeakValueDict()
Base.copy(d::WeakValueDict) = WeakValueDict(d)

WeakValueDict(ps::Pair{K,V}...)           where {K,V} = WeakValueDict{K,V}(ps)
WeakValueDict(ps::Pair{K}...)             where {K}   = WeakValueDict{K,Any}(ps)
WeakValueDict(ps::(Pair{K,V} where K)...) where {V}   = WeakValueDict{Any,V}(ps)
WeakValueDict(ps::Pair...)                            = WeakValueDict{Any,Any}(ps)

function WeakValueDict(kv)
    try
        Base.dict_with_eltype((K, V) -> WeakValueDict{K, V}, kv, eltype(kv))
    catch
        if !Base.isiterable(typeof(kv)) || !all(x->isa(x,Union{Tuple,Pair}),kv)
            throw(ArgumentError("WeakValueDict(kv): kv needs to be an iterator of tuples or pairs"))
        else
            rethrow()
        end
    end
end

function _cleanup_locked(h::WeakValueDict)
    if h.dirty
        h.dirty = false
        idx = Base.skip_deleted_floor!(h.ht)
        while idx != 0
            if h.ht.vals[idx].value === nothing
                Base._delete!(h.ht, idx)
            end
            idx = Base.skip_deleted(h.ht, idx + 1)
        end
    end
    return h
end

Base.sizehint!(d::WeakValueDict, newsz) = sizehint!(d.ht, newsz)
empty(d::WeakValueDict, ::Type{K}, ::Type{V}) where {K, V} = WeakValueDict{K, V}()

IteratorSize(::Type{<:WeakValueDict}) = SizeUnknown()

islocked(wkh::WeakValueDict) = Base.islocked(wkh.lock)
lock(f, wkh::WeakValueDict) = Base.lock(f, wkh.lock)
trylock(f, wkh::WeakValueDict) = Base.trylock(f, wkh.lock)

function Base.setindex!(wkh::WeakValueDict{K}, v, key) where K
    !isa(key, K) && throw(ArgumentError("$(Base.limitrepr(key)) is not a valid key for type $K"))
    # 'nothing' is not valid both because 'finalizer' will reject it,
    # and because we therefore use it as a sentinel value
    v === nothing && throw(ArgumentError("`nothing` is not a valid WeakValueDict key"))
    lock(wkh) do
        _cleanup_locked(wkh)
        k = getkey(wkh.ht, key, nothing)
        # The object gets the finalizer no matter if the key is already there or not.
        finalizer(wkh.finalizer, v)
        if k === nothing
            k = key
        end
        wkh.ht[k] = WeakRef(v)
    end
    return wkh
end
function get!(wkh::WeakValueDict{K}, key, default) where {K}
    v = lock(wkh) do
        if key !== nothing && haskey(wkh.ht, key)
            wkh.ht[key].value
        else
            wkh[key] = WeakRef(default)
            default
        end
    end
    return v
end
function Base.get!(default::Base.Callable, wkh::WeakValueDict{K}, key) where {K}
    v = lock(wkh) do
        _cleanup_locked(wkh)
        if key !== nothing && haskey(wkh.ht, key)
            wkh.ht[key].value
        else
            wkh[key] = default()
        end
    end
    return v
end

function Base.getkey(wkh::WeakValueDict{K}, kk, default) where K
    k = lock(wkh) do
        _cleanup_locked(wkh)
        getkey(wkh.ht, kk, nothing)
    end
    return k === nothing ? default : k::K
end

map!(f, iter::Base.ValueIterator{<:WeakValueDict})= map!(f, values(iter.dict.ht))

function Base.get(wkh::WeakValueDict{K, V}, key, default) where {K, V}
    key === nothing && throw(KeyError(nothing))
    lock(wkh) do
        _cleanup_locked(wkh)
        x = get(wkh.ht, key, default)
        if x === default
            return x
        else
            return x.value::V
        end
    end
end
function Base.get(default::Base.Callable, wkh::WeakValueDict{K}, key) where {K}
    key === nothing && throw(KeyError(nothing))
    lock(wkh) do
        return get(default, wkh.ht, key)
    end
end
function Base.pop!(wkh::WeakValueDict{K}, key) where {K}
    key === nothing && throw(KeyError(nothing))
    lock(wkh) do
        return pop!(wkh.ht, key).value
    end
end
function Base.pop!(wkh::WeakValueDict{K}, key, default) where {K}
    key === nothing && return default
    lock(wkh) do
        return pop!(wkh.ht, key, default)
    end
end
function Base.delete!(wkh::WeakValueDict, key)
    key === nothing && return wkh
    lock(wkh) do
        delete!(wkh.ht, key)
    end
    return wkh
end
function Base.empty!(wkh::WeakValueDict)
    lock(wkh) do
        empty!(wkh.ht)
    end
    return wkh
end
function Base.haskey(wkh::WeakValueDict{K}, key) where {K}
    key === nothing && return false
    lock(wkh) do
        return haskey(wkh.ht, key)
    end
end
function Base.getindex(wkh::WeakValueDict{K}, key) where {K}
    key === nothing && throw(KeyError(nothing))
    lock(wkh) do
        return getindex(wkh.ht, key).value
    end
end

Base.isempty(wkh::WeakValueDict) = length(wkh) == 0
function Base.length(t::WeakValueDict)
    lock(t) do
        _cleanup_locked(t)
        return length(t.ht)
    end
end

function Base.iterate(t::WeakValueDict{K,V}, state...) where {K, V}
    return lock(t) do
        while true
            y = iterate(t.ht, state...)
            y === nothing && return nothing
            wkv, state = y
            k = wkv[1]
            v = wkv[2].value
            if VERSION >= v"1.6"
              GC.safepoint() # ensure `v` is now gc-rooted
            end
            k === nothing && continue # indicates `k` is scheduled for deletion
            kv = Pair{K,V}(k::K, v)
            return (kv, state)
        end
    end
end

filter!(f, d::WeakValueDict) = Base.filter_in_one_pass!(f, d)
