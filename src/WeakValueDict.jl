# Weak value containers

import Base: isempty, setindex!, getkey, length, iterate, empty, delete!, pop!,
       get, sizehint!, copy, empty!, haskey, getindex

################################################################################
#
# WeakValueCache
# This section is a modification of base/dict.jl from the Julia project to
# support disappearing entries. It might work if isfinalized can be implemented.
#
################################################################################

#=
The weak value containers cannot be useful without special modifications to the
value types. Julia's original sin is that the WeakRef type is broken by design
(https://github.com/JuliaLang/julia/issues/35169). A typical mutable object
may have a life like

   A. born
   B. in use and modified by program
   C. declared dead
   D. some finalizers may be called on the object (or one of its members)
   E. actually dead: object's memory is collected

The problem is that .value of a WeakRef can be an object between D and E. Since
the finalizers generally put the object into an invalid state, this is useless
(the stated rationale for this terrible behaviour is that finalizers are capable
of resurrecting objects from the dead, so state D -> B is possible).

To have a somewhat working WeakValue dictionary or cache, it is necessary to
know if the object has transitioned past D. For this, we make the assumption
the finalizers attached do not resurrect objects from the dead, and the
assumption that the validity of the objects can be determined by a method such
as `isfinalized`. Examples:

   1. NmodRing: here we have

         mutable struct NmodRing <: Ring
            n::UInt
            ninv::UInt
         end

      and this is fine

         isfinalized(x::NmodRing) = false

   2. FmpzModRing: already things get bad because here we have

         mutable struct FmpzModRing <: Ring
            n::fmpz
            ninv::fmpz_mod_ctx_struct
            isfinalized::Bool          # possible addition to the struct
         end

      and the fictitious

         isfinalized(x::FmpzModRing) = x.isfinalized

      When this is declared dead, the object itself and the two members can be
      finalized in any order. We could try not attaching a finalizer to the
      .ninv, and have the finalizer for the FmpzModRing clean up the
      .ninv member and set the .isfinalized member to true, but there is no way
      to know if the finalizer on the .n has run or not, and thus the validity
      of the FmpzModRing object cannot be determined without special treatment
      of the .n member

   3. For the arbitrarily nested rings in AA, one would have to essentially
      reimplement the GC, which is not going to happen.
=#

# is x an invalid object? i.e. have finalizers on x or any part of x run?
function isfinalized(x)
   # v0.22 tests might pass with isfinalized(x) = false and
   # CacheDictType = WeakValueCache, but this only by luck.
   error("not implemented")
   return false
end

# for the .value member of a WeakRef
function isdead(x)
   return x == nothing || isfinalized(x)
end


"""
`WeakValueCache()` constructs a hash table whose values may be collected at
anytime and behave as if they have disappeared from the table once collected.
Putting something into a weak value cache `d[k] = v` and not modifying the
entry only necessarily means that `d[k] == v` at a later point in time if `v`
is loosely _in use_ in between. `WeakValueCache` does not work unless the
predicate `isfinalized` is properly implemented for the type of the values.
"""
mutable struct WeakValueCache{K, V}
   slots::Vector{UInt8}
   keys::Vector{K}
   vals::Vector{WeakRef}
   nmissing::Int
   nfilled::Int
   maxprobe::Int
   age::UInt64

   function WeakValueCache{K, V}() where V where K
      n = 16
      new(zeros(UInt8,n), Vector{K}(undef, n), Vector{WeakRef}(undef, n), 0, 0, 0, 0)
   end
end

function Base.show(io::IO, h::WeakValueCache)
   print(io, "\n")
   for i in 1:max(length(h.slots), length(h.keys), length(h.vals))
      show(io, i); print(io,": ");
      isassigned(h.slots, i) ? show(io,h.slots[i]) : print(io, Base.undef_ref_str); print(io, ", ")
      isassigned(h.keys, i) ? show(io, h.keys[i]) : print(io, Base.undef_ref_str); print(io, ", ")
      if isassigned(h.vals, i)
         show(io, h.vals[i])
         x = h.vals[i].value
         if !isdead(x) && !isimmutable(x)
            print(io, " @ ")
            show(io, pointer_from_objref(x))
         end
      else
         print(io, Base.undef_ref_str)
      end
      print(io, ", ")
      print(io, "\n")
   end
   print(io," nfilled: "); show(io, h.nfilled); print(io, "\n")
   print(io,"nmissing: "); show(io, h.nmissing); print(io, "\n")
   print(io,"maxprobe: "); show(io, h.maxprobe); print(io, "\n")
   print(io,"age: "); show(io, h.age); print(io, "\n")
end

_tablesz(x::Integer) = x < 16 ? 16 : one(x)<<((sizeof(x)<<3)-leading_zeros(x-1))

hashindex(key, sz) = (((hash(key)::UInt % Int) & (sz-1)) + 1)::Int

isslotempty(h::WeakValueCache, i::Int) = h.slots[i] == 0x0
isslotfilled(h::WeakValueCache, i::Int) = h.slots[i] == 0x1
isslotmissing(h::WeakValueCache, i::Int) = h.slots[i] == 0x2

function _isconsistent(h::WeakValueCache)
   nfilled = nmissing = 0
   for i in 1:length(h.slots)
      nfilled += isslotfilled(h, i)
      nmissing += isslotmissing(h, i)
   end
   return nfilled == h.nfilled && nmissing == h.nmissing
end

# missing or filled or empty -> filled
function _setindex!(h::WeakValueCache, v, key, index)
   @assert _isconsistent(h)
   h.keys[index] = key
   h.vals[index] = v
   if isslotfilled(h, index)
      return h
   elseif isslotmissing(h, index)
      h.nmissing -= 1
   end
   h.nfilled += 1
   h.slots[index] = 0x1
   h.age += 1
   sz = length(h.keys)
   if h.nmissing*4 > sz*3 || h.nfilled*3 > sz*2
      rehash!(h, h.nfilled > 64000 ? h.nfilled*2 : h.nfilled*4)
   end
   return h
end

# filled -> missing
function _deleteindex!(h::WeakValueCache, index) where {K,V}
   @assert _isconsistent(h)
   @assert isslotfilled(h, index)
   h.nfilled -= 1
   h.nmissing += 1
   h.slots[index] = 0x2
   h.age += 1
   Base._unsetindex!(h.keys, index)
   Base._unsetindex!(h.vals, index)
   return h
end

function rehash!(h::WeakValueCache{K,V}, newsz = length(h.keys)) where V where K
   olds = h.slots
   oldk = h.keys
   oldv = h.vals
   sz = length(olds)
   newsz = _tablesz(newsz)
   h.age += 1
   if h.nfilled == 0
      resize!(h.slots, newsz)
      fill!(h.slots, 0x0)
      resize!(h.keys, newsz)
      resize!(h.vals, newsz)
      h.nmissing = 0
      return h
   end
   slots = zeros(UInt8, newsz)
   keys = Vector{K}(undef, newsz)
   vals = Vector{WeakRef}(undef, newsz)
   count = 0
   maxprobe = 0
   for i = 1:sz
      if olds[i] == 0x1
         k = oldk[i]
         v = oldv[i]
         if !isdead(v.value)
            index0 = index = hashindex(k, newsz)
            while slots[index] != 0x0
               index = (index & (newsz-1)) + 1
            end
            probe = (index - index0) & (newsz-1)
            maxprobe = max(maxprobe, probe)
            slots[index] = 0x1
            keys[index] = k
            vals[index] = v
            count += 1
         end
      end
   end
   h.slots = slots
   h.keys = keys
   h.vals = vals
   h.nfilled = count
   h.nmissing = 0
   h.maxprobe = maxprobe
   return h
end

function empty!(h::WeakValueCache{K,V}) where V where K
   fill!(h.slots, 0x0)
   sz = length(h.slots)
   empty!(h.keys)
   empty!(h.vals)
   resize!(h.keys, sz)
   resize!(h.vals, sz)
   h.nmissing = 0
   h.nfilled = 0
   h.age += 1
   return h
end

# get the index where a key is stored, or -1 if not present
function ht_keyindex(h::WeakValueCache{K,V}, key) where V where K
   sz = length(h.keys)
   iter = 0
   index = hashindex(key, sz)
   while iter <= h.maxprobe
      if isslotempty(h, index)
         break
      elseif isslotfilled(h, index)
         x = h.vals[index].value
         if isdead(x)
            # filled -> missing
            _deleteindex!(h, index)
         elseif key === h.keys[index] || isequal(key, h.keys[index])
            return index
         end
      end
      index = (index & (sz-1)) + 1
      iter += 1
   end
   return -1
end

# get the index where a key is stored, or -pos if not present
# and the key would be inserted at pos
# This version is for use by setindex! and get!
function ht_keyindex2!(h::WeakValueCache{K,V}, key) where V where K
   sz = length(h.keys)
   iter = 0
   index = hashindex(key, sz)
   avail = 0
   while iter <= h.maxprobe
      if isslotempty(h, index)
         return avail < 0 ? avail : -index
      elseif isslotmissing(h, index)
         if avail == 0
            # found an available slot, but need to keep scanning
            # in case "key" already exists in a later collided slot.
            avail = -index
         end
      else
         x = h.vals[index].value
         if isdead(x)
            # filled -> missing
            _deleteindex!(h, index)
            if avail == 0
               # ditto
               avail = -index
            end
         elseif key === h.keys[index] || isequal(key, h.keys[index])
            return index
         end
      end
      index = (index & (sz-1)) + 1
      iter += 1
   end
   avail < 0 && return avail
   # Check if key is not present, may need to keep searching to find slot
   while iter <= max(16, sz>>6)
      if isslotfilled(h,index)
         x = h.vals[index].value
         if isdead(x)
            h.maxprobe = iter
            return -index
         end
      else
         h.maxprobe = iter
         return -index
      end
      index = (index & (sz-1)) + 1
      iter += 1
   end
   rehash!(h, h.nfilled > 64000 ? sz*2 : sz*4)
   return ht_keyindex2!(h, key)
end

### delete! ####

function delete!(h::WeakValueCache, key)
   index = ht_keyindex(h, key)
   if index > 0
      _deleteindex!(h, index)
   end
   return h
end

### getindex / haskey ####

function haskey(h::WeakValueCache, key)
   index = ht_keyindex(h, key)
   return index > 0 && !isdead(h.vals[index].value)
end

function getindex(h::WeakValueCache{K, V}, key) where {K, V}
   index = ht_keyindex(h, key)
   if index > 0
      x = h.vals[index].value
      if !isdead(x)
         return x::V
      end
   end
   throw(KeyError(key))
end

### get! ####

function Base.get!(default::Base.Callable, h::WeakValueCache{K, V}, key0) where {K, V}
   key = convert(K, key0)
   isequal(key, key0) || throw(ArgumentError("$key0 is not a valid key for type $K"))
   return get!(default, h, key)
end

function Base.get!(default::Base.Callable, h::WeakValueCache{K, V}, key::K) where {K, V}
   index = ht_keyindex2!(h, key)
   if index > 0
      x = h.vals[index].value
      if !isdead(x)
         return x::V
      end
   else
      index = -index
   end
   age0 = h.age
   x = convert(V, default())
   if h.age == age0
      _setindex!(h, WeakRef(x), key, index)
   else
      # the call to default could have changed the internals of h
      setindex!(h, x, key)
   end
   return x
end

### setindex! ####

function setindex!(h::WeakValueCache{K, V}, v0, key0) where {K, V}
   key = convert(K, key0)
   isequal(key, key0) || throw(ArgumentError("$key0 is not a valid key for type $K"))
   setindex!(h, v0, key)
end

function setindex!(h::WeakValueCache{K, V}, v0, key::K) where V where K
   x = convert(V, v0)
   index = ht_keyindex2!(h, key)
   _setindex!(h, WeakRef(x), key, abs(index))
   return h
end


################################################################################
#
# WeakValueDict
# This section is a modification of base/weakkeydict.jl from the Julia project,
#  and cannot be made to work without substantial changes.
#
################################################################################

#=
This WeakValueDict is broken because it does not attempt to recognize invalid
objects. The locks here are used to prevent finalizers from running while the
internal ht is being modified/inspected, but this neither useful nor necessary.
=#

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
