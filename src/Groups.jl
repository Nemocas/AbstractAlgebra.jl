###############################################################################
#
#   Groups.jl : AbstractAlgebra groups
#
###############################################################################

###############################################################################
#
#   Fundamental operations to be implemented for all group types
#
###############################################################################

Base.one(G::Group) = throw(NotImplementedError(:parent, g))

Base.parent(g::GroupElem) = throw(NotImplementedError(:parent, g))

Base.:*(g::T, h::T) where {T<:GroupElem} = throw(NotImplementedError(:*, g, h))

Base.inv(g::GroupElem) = throw(NotImplementedError(:inv, g))


# also should implement:
#   parent_type(::Type{MyGroupElem})
#   elem_type(::Type{MyGroup})
#
# ... and many more of course (among those below and beyond), such as
#   gens(::Group)
#   ngens(::Group)
#   order(::Group)
#   order(::GroupElem)


###############################################################################
#
#   Orders of groups and group elements
#
###############################################################################

struct InfiniteOrderError{T} <: Exception
    x::T
    InfiniteOrderError(g::Union{GroupElem, Group}) = new{typeof(g)}(g)
end

function Base.showerror(io::IO, err::InfiniteOrderError{T}) where T
    println(io, "Infinite order exception with ", err.x)
    print(io, "order will only return a value when it is finite. ")
    f = if T <: Group
        "is_finite(G)"
    elseif T <: GroupElem
        "is_finiteorder(g)"
    end
    print(io, "You should check with `$f` first.")
end

@doc raw"""
    order(::Type{T} = BigInt, G::Group) where T

Return the order of $G$ as an instance of `T`.
If $G$ is of infinite order, an `InfiniteOrderError` exception will be thrown.
Use `is_finite(G)` to avoid this kind of exception.
If the order does not fit into type `T`, an `InexactError` exception will be thrown.
"""
function order(::Type{T}, G::Group) where T
    throw(NotImplementedError(:order, G))
end

@doc raw"""
    order(::Type{T} = BigInt, g::GroupElem) where T

Return the order of $g$ as an instance of `T`.
If $g$ is of infinite order, an `InfiniteOrderError` exception will be thrown.
Use `is_finiteorder(G)` to avoid this kind of exception.
If the order does not fit into type `T`, an `InexactError` exception will be thrown.
"""
function order(::Type{T}, g::GroupElem) where T
    throw(NotImplementedError(:order, g))
end

# if no return type has been specified, default to `BigInt`
order(G::Group) = order(BigInt, G)
order(g::GroupElem) = order(BigInt, g)

"""
    is_finiteorder(g::GroupElem)

Return `true` if `g` is of finite order, possibly without computing it.
"""
function is_finiteorder(g::GroupElem)
    is_finite(parent(g)) && return true
    throw(NotImplementedError(:is_finiteorder, g))
end


###############################################################################
#
#   Properties
#
###############################################################################

@doc raw"""
    is_trivial(G::Group)

Test whether group $G$ is trivial.

The default implementation first checks whether the group has known
generators via `has_gens`, and if so, uses `gens` to test whether all
generators are the identity via `isone`. If no generators are available,
it uses `is_finite` and `order`.
"""
function is_trivial(G::Group)
    has_gens(G) && return all(isone, gens(G))
    is_finite(G) && return isone(order(G))
    return false
end

###############################################################################
#
#   Generating sets
#
###############################################################################

@doc raw"""
    has_gens(G::Group)

Test whether group $G$ has a generating set.

Many algorithms for groups make use of a (finite) generating set. But
some groups don't have one (e.g. because none has been computed yet, or
the group is not finitely generated). This test functions allows writing
algorithms that deal with this case gracefully.

Note that the result of this function when applied to a specific group
instance can change over time, as side effect of a generating set becoming
available for the group.

The default implementation returns `true`.
"""
has_gens(G::Group) = true

gens(G::Group) = throw(NotImplementedError(:gens, G))
number_of_generators(G::Group) = throw(NotImplementedError(:number_of_generators, G))

###############################################################################
#
#   Arithmetic
#
###############################################################################

Base.one(g::GroupElem) = one(parent(g))
Base.similar(g::GroupElem) = one(g)
Base.isone(g::GroupElem) = g == one(g)


"""
    conj(g::T, h::T) where {T <: GroupElem}

Return the conjugation of `g` by `h`, i.e. `inv(h)*g*h`.
"""
Base.conj(g::T, h::T) where {T<:GroupElem} = h\g*h

"""
    ^(g::T, h::T) where {T <: GroupElem}

Alias for [`conj`](@ref conj).
"""
Base.:(^)(g::T, h::T) where {T<:GroupElem} = conj(g, h)

"""
    comm(g::T, h::T, k::T...) where {T <: GroupElem}

Return the left associative iterated commutator ``[[g, h], ...]``, where
``[g, h] = g^{-1} h^{-1} g h``.
"""
function comm(g::T, h::T, k::T...) where {T<:GroupElem}
    res = comm!(similar(g), g, h)
    for l in k
        res = comm!(res, res, l)
    end
    return res
end

function Base.:(/)(g::T, h::T) where {T<:GroupElem}
    return g*inv(h)
end

function Base.:(\)(g::T, h::T) where {T<:GroupElem}
    return inv(g)*h
end


Base.literal_pow(::typeof(^), g::GroupElem, ::Val{-1}) = inv(g)

function Base.:(^)(g::GroupElem, n::Integer)
    n < 0 && return inv(g)^-n
    return Base.power_by_squaring(g, n)
end

################################################################################
# Mutable API where modifications are recommended for performance reasons
################################################################################

"""
    one!(g::GroupElem)

Return `one(g)`, possibly modifying `g`.
"""
one!(g::GroupElem) = one(parent(g))

"""
    inv!(out::T, g::T) where {GEl <: GroupElem}

Return `inv(g)`, possibly modifying `out`. Aliasing of `g` with `out` is
allowed.
"""
inv!(out::T, g::T) where {T<:GroupElem} = inv(g)

"""
    mul!(out::T, g::T, h::T) where {GEl <: GroupElem}

Return `g*h`, possibly modifying `out`. Aliasing of `g` or `h` with `out` is
allowed.
"""
mul!(out::T, g::T, h::T) where {T<:GroupElem} = g * h

"""
    div_right!(out::T, g::T, h::T) where {GEl <: GroupElem}

Return `g*inv(h)`, possibly modifying `out`. Aliasing of `g` or `h` with `out`
is allowed.
"""
function div_right!(out::T, g::T, h::T) where {T<:GroupElem}
    tmp = (out === g) ? inv(h) : inv!(out, h)
    return mul!(out, g, tmp)
end

"""
    div_left!(out::T, g::T, h::T) where {GEl <: GroupElem}

Return `inv(h)*g`, possibly modifying `out`. Aliasing of `g` or `h` with `out`
is allowed.
"""
function div_left!(out::T, g::T, h::T) where {T<:GroupElem}
    tmp = (out === g) ? inv(h) : inv!(out, h)
    return mul!(out, tmp, g)
end

"""
    conj!(out::T, g::T, h::T) where {GEl <: GroupElem}

Return `inv(h)*g*h`, possibly modifying `out`. Aliasing of `g` or `h` with
`out` is allowed.
"""
function conj!(out::T, g::T, h::T) where {T<:GroupElem}
    tmp = (out === g || out === h) ? inv(h) : inv!(out, h)
    tmp = (out === h) ? tmp * g : mul!(out, tmp, g)
    return mul!(out, tmp, h)
end

"""
    comm!(out::T, g::T, h::T) where {GEl <: GroupElem}

Return `inv(g)*inv(h)*g*h`, possibly modifying `out`. Aliasing of `g` or `h`
with `out` is allowed.
"""
function comm!(out::T, g::T, h::T) where {T<:GroupElem}
    # TODO: can we make comm! with 3 arguments without allocation??
    tmp = (out === g) ? conj(g, h) : conj!(out, g, h)
    return div_left!(out, tmp, g)
end
