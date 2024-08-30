# Data type and parent object methods

parent_type(::Type{PermGroupRingElem{T}}) where {T<:RingElement} = PermGroupRing{T}

elem_type(::Type{PermGroupRing{T}}) where {T<:RingElement} = PermGroupRingElem{T}

base_ring_type(::Type{PermGroupRing{T}}) where {T<:RingElement} = parent_type(T)

base_ring(R::PermGroupRing) = R.base_ring::base_ring_type(R)

parent(f::PermGroupRingElem) = f.parent

is_domain_type(::Type{PermGroupRingElem{T}}) where {T<:RingElement} = is_domain_type(T)

is_exact_type(::Type{PermGroupRingElem{T}}) where {T<:RingElement} = is_exact_type(T)

function deepcopy_internal(f::PermGroupRingElem{T}, dict::IdDict) where {T<:RingElement}
    r = PermGroupRingElem{T}(deepcopy_internal(f.coeffs, dict))
    r.parent = f.parent # parent should not be deepcopied
    return r
end

# Basic manipulation

zero(R::PermGroupRing) = R()

one(R::PermGroupRing) = R(1)

iszero(f::PermGroupRingElem) = f == 0

isone(f::PermGroupRingElem) = f == 1

# String I/O

function show(io::IO, R::PermGroupRing)
    print(io, "Permutation group ring over ")
    show(io, base_ring(R))
 end

 function show(io::IO, f::PermGroupRingElem)
    print(io, f.coeffs)
 end

# Arithmetic functions

function -(a::PermGroupRingElem{T}) where {T<:RingElement}
    r = parent(a)()
    for (k, v) in a.coeffs
        r.coeffs[k] = -v
    end
    return r
end

function +(a::PermGroupRingElem{T}, b::PermGroupRingElem{T}) where {T<:RingElement}
    r = parent(a)()
    for (k, v) in a.coeffs
        r.coeffs[k] = v
    end
    for (k, v) in b.coeffs
        if haskey(r.coeffs, k)
            r.coeffs[k] += v
            if r.coeffs[k] == 0
                delete!(r.coeffs, k)
            end
        else
            r.coeffs[k] = v
        end
    end
    return r
end

-(a::PermGroupRingElem{T}, b::PermGroupRingElem{T}) where {T<:RingElement} = a + (-b)

function *(a::PermGroupRingElem{T}, b::PermGroupRingElem{T}) where {T<:RingElement}
    r = parent(a)()
    for (k1, v1) in a.coeffs
        for (k2, v2) in b.coeffs
            k = k1 * k2
            if haskey(r.coeffs, k)
                r.coeffs[k] += v1 * v2
            else
                r.coeffs[k] = v1 * v2
            end
        end
    end
    filter!(x -> x[2] != 0, r.coeffs)
    return r
end

function ^(a::PermGroupRingElem{T}, n::Int) where {T<:RingElement}
    if n == 0
        return one(parent(a))
    elseif n == 1
        return a
    elseif n < 0
        DomainError(n)
    end
    r = *(repeat([a], n)...)
    return r
end

# Ad hoc arithmetic functions

*(a::PermGroupRingElem{T}, n::T) where {T<:RingElement} = a * parent(a)(n)

*(n::T, a::PermGroupRingElem{T}) where {T<:RingElement} = a * n

+(a::PermGroupRingElem{T}, n::T) where {T<:RingElement} = a + parent(a)(n)

+(n::T, a::PermGroupRingElem{T}) where {T<:RingElement} = a + n

*(a::PermGroupRingElem{T}, p::Perm) where {T<:RingElement} = a * parent(a)(p)

*(p::Perm, a::PermGroupRingElem{T}) where {T<:RingElement} = parent(a)(p) * a

+(a::PermGroupRingElem{T}, p::Perm) where {T<:RingElement} = a + parent(a)(p)

+(p::Perm, a::PermGroupRingElem{T}) where {T<:RingElement} = a + p

# Comparison

function ==(a::PermGroupRingElem{T}, b::PermGroupRingElem{T}) where {T<:RingElement}
    if length(a.coeffs) != length(b.coeffs)
        return false
    end
    for (k, v) in a.coeffs
        if !haskey(b.coeffs, k) || b.coeffs[k] != v
            return false
        end
    end
    return true
end

# Ad hoc comparison

==(a::PermGroupRingElem{T}, n::T) where {T<:RingElement} = length(a.coeffs) == 1 && a.coeffs[Perm(parent(a).l)] == n

==(n::T, a::PermGroupRingElem{T}) where {T<:RingElement} = a == n

==(a::PermGroupRingElem{T}, p::Perm) where {T<:RingElement} = length(a.coeffs) == 1 && a.coeffs[p] == base_ring(parent(a))(1)

==(p::Perm, a::PermGroupRingElem{T}) where {T<:RingElement} = a == p

# TODO Some functionality are expected by ring interfaces but not necessary for ECC construction,
# including hash, exact division, random generation, promotion rules

# Constructors by overloading the call syntax for parent objects

function (R::PermGroupRing{T})() where {T<:RingElement}
    r = PermGroupRingElem{T}()
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(coeffs::Dict{<:Perm,<:Union{T,Integer,Rational,AbstractFloat}}) where {T<:RingElement}
    if valtype(coeffs) == T
        for (k,v) in coeffs
            length(k.d) == R.l || error("Invalid permutation length")
            parent(v) == R || error("Unable to coerce a group ring element")
        end
    else
        coeffs = Dict(k => base_ring(R)(v) for (k, v) in coeffs)
    end
    r = PermGroupRingElem{T}(coeffs)
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(n::Union{Integer,Rational,AbstractFloat}) where {T<:RingElement}
    coeffs = iszero(n) ? Dict{Perm,T}() : Dict(Perm(R.l) => base_ring(R)(n))
    r = PermGroupRingElem{T}(coeffs)
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(n::T) where {T<:RingElement}
    base_ring(R) == parent(n) || error("Unable to coerce a group ring element")
    r = PermGroupRingElem{T}(Dict(Perm(R.l) => n))
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(p::Perm) where {T<:RingElement}
    length(p.d) == R.l || error("Invalid permutation length")
    r = PermGroupRingElem{T}(Dict(p => one(base_ring(R))))
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(f::PermGroupRingElem{T}) where {T<:RingElement}
    parent(f) == R || error("Unable to coerce a group ring element")
    return f
end

@doc raw"""
Permutation group ring constructor.

See also: [`perm_group_ring`](@ref).
"""
function perm_group_ring(R::Ring, l::Int, cached::Bool=true)
    T = elem_type(R)
    return PermGroupRing{T}(R, l, cached)
end
