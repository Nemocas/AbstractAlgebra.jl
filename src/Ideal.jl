###############################################################################
#
#   Ideal.jl : Generic functionality for two-sided ideals
#
#
#  A the very least, an implementation of an Ideal subtype should provide the
#  following methods:
#  - ideal_type(::Type{MyRingType}) = MyIdealType
#  - ideal(R::MyRingType, xs::Vector{MyRingElemType})::MyIdealType
#  - base_ring(I::MyIdealType)
#  - gen(I::MyIdealType, k::Int)
#  - gens(I::MyIdealType)
#  - ngens(I::MyIdealType)
#  - Base.in(v::MyRingElemType, I::MyIdealType)
#  - ...
#
# Many other functions are then automatically derived from these.
###############################################################################

###############################################################################
#
#   Type and parent functions
#
###############################################################################

# fundamental interface, to be documented
ideal_type(x) = ideal_type(typeof(x))
ideal_type(T::DataType) = throw(MethodError(ideal_type, (T,)))

#
parent(I::Ideal) = DefaultIdealSet(base_ring(I))

parent_type(::Type{<:Ideal{T}}) where {T} = DefaultIdealSet{T}

#
base_ring(S::DefaultIdealSet) = S.base_ring::base_ring_type(S)

base_ring_type(::Type{<:IdealSet{T}}) where {T} = parent_type(T)

elem_type(::Type{<:IdealSet{T}}) where {T} = ideal_type(parent_type(T))

###############################################################################
#
#   Ideal constructors
#
###############################################################################

# All constructors ultimately delegate to a method
#   ideal(R::T, xs::Vector{U}) where T <: NCRing
# and U === elem_type(T)
ideal(R::T, xs::Vector{S}) where {T <: NCRing, S <: NCRingElement} = ideal_type(T)(R, xs)

# the following helper enables things like `ideal(R, [])` or `ideal(R, [1])`
# the type check ensures we don't run into an infinite recursion
function ideal(R::NCRing, xs::AbstractVector{T}; kw...) where T<:NCRingElement
  xs isa Vector{elem_type(R)} && error("ideals unsupported for ring $R")
  return ideal(R, elem_type(R)[R(x) for x in xs]; kw...)
end

function ideal(R::NCRing, x, y...; kw...)
  return ideal(R, elem_type(R)[R(z) for z in [x, y...]]; kw...)
end

function ideal(x::NCRingElement; kw...)
  return ideal(parent(x), x; kw...)
end

function ideal(xs::AbstractVector{T}; kw...) where T<:NCRingElement
  @req !is_empty(xs) "Empty collection, cannot determine parent ring. Try ideal(ring, xs) instead of ideal(xs)"
  return ideal(parent(xs[1]), xs; kw...)
end

function Base.similar(I::T, xs::Vector) where {T <: Ideal}
  R = base_ring(I)
  @assert T === ideal_type(R)
  return ideal(R, xs)
end

###############################################################################
#
#   Basic predicates
#
###############################################################################

iszero(I::Ideal) = all(iszero, gens(I))

@doc raw"""
    Base.issubset(I::T, J::T) where {T <: Ideal}

Return `true` if the ideal `I` is a subset of the ideal `J`.
"""
function Base.issubset(I::T, J::T) where {T <: Ideal}
  I === J && return true
  check_base_ring(I, J)
  return all(in(J), gens(I))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function Base.:(==)(I::T, J::T) where {T <: Ideal}
  return is_subset(I, J) && is_subset(J, I)
end

function Base.:hash(I::T, h::UInt) where {T <: Ideal}
  h = hash(base_ring(I), h)
  return h
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function Base.:+(I::T, J::T) where {T <: Ideal}
  check_base_ring(I, J)
  return similar(I, vcat(gens(I), gens(J)))
end

function Base.:*(I::T, J::T) where {T <: Ideal}
  check_base_ring(I, J)
  return similar(I, [x*y for x in gens(I) for y in gens(J)])
end

###############################################################################
#
#   Ad hoc binary operations
#
###############################################################################

function *(R::Ring, x::RingElement)
  return ideal(R, x)
end

function *(x::RingElement, R::Ring)
  return ideal(R, x)
end

function *(I::Ideal{T}, p::T) where T <: RingElement
  iszero(p) && return similar(I, T[])
  return similar(I, [v*p for v in gens(I)])
end

function *(p::T, I::Ideal{T}) where T <: RingElement
  return I*p
end

function *(I::Ideal{T}, p::S) where {S <: RingElement, T <: RingElement}
  iszero(p*one(R)) && return similar(I, T[])
  return similar(I, [v*p for v in gens(I)])
end

function *(p::S, I::Ideal{T}) where {S <: RingElement, T <: RingElement}
  return I*p
end
