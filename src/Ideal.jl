###############################################################################
#
#   Ideal constructor
#
###############################################################################

# We assume that the function
#   ideal(R::T, xs::Vector{U})
# with U === elem_type(T) is implemented by anyone implementing ideals
# for AbstractAlgebra rings.
# The functions in this file extend the interface for `ideal`.

# the following helper enables things like `ideal(R, [])` or `ideal(R, [1])`
# the type check ensures we don't run into an infinite recursion
function ideal(R::NCRing, xs::AbstractVector{T}; kw...) where T<:NCRingElement
  xs isa Vector{elem_type(R)} && error("ideals unsupported for ring $R")
  return ideal(R, elem_type(R)[R(x) for x in xs]; kw...)
end

function ideal(R::NCRing, x, y...; kw...)
  return ideal(R, elem_type(R)[R(z) for z in [x, y...]]; kw...)
end

function *(R::Ring, x::RingElement)
  return ideal(R, x)
end

function *(x::RingElement, R::Ring)
  return ideal(R, x)
end

function ideal(x::NCRingElement; kw...)
  return ideal(parent(x), x; kw...)
end

function ideal(xs::AbstractVector{T}; kw...) where T<:NCRingElement
  !is_empty(xs) || throw(ArgumentError("Empty collection, cannot determine parent ring. Try ideal(ring, xs) instead of ideal(xs)"))
  return ideal(parent(xs[1]), xs; kw...)
end

iszero(I::Ideal) = all(iszero, gens(I))

base_ring_type(::Type{<:IdealSet{T}}) where T <: RingElement = parent_type(T)
