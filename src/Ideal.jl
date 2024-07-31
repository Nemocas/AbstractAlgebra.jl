###############################################################################
#
#   Ideal constructor
#
###############################################################################

# We assume that the function
#   ideal(R::Ring, xs::AbstractVector{T}) where T<:RingElement
# is implemented by anyone implementing ideals for AbstractAlgebra rings.
# The functions in this file extend the interface for `ideal`.

function ideal(R::Ring, x::RingElement)
  return ideal(R, elem_type(R)[R(x)])
end

function *(R::Ring, x::RingElement)
  return ideal(R, x)
end

function *(x::RingElement, R::Ring)
  return ideal(R, x)
end

function ideal(x::RingElement)
  return ideal(parent(x), x)
end

function ideal(xs::AbstractVector{T}) where T<:RingElement
  !is_empty(xs) || throw(ArgumentError("Empty collection, cannot determine parent ring. Try ideal(ring, xs) instead of ideal(xs)"))
  return ideal(parent(xs[1]), xs)
end

iszero(I::Ideal) = all(iszero, gens(I))

base_ring_type(::Type{<:IdealSet{T}}) where T <: RingElement = parent_type(T)
