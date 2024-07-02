###############################################################################
#
#   Ideal constructor
#
###############################################################################

# We assume that the functions
#   ideal(R::Ring, x::RingElement)
#   ideal(R::Ring, xs::AbstractVector)
# are implemented by anyone implementing ideals for AbstractAlgebra rings.
# The functions in this file extend the interface for `ideal`.

function *(R::Ring, x::Any)
  return ideal(R, x)
end

function *(x::Any, R::Ring)
  return ideal(R, x)
end

function ideal(R::Ring, x::Any)
  return ideal(R, R(x))
end

function ideal(x::RingElement)
  return ideal(parent(x), x)
end

function ideal(xs::AbstractVector{T}) where T<:RingElement
  !is_empty(xs) || throw(ArgumentError("Empty collection, cannot determine parent ring, try ideal(ring, xs) instead of ideal(xs)"))
  return ideal(parent(xs[1]), xs)
end

iszero(I::Ideal) = all(iszero, gens(I))

base_ring_type(::Type{<:IdealSet{T}}) where T <: RingElement = parent_type(T)
