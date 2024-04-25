###############################################################################
#
#   Ideal constructor
#
###############################################################################

function ideal
end

function *(x::RingElement, R::Ring)
  return ideal(R, x)
end

function *(R::Ring, x::RingElement)
  return ideal(R, x)
end

iszero(I::Ideal) = all(iszero, gens(I))

base_ring_type(::Type{<:IdealSet{T}}) where T <: RingElement = parent_type(T)
