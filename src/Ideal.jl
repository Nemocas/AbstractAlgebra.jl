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
function ideal(R::Ring, xs::AbstractVector{T}) where T<:RingElement
  xs isa Vector{elem_type(R)} && error("ideals unsupported for ring $R")
  return ideal(R, elem_type(R)[R(x) for x in xs])
end

function ideal(R::Ring, x, y...)
  return ideal(R, elem_type(R)[R(z) for z in [x, y...]])
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
