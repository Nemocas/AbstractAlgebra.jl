struct Integers{T <: Integer} <: Ring
end

struct Rationals{T <: Integer} <: Field
end

const RingElement   = Union{RingElem,   Integer, Rational, AbstractFloat}

const FieldElement = Union{FieldElem, Rational, AbstractFloat}
