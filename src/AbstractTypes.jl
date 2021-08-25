abstract type Ring end

abstract type Field <: Ring end

abstract type RingElem end

abstract type FieldElem <: RingElem end

abstract type PolyRing{T} <: Ring end

abstract type FracField{T} <: Field end

abstract type PolyElem{T} <: RingElem end

abstract type FracElem{T} <: FieldElem end

promote_rule(T, U) = Union{}

promote_rule(a::Type{S}, b::Type{T}) where {S <: Real, T <: Real} = Base.promote_rule(a, b)
