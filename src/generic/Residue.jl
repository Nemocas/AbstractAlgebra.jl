###############################################################################
#
#   Residue.jl : generic residue rings (modulo a principal ideal)
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{ResidueRingElem{T}}) where T <: RingElement = ResidueRing{T}

elem_type(::Type{ResidueRing{T}}) where {T <: RingElement} = ResidueRingElem{T}

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{ResidueRingElem{T}}, ::Type{ResidueRingElem{T}}) where T <: RingElement = ResidueRingElem{T}

function promote_rule(::Type{ResidueRingElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? ResidueRingElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::ResidueRing{T})(b::RingElement) where {T <: RingElement}
   return a(base_ring(a)(b))
end

function (a::ResidueRing{T})() where {T <: RingElement}
   z = ResidueRingElem{T}(zero(base_ring(a)))
   z.parent = a
   return z
end

function (a::ResidueRing{T})(b::Integer) where {T <: RingElement}
   z = ResidueRingElem{T}(mod(base_ring(a)(b), modulus(a)))
   z.parent = a
   return z
end

function (a::ResidueRing{T})(b::T) where {T <: RingElem}
   base_ring(a) != parent(b) && error("Operation on incompatible objects")
   z = ResidueRingElem{T}(mod(b, modulus(a)))
   z.parent = a
   return z
end

function (a::ResidueRing{T})(b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   a != parent(b) && error("Operation on incompatible objects")
   return b
end

