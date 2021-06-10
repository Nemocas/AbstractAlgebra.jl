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

parent_type(::Type{Res{T}}) where T <: RingElement = ResRing{T}

elem_type(::Type{ResRing{T}}) where {T <: RingElement} = Res{T}

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{Res{T}}, ::Type{Res{T}}) where T <: RingElement = Res{T}

function promote_rule(::Type{Res{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? Res{T} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::ResRing{T})(b::RingElement) where {T <: RingElement}
   return a(base_ring(a)(b))
end

function (a::ResRing{T})() where {T <: RingElement}
   z = Res{T}(zero(base_ring(a)))
   z.parent = a
   return z
end

function (a::ResRing{T})(b::Integer) where {T <: RingElement}
   z = Res{T}(mod(base_ring(a)(b), modulus(a)))
   z.parent = a
   return z
end

function (a::ResRing{T})(b::T) where {T <: RingElem}
   base_ring(a) != parent(b) && error("Operation on incompatible objects")
   z = Res{T}(mod(b, modulus(a)))
   z.parent = a
   return z
end

function (a::ResRing{T})(b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   a != parent(b) && error("Operation on incompatible objects")
   return b
end

