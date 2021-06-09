###############################################################################
#
#   ResidueField.jl : generic residue fields (modulo a principal ideal)
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{ResF{T}}) where T <: RingElement = ResField{T}

elem_type(::Type{ResField{T}}) where {T <: RingElement} = ResF{T}

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{ResF{T}}, ::Type{ResF{T}}) where T <: RingElement = ResF{T}

function promote_rule(::Type{ResF{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? ResF{T} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::ResField{T})(b::RingElement) where {T <: RingElement}
   return a(base_ring(a)(b))
end

function (a::ResField{T})() where {T <: RingElement}
   z = ResF{T}(zero(base_ring(a)))
   z.parent = a
   return z
end

function (a::ResField{T})(b::Integer) where {T <: RingElement}
   z = ResF{T}(mod(base_ring(a)(b), modulus(a)))
   z.parent = a
   return z
end

function (a::ResField{T})(b::T) where {T <: RingElem}
   base_ring(a) != parent(b) && error("Operation on incompatible objects")
   z = ResF{T}(mod(b, modulus(a)))
   z.parent = a
   return z
end

function (a::ResField{T})(b::AbstractAlgebra.ResFieldElem{T}) where {T <: RingElement}
   a != parent(b) && error("Operation on incompatible objects")
   return b
end
