###############################################################################
#
#   residue_field.jl : generic residue fields (modulo a principal ideal)
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{ResidueFieldElem{T}}) where T <: RingElement = ResidueField{T}

elem_type(::Type{ResidueField{T}}) where {T <: RingElement} = ResidueFieldElem{T}

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{ResidueFieldElem{T}}, ::Type{ResidueFieldElem{T}}) where T <: RingElement = ResidueFieldElem{T}

function promote_rule(::Type{ResidueFieldElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? ResidueFieldElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::ResidueField{T})(b::RingElement) where {T <: RingElement}
   return a(base_ring(a)(b))
end

function (a::ResidueField{T})() where {T <: RingElement}
   z = ResidueFieldElem{T}(zero(base_ring(a)))
   z.parent = a
   return z
end

function (a::ResidueField{T})(b::Integer) where {T <: RingElement}
   z = ResidueFieldElem{T}(mod(base_ring(a)(b), modulus(a)))
   z.parent = a
   return z
end

function (a::ResidueField{T})(b::T) where {T <: RingElem}
   base_ring(a) != parent(b) && error("Operation on incompatible objects")
   z = ResidueFieldElem{T}(mod(b, modulus(a)))
   z.parent = a
   return z
end

function (a::ResidueField{T})(b::AbstractAlgebra.ResFieldElem{T}) where {T <: RingElement}
   a != parent(b) && error("Operation on incompatible objects")
   return b
end
