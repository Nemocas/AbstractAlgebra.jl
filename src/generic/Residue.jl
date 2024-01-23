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

parent_type(::Type{EuclideanRingResidueRingElem{T}}) where T <: RingElement = EuclideanRingResidueRing{T}

elem_type(::Type{EuclideanRingResidueRing{T}}) where {T <: RingElement} = EuclideanRingResidueRingElem{T}

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{EuclideanRingResidueRingElem{T}}, ::Type{EuclideanRingResidueRingElem{T}}) where T <: RingElement = EuclideanRingResidueRingElem{T}

function promote_rule(::Type{EuclideanRingResidueRingElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? EuclideanRingResidueRingElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::EuclideanRingResidueRing{T})(b::RingElement) where {T <: RingElement}
   return a(base_ring(a)(b))
end

function (a::EuclideanRingResidueRing{T})() where {T <: RingElement}
   z = EuclideanRingResidueRingElem{T}(zero(base_ring(a)))
   z.parent = a
   return z
end

function (a::EuclideanRingResidueRing{T})(b::Integer) where {T <: RingElement}
   z = EuclideanRingResidueRingElem{T}(mod(base_ring(a)(b), modulus(a)))
   z.parent = a
   return z
end

function (a::EuclideanRingResidueRing{T})(b::T) where {T <: RingElem}
   base_ring(a) != parent(b) && error("Operation on incompatible objects")
   z = EuclideanRingResidueRingElem{T}(mod(b, modulus(a)))
   z.parent = a
   return z
end

function (a::EuclideanRingResidueRing{T})(b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   a != parent(b) && error("Operation on incompatible objects")
   return b
end

################################################################################
#
#  Map
#
################################################################################

domain(f::EuclideanRingResidueMap) = f.domain

codomain(f::EuclideanRingResidueMap) = f.codomain

function image(f::EuclideanRingResidueMap, a)
  return codomain(f)(a)
end

(f::EuclideanRingResidueMap)(a) = image(f, a)

function preimage(f::EuclideanRingResidueMap, a)
  parent(a) != codomain(f) && error("Not an element of the codomain")
  return data(a)
end
