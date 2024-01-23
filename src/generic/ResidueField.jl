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

parent_type(::Type{EuclideanRingResidueFieldElem{T}}) where T <: RingElement = EuclideanRingResidueField{T}

elem_type(::Type{EuclideanRingResidueField{T}}) where {T <: RingElement} = EuclideanRingResidueFieldElem{T}

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{EuclideanRingResidueFieldElem{T}}, ::Type{EuclideanRingResidueFieldElem{T}}) where T <: RingElement = EuclideanRingResidueFieldElem{T}

function promote_rule(::Type{EuclideanRingResidueFieldElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? EuclideanRingResidueFieldElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::EuclideanRingResidueField{T})(b::RingElement) where {T <: RingElement}
   return a(base_ring(a)(b))
end

function (a::EuclideanRingResidueField{T})() where {T <: RingElement}
   z = EuclideanRingResidueFieldElem{T}(zero(base_ring(a)))
   z.parent = a
   return z
end

function (a::EuclideanRingResidueField{T})(b::Integer) where {T <: RingElement}
   z = EuclideanRingResidueFieldElem{T}(mod(base_ring(a)(b), modulus(a)))
   z.parent = a
   return z
end

function (a::EuclideanRingResidueField{T})(b::T) where {T <: RingElem}
   base_ring(a) != parent(b) && error("Operation on incompatible objects")
   z = EuclideanRingResidueFieldElem{T}(mod(b, modulus(a)))
   z.parent = a
   return z
end

function (a::EuclideanRingResidueField{T})(b::AbstractAlgebra.ResFieldElem{T}) where {T <: RingElement}
   a != parent(b) && error("Operation on incompatible objects")
   return b
end

################################################################################
#
#  Random random functionality
#
################################################################################

function RandomExtensions.make(S::EuclideanRingResidueField{Generic.Poly{Rational{BigInt}}}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1])
   else
      n = degree(S.modulus)
      Make(S, make(base_ring(S), n - 1:n - 1, vs...))
   end
end
