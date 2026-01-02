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
#   Parent object call overload
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

function (a::EuclideanRingResidueRing)(b::T) where {T <: Union{Rational, FracElem}}
   return inv(a(denominator(b))) * a(numerator(b))
end

function (a::EuclideanRingResidueRing{T})(b::T) where {T <: RingElem}
   base_ring(a) !== parent(b) && error("Operation on incompatible objects")
   z = EuclideanRingResidueRingElem{T}(mod(b, modulus(a)))
   z.parent = a
   return z
end

function (a::EuclideanRingResidueRing{T})(b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   a !== parent(b) && error("Operation on incompatible objects")
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
  parent(a) !== domain(f) && error("Not an element of the domain")
  return codomain(f)(a)
end

(f::EuclideanRingResidueMap)(a) = image(f, a)

function preimage(f::EuclideanRingResidueMap, a)
  parent(a) != codomain(f) && error("Not an element of the codomain")
  return lift(a)
end

###############################################################################
#
#   Some functions for residue rings of polynomial rings
#
###############################################################################

function gen(R::Union{EuclideanRingResidueRing{T}, EuclideanRingResidueField{T}}) where {T<:PolyRingElem}
   return R(gen(base_ring(R)))
end

# TODO: the names `gen` (algebra generator) and `gens` (module generators) are
# very unfortunate
function gens(R::Union{EuclideanRingResidueRing{T}, EuclideanRingResidueField{T}}) where {T<:PolyRingElem} ## probably needs more cases
   ## as the other residue functions
   g = gen(R)
   r = Vector{typeof(g)}()
   push!(r, one(R))
   if degree(modulus(R)) == 1
      return r
   end
   push!(r, g)
   for i = 2:degree(modulus(R))-1
      push!(r, r[end] * g)
   end
   return r
end

function characteristic(R::Union{EuclideanRingResidueRing{T}, EuclideanRingResidueField{T}}) where {T<:PolyRingElem}
   return characteristic(base_ring(R))
end
is_known(::typeof(characteristic), R::Union{EuclideanRingResidueRing{T}, EuclideanRingResidueField{T}}) where {T<:PolyRingElem} = is_known(characteristic, base_ring(R))

function size(R::Union{EuclideanRingResidueRing{T}, EuclideanRingResidueField{T}}) where {T<:PolyRingElem}
   return size(base_ring(base_ring(R)))^degree(modulus(R))
end

function rand(R::Union{EuclideanRingResidueRing{T}, EuclideanRingResidueField{T}}) where {T<:PolyRingElem}
   r = rand(base_ring(base_ring(R)))
   g = gen(R)
   for i = 1:degree(modulus(R))
      r = r * g + rand(base_ring(base_ring(R)))
   end
   return r
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::T) where {T <: EuclideanRingResidueRingElem}
   a.data = zero!(a.data)
   return a
end

function mul!(c::T, a::T, b::T) where {T <: EuclideanRingResidueRingElem}
   c.data = mod(data(a)*data(b), modulus(a))
   return c
end

function add!(c::T, a::T, b::T) where {T <: EuclideanRingResidueRingElem}
   c.data = mod(data(a) + data(b), modulus(a))
   return c
end

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(R::EuclideanRingResidueRing)
  return R(ConformanceTests.generate_element(base_ring(R)))
end
