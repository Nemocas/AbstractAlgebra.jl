###############################################################################
#
#   Fraction.jl : generic fraction fields
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{Frac{T}}) where T <: RingElem = FracField{T}

elem_type(::Type{FracField{T}}) where {T <: RingElem} = Frac{T}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.numerator(a::Frac, canonicalise::Bool=true)
   if canonicalise
      u = canonical_unit(a.den)
      return divexact(a.num, u)
   else
      return a.num
   end
end

function Base.denominator(a::Frac, canonicalise::Bool=true)
   if canonicalise
      u = canonical_unit(a.den)
      return divexact(a.den, u)
   else
      return a.den
   end
end

function deepcopy_internal(a::Frac{T}, dict::IdDict) where {T <: RingElem}
   v = Frac{T}(deepcopy(numerator(a, false)), deepcopy(denominator(a, false)))
   v.parent = parent(a)
   return v
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{Frac{T}}, ::Type{Frac{T}}) where T <: RingElement = Frac{T}
promote_rule(::Type{Frac{T}}, ::Type{Frac{T}}) where T <: RingElem = Frac{T}

function promote_rule(::Type{Frac{T}}, ::Type{U}) where {T <: RingElem, U <: RingElem}
   promote_rule(T, U) == T ? Frac{T} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::FracField{T})(b::RingElement) where {T <: RingElement}
   return a(base_ring(a)(b))
end

function (a::FracField{T})() where {T <: RingElement}
   z = Frac{T}(zero(base_ring(a)), one(base_ring(a)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::T) where {T <: RingElement}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   z = Frac{T}(b, one(base_ring(a)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::T, c::T) where {T <: RingElement}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   parent(c) != base_ring(a) && error("Could not coerce to fraction")
   z = Frac{T}(b, c)
   z.parent = a
   return z
end

function (a::FracField{T})(b::T, c::Union{Integer, Rational, AbstractFloat}) where {T <: RingElement}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   z = Frac{T}(b, base_ring(a)(c))
   z.parent = a
   return z
end

function (a::FracField{T})(b::Union{Integer, Rational, AbstractFloat}, c::T) where {T <: RingElement}
   parent(c) != base_ring(a) && error("Could not coerce to fraction")
   z = Frac{T}(base_ring(a)(b), c)
   z.parent = a
   return z
end

function (a::FracField{T})(b::Union{Integer, Rational, AbstractFloat}) where {T <: RingElement}
   z = Frac{T}(base_ring(a)(b), one(base_ring(a)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::Integer, c::Integer) where {T <: RingElement}
   z = Frac{T}(base_ring(a)(b), base_ring(a)(c))
   z.parent = a
   return z
end

function (a::FracField{T})(b::Frac{T}) where {T <: RingElement}
   a != parent(b) && error("Could not coerce to fraction")
   return b
end

###############################################################################
#
#   FractionField constructor
#
###############################################################################

function FractionField(R::AbstractAlgebra.Ring; cached=true)
   R2 = R
   T = elem_type(R)

   return FracField{T}(R, cached)
end
