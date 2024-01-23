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

parent_type(::Type{FracFieldElem{T}}) where T <: RingElem = FracField{T}

elem_type(::Type{FracField{T}}) where {T <: RingElem} = FracFieldElem{T}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.numerator(a::FracFieldElem, canonicalise::Bool=true)
   if canonicalise
      u = canonical_unit(a.den)
      return divexact(a.num, u)
   else
      return a.num
   end
end

function Base.denominator(a::FracFieldElem, canonicalise::Bool=true)
   if canonicalise
      u = canonical_unit(a.den)
      return divexact(a.den, u)
   else
      return a.den
   end
end

function deepcopy_internal(a::FracFieldElem{T}, dict::IdDict) where {T <: RingElem}
   v = FracFieldElem{T}(deepcopy_internal(numerator(a, false), dict),
               deepcopy_internal(denominator(a, false), dict))
   v.parent = parent(a)
   return v
end

number_of_generators(F::FracField) = number_of_generators(base_ring(F))

gen(F::FracField) = F(gen(base_ring(F)))

gen(F::FracField, i::Int) = F(gen(base_ring(F), i))

gens(F::FracField) = F.(gens(base_ring(F)))


###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{FracFieldElem{T}}, ::Type{FracFieldElem{T}}) where T <: RingElement = FracFieldElem{T}
promote_rule(::Type{FracFieldElem{T}}, ::Type{FracFieldElem{T}}) where T <: RingElem = FracFieldElem{T}

function promote_rule(::Type{FracFieldElem{T}}, ::Type{U}) where {T <: RingElem, U <: RingElem}
   promote_rule(T, U) == T ? FracFieldElem{T} : Union{}
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
   z = FracFieldElem{T}(zero(base_ring(a)), one(base_ring(a)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::T) where {T <: RingElement}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   z = FracFieldElem{T}(b, one(base_ring(a)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::T, c::T) where {T <: RingElement}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   parent(c) != base_ring(a) && error("Could not coerce to fraction")
   z = FracFieldElem{T}(b, c)
   z.parent = a
   return z
end

function (a::FracField{T})(b::T, c::T) where {U <: FieldElem, T <: PolyRingElem{U}}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   parent(c) != base_ring(a) && error("Could not coerce to fraction")
   u = canonical_unit(c)
   if !isone(u)
      b = divexact(b, u)
      c = divexact(c, u)
   end
   z = FracFieldElem{T}(b, c)
   z.parent = a
   return z
end

function (a::FracField{T})(b::T, c::Union{Integer, Rational, AbstractFloat}) where {T <: RingElement}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   z = FracFieldElem{T}(b, base_ring(a)(c))
   z.parent = a
   return z
end

function (a::FracField{T})(b::T, c::Rational) where {U <: FieldElem, T <: PolyRingElem{U}}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   b *= inv(c)
   z = FracFieldElem{T}(b, one(base_ring(a)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::Union{Integer, Rational, AbstractFloat}, c::T) where {T <: RingElement}
   parent(c) != base_ring(a) && error("Could not coerce to fraction")
   z = FracFieldElem{T}(base_ring(a)(b), c)
   z.parent = a
   return z
end

function (a::FracField{T})(b::Union{Integer, Rational}, c::T) where {U <: FieldElem, T <: PolyRingElem{U}}
   parent(c) != base_ring(a) && error("Could not coerce to fraction")
   b = base_ring(a)(b)
   u = canonical_unit(c)
   if !isone(u)
      b = divexact(b, u)
      c = divexact(c, u)
   end
   z = FracFieldElem{T}(b, c)
   z.parent = a
   return z
end

function (a::FracField{T})(b::Union{Integer, AbstractFloat}) where {T <: RingElement}
   z = FracFieldElem{T}(base_ring(a)(b), one(base_ring(a)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::Rational) where {T <: RingElement}
   z = FracFieldElem{T}(base_ring(a)(numerator(b, false)),
               base_ring(a)(denominator(b, false)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::Integer, c::Integer) where {T <: RingElement}
   z = FracFieldElem{T}(base_ring(a)(b), base_ring(a)(c))
   z.parent = a
   return z
end

function (a::FracField{T})(b::FracFieldElem{T}) where {T <: RingElement}
   a != parent(b) && error("Could not coerce to fraction")
   return b
end

###############################################################################
#
#   fraction_field constructor
#
###############################################################################

function fraction_field(R::AbstractAlgebra.Ring; cached::Bool=true)
   T = elem_type(R)

   return FracField{T}(R, cached)
end
