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

parent(a::FracFieldElem) = a.parent

parent_type(::Type{FracFieldElem{T}}) where T <: RingElem = FracField{T}

elem_type(::Type{FracField{T}}) where {T <: RingElem} = FracFieldElem{T}

base_ring(a::FracField{T}) where T <: RingElem = a.base_ring::base_ring_type(a)

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

promote_rule(::Type{FracFieldElem{T}}, ::Type{FracFieldElem{T}}) where T <: RingElem = FracFieldElem{T}

function promote_rule(::Type{FracFieldElem{T}}, ::Type{U}) where {T <: RingElem, U <: RingElem}
   promote_rule(T, U) == T ? FracFieldElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::FracField{T})(b::RingElement) where {T <: RingElem}
   return a(base_ring(a)(b))
end

function (a::FracField{T})() where {T <: RingElem}
   z = FracFieldElem{T}(zero(base_ring(a)), one(base_ring(a)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::T) where {T <: RingElem}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   z = FracFieldElem{T}(b, one(base_ring(a)))
   z.parent = a
   return z
end

function _reduce_fraction(x, y)
   iszero(y) && throw(DivideError())
   g = gcd(x, y)
   return divexact(x, g), divexact(y, g)
end

function (a::FracField{T})(b::T, c::T; reduce::Bool = false) where {T <: RingElem}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   parent(c) != base_ring(a) && error("Could not coerce to fraction")
   if reduce
     b, c = _reduce_fraction(b, c)
   end
   z = FracFieldElem{T}(b, c)
   z.parent = a
   return z
end

function (a::FracField{T})(b::T, c::T; reduce::Bool = false) where {U <: FieldElem, T <: PolyRingElem{U}}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   parent(c) != base_ring(a) && error("Could not coerce to fraction")

   if reduce
      b, c = _reduce_fraction(b, c)
   else
      u = canonical_unit(c)
      if !isone(u)
         b = divexact(b, u)
         c = divexact(c, u)
      end
   end
   z = FracFieldElem{T}(b, c)
   z.parent = a
   return z
end

function (a::FracField{T})(b::T, c::JuliaRingElement) where {T <: RingElem}
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

function (a::FracField{T})(b::JuliaRingElement, c::T) where {T <: RingElem}
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

function (a::FracField{T})(b::Union{Integer, AbstractFloat}) where {T <: RingElem}
   z = FracFieldElem{T}(base_ring(a)(b), one(base_ring(a)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::Rational) where {T <: RingElem}
   z = FracFieldElem{T}(base_ring(a)(numerator(b, false)),
               base_ring(a)(denominator(b, false)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::Integer, c::Integer) where {T <: RingElem}
   z = FracFieldElem{T}(base_ring(a)(b), base_ring(a)(c))
   z.parent = a
   return z
end

function (a::FracField{T})(b::FracFieldElem{T}) where {T <: RingElem}
   a != parent(b) && error("Could not coerce to fraction")
   return b
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(c::FracFieldElem)
   c.num = zero!(c.num)
   if !isone(c.den)
      c.den = one(base_ring(c))
   end
   return c
end

function mul!(c::FracFieldElem{T}, a::FracFieldElem{T}, b::FracFieldElem{T}) where {T <: RingElem}
   n1 = numerator(a, false)
   d2 = denominator(b, false)
   n2 = numerator(b, false)
   d1 = denominator(a, false)
   if d1 == d2
      c.num = n1*n2
      c.den = d1*d2
   elseif isone(d1)
      gd = gcd(n1, d2)
      if isone(gd)
         c.num = n1*n2
         c.den = deepcopy(d2)
      else
         c.num = divexact(n1, gd)*n2
         c.den = divexact(d2, gd)
      end
   elseif isone(d2)
      gd = gcd(n2, d1)
      if isone(gd)
         c.num = n2*n1
         c.den = deepcopy(d1)
      else
         c.num = divexact(n2, gd)*n1
         c.den = divexact(d1, gd)
      end
   else
      g1 = gcd(n1, d2)
      g2 = gcd(n2, d1)
      if !isone(g1)
         n1 = divexact(n1, g1)
         d2 = divexact(d2, g1)
      end
      if !isone(g2)
         n2 = divexact(n2, g2)
         d1 = divexact(d1, g2)
      end
      c.num = n1*n2
      c.den = d1*d2
   end
   return c
end

function add!(c::FracFieldElem{T}, a::FracFieldElem{T}, b::FracFieldElem{T}) where {T <: RingElem}
   d1 = denominator(a, false)
   d2 = denominator(b, false)
   n1 = numerator(a, false)
   n2 = numerator(b, false)
   if d1 == d2
      c.num = n1 + n2
      if isone(d1)
         c.den = deepcopy(d1)
      else
         gd = gcd(c.num, d1)
         if isone(gd)
            c.den = deepcopy(d1)
         else
            c.num = divexact(c.num, gd)
            c.den = divexact(d1, gd)
         end
      end
   elseif isone(d1)
      c.num = n1*d2 + n2
      c.den = deepcopy(d2)
   elseif isone(d2)
      c.num = n1 + n2*d1
      c.den = deepcopy(d1)
   else
      gd = gcd(d1, d2)
      if isone(gd)
         c.num = n1*d2 + n2*d1
         c.den = d1*d2
      else
         q1 = divexact(d1, gd)
         q2 = divexact(d2, gd)
         c.num = q1*n2 + q2*n1
         t = gcd(c.num, gd)
         if isone(t)
            c.den = q2*d1
         else
            gd = divexact(d1, t)
            c.num = divexact(c.num, t)
            c.den = gd*q2
         end
      end
   end
   return c
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

function ConformanceTests.generate_element(R::FracField{T}) where {T <: RingElem}
  num = ConformanceTests.generate_element(base_ring(R))
  den = ConformanceTests.generate_element(base_ring(R))
  if iszero(den)
    den = one!(den)
  end
  return R(num, den)
end
