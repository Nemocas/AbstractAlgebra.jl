###############################################################################
#
#   Fraction.jl : generic fraction fields
#
###############################################################################

export FractionField, Fraction, num, den

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

elem_type{T <: RingElem}(::FractionField{T}) = Fraction{T}

base_ring(a::FractionField) = a.base_ring

base_ring(a::Fraction) = base_ring(parent(a))

parent(a::Fraction) = a.parent

function check_parent(a::Fraction, b::Fraction)
   parent(a) != parent(b) && error("Incompatible rings in fraction field operation")
end

###############################################################################
#
#   Constructors
#
###############################################################################

function //{T <: RingElem}(x::T, y::T)
   y == 0 && throw(DivideError())
   g = gcd(x, y)
   z = Fraction{T}(divexact(x, g), divexact(y, g))
   try
      z.parent = FractionDict[R]
   catch
      z.parent = FractionField(parent(x))
   end
   return z
end

//{T <: RingElem}(x::T, y::Integer) = x//parent(x)(y)

//{T <: RingElem}(x::Integer, y::T) = parent(y)(x)//y

# disambiguation
//{T <: RingElem}(x::Fraction{T}, y::Fraction{T}) = divexact(x, y)

//{T <: RingElem}(x::T, y::Fraction{T}) = parent(y)(x)//y

//{T <: RingElem}(x::Fraction{T}, y::T) = x//parent(x)(y)

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function hash(a::Fraction)
   h = 0x8a30b0d963237dd5
   return h $ hash(a.num) $ hash(a.den)
end

function num(a::Fraction)
   u = canonical_unit(a.den)
   return divexact(a.num, u)
end

function den(a::Fraction)
   u = canonical_unit(a.den)
   return divexact(a.den, u)
end

zero(a::FractionField) = parent(a)(0)

one(a::FractionField) = parent(a)(1)

iszero(a::Fraction) = iszero(num(a))

isone(a::Fraction) = num(a) == den(a)

isunit(a::Fraction) = num(a) != 0

function deepcopy{T <: RingElem}(a::Fraction{T})
   v = Fraction{T}(deepcopy(num(a)), deepcopy(den(a)))
   v.parent = parent(a)
   return v
end 

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::Fraction) = a

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, x::Fraction)
   u = canonical_unit(x.den)
   num = divexact(x.num, u)
   den = divexact(x.den, u);
   if den != 1 && needs_parentheses(num)
      print(io, "(")
   end
   print(io, num)
   if den != 1
      if needs_parentheses(num)
         print(io, ")")
      end
      print(io, "//")
      if needs_parentheses(den)
         print(io, "(")
      end
      print(io, den)
      if needs_parentheses(den)
         print(io, ")")
      end
   end
end

function show(io::IO, a::FractionField)
   print(io, "Fraction field of ", base_ring(a))
end

needs_parentheses(x::Fraction) = x.den == 1 && needs_parentheses(x.num)

is_negative(x::Fraction) = !needs_parentheses(x.num) && is_negative(x.num)

show_minus_one{T <: RingElem}(::Type{Fraction{T}}) = show_minus_one(T)

###############################################################################
#
#   Unary operators
#
###############################################################################

function -{T <: RingElem}(a::Fraction{T})
   z = Fraction{T}(-a.num, a.den)
   z.parent = parent(a)
   return z
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +{T <: RingElem}(a::Fraction{T}, b::Fraction{T})
   check_parent(a, b)
   num = a.num*b.den + b.num*a.den
   den = a.den*b.den
   g = gcd(num, den)
   z = Fraction{T}(divexact(num, g), divexact(den, g))
   z.parent = parent(a)
   return z
end

function -{T <: RingElem}(a::Fraction{T}, b::Fraction{T})
   check_parent(a, b)
   num = a.num*b.den - b.num*a.den
   den = a.den*b.den
   g = gcd(num, den)
   z = Fraction{T}(divexact(num, g), divexact(den, g))
   z.parent = parent(a)
   return z
end

function *{T <: RingElem}(a::Fraction{T}, b::Fraction{T})
   check_parent(a, b)
   g1 = gcd(a.num, b.den)
   g2 = gcd(b.num, a.den)
   num = divexact(a.num, g1)*divexact(b.num, g2)
   den = divexact(a.den, g2)*divexact(b.den, g1)
   z = Fraction{T}(num, den)
   z.parent = parent(a)
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *{T <: RingElem}(a::Fraction{T}, b::Integer)
   c = base_ring(a)(b)
   g = gcd(a.den, c)
   num = a.num*divexact(c, g)
   den = divexact(a.den, g)
   z = Fraction{T}(num, den)
   z.parent = parent(a)
   return z
end

function *{T <: RingElem}(a::Integer, b::Fraction{T})
   c = base_ring(b)(a)
   g = gcd(b.den, c)
   num = b.num*divexact(c, g)
   den = divexact(b.den, g)
   z = Fraction{T}(num, den)
   z.parent = parent(b)
   return z
end

function *{T <: RingElem}(a::Fraction{T}, b::T)
   g = gcd(a.den, b)
   num = a.num*divexact(b, g)
   den = divexact(a.den, g)
   z = Fraction{T}(num, den)
   z.parent = parent(a)
   return z

end

function *{T <: RingElem}(a::T, b::Fraction{T})
   g = gcd(b.den, a)
   num = b.num*divexact(a, g)
   den = divexact(b.den, g)
   z = Fraction{T}(num, den)
   z.parent = parent(b)
   return z
end

function +{T <: RingElem}(a::Fraction{T}, b::Integer)
   num = a.num + a.den*b
   den = a.den
   g = gcd(num, den)
   z = Fraction{T}(divexact(num, g), divexact(den, g))
   z.parent = parent(a)
   return z
end

function -{T <: RingElem}(a::Fraction{T}, b::Integer)
   num = a.num - a.den*b
   den = a.den
   g = gcd(num, den)
   z = Fraction{T}(divexact(num, g), divexact(den, g))
   z.parent = parent(a)
   return z
end

+{T <: RingElem}(a::Integer, b::Fraction{T}) = b + a

function -{T <: RingElem}(a::Integer, b::Fraction{T})
   num = a*b.den - b.num
   den = b.den
   g = gcd(num, den)
   z = Fraction{T}(divexact(num, g), divexact(den, g))
   z.parent = parent(b)
   return z
end

function +{T <: RingElem}(a::Fraction{T}, b::T)
   num = a.num + a.den*b
   den = a.den
   g = gcd(num, den)
   z = Fraction{T}(divexact(num, g), divexact(den, g))
   z.parent = parent(a)
   return z
end

function -{T <: RingElem}(a::Fraction{T}, b::T)
   num = a.num - a.den*b
   den = a.den
   g = gcd(num, den)
   z = Fraction{T}(divexact(num, g), divexact(den, g))
   z.parent = parent(a)
   return z
end

+{T <: RingElem}(a::T, b::Fraction{T}) = b + a

function -{T <: RingElem}(a::T, b::Fraction{T})
   num = a*b.den - b.num
   den = b.den
   g = gcd(num, den)
   z = Fraction{T}(divexact(num, g), divexact(den, g))
   z.parent = parent(b)
   return z
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function =={T <: RingElem}(x::Fraction{T}, y::Fraction{T})
   check_parent(x, y)
   return (x.den == y.den && x.num == y.num) || (x.num*y.den == x.den*y.num)
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

function ==(x::Fraction, y::Integer)
   return (x.den == 1 && x.num == y) || (x.num == x.den*y)
end

==(x::Integer, y::Fraction) = y == x

function =={T <: RingElem}(x::Fraction{T}, y::T)
   return (x.den == 1 && x.num == y) || (x.num == x.den*y)
end

=={T <: RingElem}(x::T, y::Fraction{T}) = y == x


###############################################################################
#
#   Inversion
#
###############################################################################

function inv{T <: RingElem}(a::Fraction{T})
   a.num == 0 && throw(DivideError())
   z = Fraction{T}(a.den, a.num)
   z.parent = parent(a)
   return z
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact{T <: RingElem}(a::Fraction{T}, b::Fraction{T})
   check_parent(a, b)
   g1 = gcd(a.num, b.num)
   g2 = gcd(b.den, a.den)
   num = divexact(a.num, g1)*divexact(b.den, g2)
   den = divexact(a.den, g2)*divexact(b.num, g1)
   z = Fraction{T}(num, den)
   z.parent = parent(a)
   return z
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact{T <: RingElem}(a::Fraction{T}, b::Integer)
   b == 0 && throw(DivideError())
   c = base_ring(a)(b)
   g = gcd(a.num, c)
   num = divexact(a.num, g)
   den = a.den*divexact(c, g)
   z = Fraction{T}(num, den)
   z.parent = parent(a)
   return z
end

function divexact{T <: RingElem}(a::Integer, b::Fraction{T})
   b == 0 && throw(DivideError())
   c = base_ring(b)(a)
   g = gcd(b.num, c)
   num = b.den*divexact(c, g)
   den = divexact(b.num, g)
   z = Fraction{T}(num, den)
   z.parent = parent(b)
   return z
end

# remove ambiguity
divexact{T <: RingElem}(a::Fraction{T}, b::Poly{T}) = error("Not supported")

# remove ambiguity
divexact{T <: RingElem}(a::Poly{T}, b::Fraction{T}) = error("Not supported")

function divexact{T <: RingElem}(a::Fraction{T}, b::T)
   b == 0 && throw(DivideError())
   g = gcd(a.num, b)
   num = divexact(a.num, g)
   den = a.den*divexact(b, g)
   z = Fraction{T}(num, den)
   z.parent = parent(a)
   return z
end

function divexact{T <: RingElem}(a::T, b::Fraction{T})
   b == 0 && throw(DivideError())
   g = gcd(b.num, a)
   num = b.den*divexact(a, g)
   den = divexact(b.num, g)
   z = Fraction{T}(num, den)
   z.parent = parent(b)
   return z
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^{T <: RingElem}(a::Fraction{T}, b::Int)
   if b < 0
      a = inv(a)
      b = -b
   end
   z = Fraction{T}(a.num^b, a.den^b)
   z.parent = parent(a)
   return z
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd{T <: RingElem}(a::Fraction{T}, b::Fraction{T})
   check_parent(a, b)
   num = gcd(a.num*b.den, a.den*b.num)
   den = a.den*b.den
   g = gcd(num, den)
   z = Fraction{T}(divexact(num, g), divexact(den, g))
   z.parent = parent(a)
   return z
end

###############################################################################
#
#   Unsafe operators and functions
#
###############################################################################

function mul!{T <: RingElem}(c::Fraction{T}, a::Fraction{T}, b::Fraction{T})
   g1 = gcd(a.num, b.den)
   g2 = gcd(b.num, a.den)
   c.num = divexact(a.num, g1)*divexact(b.num, g2)
   c.den = divexact(a.den, g2)*divexact(b.den, g1)
end

function addeq!{T <: RingElem}(c::Fraction{T}, a::Fraction{T})
   num = c.num*a.den + a.num*c.den
   den = c.den*a.den
   g = gcd(num, den)
   c.num = divexact(num, g)
   c.den = divexact(den, g)
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

Base.promote_rule{T <: RingElem}(::Type{Fraction{T}}, ::Type{T}) = Fraction{T}

Base.promote_rule{T <: RingElem, U <: Integer}(::Type{Fraction{T}}, ::Type{U}) = Fraction{T}

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function Base.call{T <: RingElem}(a::FractionField{T}, b::RingElem)
   return a(base_ring(a)(b))
end

function Base.call{T <: RingElem}(a::FractionField{T})
   z = Fraction{T}(zero(base_ring(a)), one(base_ring(a)))
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::FractionField{T}, b::T)
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   z = Fraction{T}(b, one(base_ring(a)))
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::FractionField{T}, b::Integer)
   z = Fraction{T}(base_ring(a)(b), one(base_ring(a)))
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::FractionField{T}, b::Fraction{T})
   a != parent(b) && error("Could not coerce to fraction")
   return b
end

###############################################################################
#
#   FractionField constructor
#
###############################################################################

function FractionField(R::Ring)
   R2 = R
   T = elem_type(R)
   parent_type = Fraction{T}
   while base_ring(R2) != None
      R2 = base_ring(R2)
      T2 = elem_type(R2)
      eval(:(Base.promote_rule(::Type{$parent_type}, ::Type{$T2}) = $parent_type))
   end

   return FractionField{T}(R)
end

