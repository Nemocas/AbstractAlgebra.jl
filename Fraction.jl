export Fraction, FractionField, num, den, zero, one, gcd, divexact, mul!, addeq!, inv,
       canonical_unit, mod, divrem, needs_parentheses, is_negative, show_minus_one

import Base: convert, zero, one, show, gcd

import Rings: divexact, mul!, addeq!, inv, canonical_unit, mod, divrem, needs_parentheses,
              is_negative, show_minus_one

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

type Fraction{T <: Ring} <: Field
   num :: T
   den :: T

   Fraction(a :: T, b :: T) = new(a, b)

   Fraction() = Fraction{T}(zero(T), one(T))
   Fraction(a::Integer) = Fraction{T}(T(a), one(T))
   Fraction(a::T) = Fraction{T}(a, one(T))
   Fraction(a::Fraction{T}) = a
   Fraction{R <: Ring}(a::R) = Fraction{T}(convert(T, a), one(T))
end

###########################################################################################
#
#   Constructors
#
###########################################################################################

function /(x::Int, y::Int) 
   y == 0 && throw(DivideError())
   g = gcd(x, y)
   if y < 0
      Fraction{ZZ}(divexact(ZZ(-x), ZZ(g)), divexact(ZZ(-y), ZZ(g)))
   else
      Fraction{ZZ}(divexact(ZZ(x), ZZ(g)), divexact(ZZ(y), ZZ(g)))
   end
end

function /(x::ZZ, y::ZZ) 
   y == 0 && throw(DivideError())
   g = gcd(x, y)
   if y < 0
      Fraction{ZZ}(divexact(-x, g), divexact(-y, g))
   else
      Fraction{ZZ}(divexact(x, g), divexact(y, g))
   end
end

function /{T <: Ring, S}(x::Poly{T, S}, y::Poly{T, S})
   y == 0 && throw(DivideError())
   g = gcd(x, y)
   num = divexact(x, g)
   den = divexact(y, g)
   u = canonical_unit(den)
   Fraction{Poly{T, S}}(divexact(num, u), divexact(den, u))
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function num{T <: Ring}(a::Fraction{T})
   return a.num
end

function den{T <: Ring}(a::Fraction{T})
   return a.den
end

zero{T <: Ring}(::Type{Fraction{T}}) = Fraction{T}(0)

one{T <: Ring}(::Type{Fraction{T}}) = Fraction{T}(1)

###########################################################################################
#
#   Unary operators
#
###########################################################################################

function -{T <: Ring}(a::Fraction{T})
   Fraction{T}(-a.num, a.den)
end

###########################################################################################
#
#   Comparisons
#
###########################################################################################

=={T <: Ring}(x::Fraction{T}, y::Fraction{T}) = x.num == y.num && x.den == y.den

=={T <: Ring}(x::Fraction{T}, y::T) = x.num == y && x.den == 1

=={T <: Ring}(x::T, y::Fraction{T}) = y.num == x && y.den == 1

=={T <: Ring}(x::Fraction{T}, y::ZZ) = x.den == 1 && x.num == T(y)

=={T <: Ring}(x::Fraction{T}, y::Int) = x.den == 1 && x.num == T(y)

=={T <: Ring}(x::ZZ, y::Fraction{T}) = y.den == 1 && T(x) == y.num

=={T <: Ring}(x::Int, y::Fraction{T}) = y.den == 1 && T(x) == y.num


###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{T <: Ring}(io::IO, x::Fraction{T})
   if x.den != 1 && needs_parentheses(x.num)
      print(io, "(")
   end
   print(io, x.num)
   if x.den != 1
      if needs_parentheses(x.num)
         print(io, ")")
      end
      print(io, "/")
      if needs_parentheses(x.den)
         print(io, "(")
      end
      print(io, x.den)
      if needs_parentheses(x.den)
         print(io, ")")
      end
   end
end

function show(io::IO, x::Fraction{ZZ})
   print(io, x.num)
   if x.den != 1
      print(io, "/", x.den)
   end
end

function show{T <: Ring}(io::IO, ::Type{Fraction{T}})
   print(io, "Fraction field of ")
   show(io, T)
end

needs_parentheses{T <: Ring}(x::Fraction{T}) = false

is_negative{T <: Ring}(x::Fraction{T}) = !needs_parentheses(x.num) && is_negative(x.num)

show_minus_one{T <: Ring}(::Type{Fraction{T}}) = show_minus_one(T)

###########################################################################################
#
#   Conversions
#
###########################################################################################

Base.convert{T <: Ring}(::Type{Fraction{T}}, a::T) = Fraction{T}(a)

Base.convert{T <: Ring}(::Type{Fraction{T}}, a::Int) = Fraction{T}(a)

###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit{T}(a::Fraction{T}) = a

###########################################################################################
#
#   Binary operators and functions
#
###########################################################################################

+{T <: Ring}(a::Fraction{T}, b::Fraction{T}) = (a.num*b.den + b.num*a.den)/(a.den*b.den)

-{T <: Ring}(a::Fraction{T}, b::Fraction{T}) = (a.num*b.den - b.num*a.den)/(a.den*b.den)

function *{T <: Ring}(a::Fraction{T}, b::Fraction{T})
   g1 = gcd(a.num, b.den)
   g2 = gcd(b.num, a.den)
   num = divexact(a.num, g1)*divexact(b.num, g2)
   den = divexact(a.den, g2)*divexact(b.den, g1)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /{T <: Ring}(a::Fraction{T}, b::Fraction{T})
   g1 = gcd(a.num, b.num)
   g2 = gcd(b.den, a.den)
   num = divexact(a.num, g1)*divexact(b.den, g2)
   den = divexact(a.den, g2)*divexact(b.num, g1)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

divexact{T <: Ring}(a::Fraction{T}, b::Fraction{T}) = a/b

function gcd{T <: Ring}(a::Fraction{T}, b::Fraction{T})
   Fraction{T}(gcd(a.num, b.num), gcd(a.den, b.den))
end

###########################################################################################
#
#   Unsafe operators and functions
#
###########################################################################################

function mul!{T <: Ring}(c::Fraction{T}, a::Fraction{T}, b::Fraction{T})
   g1 = gcd(a.num, b.den)
   g2 = gcd(b.num, a.den)
   c.num = divexact(a.num, g1)*divexact(b.num, g2)
   c.den = divexact(a.den, g2)*divexact(b.den, g1)
end

function addeq!{T <: Ring}(c::Fraction{T}, a::Fraction{T})
   num = c.num*a.den + a.num*c.den
   den = c.den*a.den
   g = gcd(num, den)
   c.num = divexact(num, g)
   c.den = divexact(den, g)
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *{T <: Ring}(a::Fraction{T}, b::Int)
   c = T(b)
   g = gcd(a.den, c)
   num = a.num*divexact(c, g)
   den = divexact(a.den, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function *{T <: Ring}(a::Int, b::Fraction{T})
   c = T(a)
   g = gcd(b.den, c)
   num = b.num*divexact(c, g)
   den = divexact(b.den, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function *{T <: Ring}(a::Fraction{T}, b::ZZ)
   c = T(b)
   g = gcd(a.den, c)
   num = a.num*divexact(c, g)
   den = divexact(a.den, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function *{T <: Ring}(a::ZZ, b::Fraction{T})
   c = T(a)
   g = gcd(b.den, c)
   num = b.num*divexact(c, g)
   den = divexact(b.den, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function *{T <: Ring}(a::Fraction{T}, b::T)
   g = gcd(a.den, b)
   num = a.num*divexact(b, g)
   den = divexact(a.den, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function *{T <: Ring}(a::T, b::Fraction{T})
   g = gcd(b.den, a)
   num = b.num*divexact(a, g)
   den = divexact(b.den, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /{T <: Ring}(a::Fraction{T}, b::Int)
   c = T(b)
   g = gcd(a.num, c)
   num = divexact(a.num, g)
   den = a.den*divexact(c, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /{T <: Ring}(a::Int, b::Fraction{T})
   c = T(a)
   g = gcd(b.num, c)
   num = b.den*divexact(c, g)
   den = divexact(b.num, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /{T <: Ring}(a::Fraction{T}, b::ZZ)
   c = T(b)
   g = gcd(a.num, c)
   num = divexact(a.num, g)
   den = a.den*divexact(c, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /{T <: Ring}(a::ZZ, b::Fraction{T})
   c = T(a)
   g = gcd(b.num, c)
   num = b.den*divexact(c, g)
   den = divexact(b.num, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /{T <: Ring}(a::Fraction{T}, b::T)
   g = gcd(a.num, b)
   num = divexact(a.num, g)
   den = a.den*divexact(b, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /{T <: Ring}(a::T, b::Fraction{T})
   g = gcd(b.num, a)
   num = b.den*divexact(a, g)
   den = divexact(b.num, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function +{T <: Ring}(a::Fraction{T}, b::Int)
   (a.num + a.den*b)/a.den
end

function +{T <: Ring}(a::Int, b::Fraction{T})
   (a*b.den + b.num)/b.den
end

function +{T <: Ring}(a::Fraction{T}, b::ZZ)
   (a.num + a.den*b)/a.den
end

function +{T <: Ring}(a::ZZ, b::Fraction{T})
   (a*b.den + b.num)/b.den
end

function +{T <: Ring}(a::Fraction{T}, b::T)
   (a.num + a.den*b)/a.den
end

function +{T <: Ring}(a::T, b::Fraction{T})
   (a*b.den + b.num)/b.den
end

function -{T <: Ring}(a::Fraction{T}, b::Int)
   (a.num - a.den*b)/a.den
end

function -{T <: Ring}(a::Int, b::Fraction{T})
   (a*b.den - b.num)/b.den
end

function -{T <: Ring}(a::Fraction{T}, b::ZZ)
   (a.num - a.den*b)/a.den
end

function -{T <: Ring}(a::ZZ, b::Fraction{T})
   (a*b.den - b.num)/b.den
end

function -{T <: Ring}(a::Fraction{T}, b::T)
   (a.num - a.den*b)/a.den
end

function -{T <: Ring}(a::T, b::Fraction{T})
   (a*b.den - b.num)/b.den
end

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{T <: Ring}(a::Fraction{T}, b::Int)
   if b < 0
      a = inv(a)
      b = -b
   end
   Fraction{T}(a.num^b, a.den^b)
end

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv{T <: Ring}(a::Fraction{T})
   a.num == 0 && throw(DivideError())
   num = a.den
   den = a.num
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

###########################################################################################
#
#   FractionField constructor
#
###########################################################################################

function FractionField{T <: Ring}(::Type{T})
   return Fraction{T}
end

function FractionField{T <: Fraction}(::Type{T})
   return T # fraction field of a field is itself
end
