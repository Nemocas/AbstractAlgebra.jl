export Fraction, FractionField, num, den, zero, one, gcd, divexact, mul!, addeq!, inv,
       canonical_unit

import Base: convert, zero, one, show, gcd

import Rings: divexact, mul!, addeq!, inv, canonical_unit

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
   c = canonical_unit(den)
   Fraction{Poly{T, S}}(divexact(num, c), divexact(den, c))
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
#   Unary operations
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

=={T}(x::Fraction{T}, y::Fraction{T}) = x.num == y.num && x.den == y.den

=={T}(x::Fraction{T}, y::ZZ) = x.den == 1 && x.num == T(y)

=={T}(x::Fraction{T}, y::Int) = x.den == 1 && x.num == T(y)

=={T}(x::ZZ, y::Fraction{T}) = y.den == 1 && T(x) == y.num

=={T}(x::Int, y::Fraction{T}) = y.den == 1 && T(x) == y.num


###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{T <: Ring}(io::IO, x::Fraction{T})
   if x.den != 1
      print(io, "(")
   end
   print(io, x.num)
   if x.den != 1
      print(io, ")/(", x.den, ")")
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
#   Binary operations and functions
#
###########################################################################################

+{T <: Ring}(a::Fraction{T}, b::Fraction{T}) = (a.num*b.den + b.num*a.den)/(a.den*b.den)

-{T <: Ring}(a::Fraction{T}, b::Fraction{T}) = (a.num*b.den - b.num*a.den)/(a.den*b.den)

function *{T <: Ring}(a::Fraction{T}, b::Fraction{T})
   g1 = gcd(a.num, b.den)
   g2 = gcd(b.num, a.den)
   Fraction{T}(divexact(a.num, g1)*divexact(b.num, g2), divexact(a.den, g2)*divexact(b.den, g1))
end

function /{T <: Ring}(a::Fraction{T}, b::Fraction{T})
   g1 = gcd(a.num, b.num)
   g2 = gcd(b.den, a.den)
   Fraction{T}(divexact(a.num, g1)*divexact(b.den, g2), divexact(a.den, g2)*divexact(b.num, g1))
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
#   Ad hoc binary operations
#
###########################################################################################

function *{T <: Ring}(a::Fraction{T}, b::Int)
   c = T(b)
   g = gcd(a.den, c)
   Fraction{T}(a.num*divexact(c, g), divexact(a.den, g))
end

function *{T <: Ring}(a::Int, b::Fraction{T})
   c = T(a)
   g = gcd(b.den, c)
   Fraction{T}(b.num*divexact(c, g), divexact(b.den, g))
end

function *{T <: Ring}(a::Fraction{T}, b::ZZ)
   c = T(b)
   g = gcd(a.den, c)
   Fraction{T}(a.num*divexact(c, g), divexact(a.den, g))
end

function *{T <: Ring}(a::ZZ, b::Fraction{T})
   c = T(a)
   g = gcd(b.den, c)
   Fraction{T}(b.num*divexact(c, g), divexact(b.den, g))
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

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{T <: Ring}(a::Fraction{T}, b::Int)
   Fraction{T}(a.num^b, a.den^b)
end

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv{T <: Ring}(a::Fraction{T})
   a.num == 0 && throw(DivideError())
   Fraction{T}(a.den, a.num)
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