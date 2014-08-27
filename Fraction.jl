export Fraction, FractionField, num, den, zero, one, gcd, divexact, mul!, addeq!, inv,
       canonical_unit, mod, divrem, needs_parentheses, is_negative, show_minus_one, QQ,
       cmp, height, height_bits, reconstruct, harmonic, dedekind_sum, next_minimal,
       next_signed_minimal, next_calkin_wilf, next_signed_calkin_wilf, isunit

import Base: convert, show, gcd, string

import Rings: divexact, mul!, addeq!, inv, canonical_unit, mod, divrem, needs_parentheses,
              is_negative, show_minus_one, zero, one, isunit

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

type fmpq
   num :: Int # can't make these a ZZ
   den :: Int

   function fmpq(a::ZZ, b::ZZ)
      d = new(0, 0)
      ccall((:fmpq_init, :libflint), Void, (Ptr{fmpq},), &d)
      ccall((:fmpq_set_fmpz_frac, :libflint), Void, (Ptr{fmpq}, Ptr{ZZ}, Ptr{ZZ}), &d, &a, &b)
      finalizer(d, _fmpq_clear_fn)
      return d
   end
end

function _fmpq_clear_fn(a::fmpq)
   ccall((:fmpq_clear, :libflint), Void, (Ptr{fmpq},), &a)
end

type FractionStruct{T <: Ring}
   num :: T
   den :: T
end

type Fraction{T <: Ring} <: Field
   data :: Union(FractionStruct, fmpq)

   Fraction(a :: T, b :: T) = new(FractionStruct(a, b))
   Fraction(a :: ZZ, b :: ZZ) = new(fmpq(a, b))

   Fraction() = Fraction{T}(zero(T), one(T))
   Fraction(a::Integer) = Fraction{T}(T(a), one(T))
   Fraction(a::T) = Fraction{T}(a, one(T))
   Fraction(a::Fraction{T}) = a
   Fraction{R <: Ring}(a::R) = Fraction{T}(convert(T, a), one(T))
end

typealias QQ Fraction{ZZ}

###########################################################################################
#
#   Constructors
#
###########################################################################################

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
   return a.data.num
end

function num(a::Fraction{ZZ})
   c = ZZ()
   ccall((:fmpq_numerator, :libflint), Void, (Ptr{ZZ}, Ptr{fmpq}), &c, &(a.data))
   return c
end

function den{T <: Ring}(a::Fraction{T})
   return a.data.den
end

function den(a::Fraction{ZZ})
   c = ZZ()
   ccall((:fmpq_denominator, :libflint), Void, (Ptr{ZZ}, Ptr{fmpq}), &c, &(a.data))
   return c
end

zero{T <: Ring}(::Type{Fraction{T}}) = Fraction{T}(0)

one{T <: Ring}(::Type{Fraction{T}}) = Fraction{T}(1)

isunit{T <: Ring}(a::Fraction{T}) = num(a) != 0

function height(a::Fraction{ZZ})
   c = ZZ()
   ccall((:fmpq_height, :libflint), Void, (Ptr{ZZ}, Ptr{fmpq}), &c, &(a.data))
   return c
end

function height_bits(a::Fraction{ZZ})
   return ccall((:fmpq_height_bits, :libflint), Int, (Ptr{fmpq},), &(a.data))
end

###########################################################################################
#
#   Unary operators
#
###########################################################################################

function -{T <: Ring}(a::Fraction{T})
   Fraction{T}(-a.data.num, a.data.data.den)
end

function -(a::Fraction{ZZ})
   c = Fraction{ZZ}()
   ccall((:fmpq_neg, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(a.data))
   return c
end

###########################################################################################
#
#   Comparisons
#
###########################################################################################

=={T <: Ring}(x::Fraction{T}, y::Fraction{T}) = x.data.num == y.data.num && x.data.data.den == y.data.data.den

function ==(a::Fraction{ZZ}, b::Fraction{ZZ})
   return bool(ccall((:fmpq_equal, :libflint), Cint, (Ptr{fmpq}, Ptr{fmpq}), &(a.data), &(b.data)))
end

=={T <: Ring}(x::Fraction{T}, y::T) = x.data.num == y && x.data.den == 1

=={T <: Ring}(x::Fraction{T}, y::ZZ) = x.data.den == 1 && x.data.num == T(y)

=={T <: Ring}(x::Fraction{T}, y::Int) = x.data.den == 1 && x.data.num == T(y)

function ==(a::Fraction{ZZ}, b::ZZ)
   return bool(ccall((:fmpq_equal_fmpz, :libflint), Cint, (Ptr{fmpq}, Ptr{ZZ}), &(a.data), &b))
end

function ==(a::Fraction{ZZ}, b::Int)
   return bool(ccall((:fmpq_equal_si, :libflint), Cint, (Ptr{fmpq}, Int), &(a.data), b))
end

=={T <: Ring}(x::T, y::Fraction{T}) = y.data.num == x && y.data.den == 1

=={T <: Ring}(x::ZZ, y::Fraction{T}) = y.data.den == 1 && T(x) == y.data.num

=={T <: Ring}(x::Int, y::Fraction{T}) = y.data.den == 1 && T(x) == y.data.num

==(a::ZZ, b::Fraction{ZZ}) = b == a

==(a::Int, b::Fraction{ZZ}) = b == ZZ(a)

function cmp(a::Fraction{ZZ}, b::Fraction{ZZ})
   return int(ccall((:fmpq_cmp, :libflint), Cint, (Ptr{fmpq}, Ptr{fmpq}), &(a.data), &(b.data)))
end

<(a::Fraction{ZZ}, b::Fraction{ZZ}) = cmp(a, b) < 0

<(a::Fraction{ZZ}, b::ZZ) = cmp(a, Fraction{ZZ}(b)) < 0

<(a::Fraction{ZZ}, b::Int) = cmp(a, Fraction{ZZ}(b)) < 0

<(a::ZZ, b::Fraction{ZZ}) = cmp(Fraction{ZZ}(a), b) < 0

<(a::Int, b::Fraction{ZZ}) = cmp(Fraction{ZZ}(a), b) < 0

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{T <: Ring}(io::IO, x::Fraction{T})
   if x.data.den != 1 && needs_parentheses(x.data.num)
      print(io, "(")
   end
   print(io, x.data.num)
   if x.data.den != 1
      if needs_parentheses(x.data.num)
         print(io, ")")
      end
      print(io, "/")
      if needs_parentheses(x.data.den)
         print(io, "(")
      end
      print(io, x.data.den)
      if needs_parentheses(x.data.den)
         print(io, ")")
      end
   end
end

function show(io::IO, a::Fraction{ZZ})
   p = ccall((:fmpq_get_str,:libflint), Ptr{Uint8}, (Ptr{Uint8}, Int, Ptr{fmpq}), C_NULL, 10, &(a.data))
   len = int(ccall(:strlen, Csize_t, (Ptr{Uint8},), p))
   print(io, ASCIIString(pointer_to_array(p, len, true)))
end

function show{T <: Ring}(io::IO, ::Type{Fraction{T}})
   print(io, "Fraction field of ")
   show(io, T)
end

needs_parentheses{T <: Ring}(x::Fraction{T}) = x.data.den == 1 && needs_parentheses(x.data.num)

needs_parentheses(x::Fraction{ZZ}) = false

is_negative{T <: Ring}(x::Fraction{T}) = !needs_parentheses(x.data.num) && is_negative(x.data.num)

is_negative(x::Fraction{ZZ}) = x < 0

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

+{T <: Ring}(a::Fraction{T}, b::Fraction{T}) = (a.data.num*b.data.den + b.data.num*a.data.den)/(a.data.den*b.data.den)

function +(a::Fraction{ZZ}, b::Fraction{ZZ})
   c = Fraction{ZZ}()
   ccall((:fmpq_add, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(a.data), &(b.data))
   return c
end

-{T <: Ring}(a::Fraction{T}, b::Fraction{T}) = (a.data.num*b.data.den - b.data.num*a.data.den)/(a.data.den*b.data.den)

function -(a::Fraction{ZZ}, b::Fraction{ZZ})
   c = Fraction{ZZ}()
   ccall((:fmpq_sub, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(a.data), &(b.data))
   return c
end

function *{T <: Ring}(a::Fraction{T}, b::Fraction{T})
   g1 = gcd(a.data.num, b.data.den)
   g2 = gcd(b.data.num, a.data.den)
   num = divexact(a.data.num, g1)*divexact(b.data.num, g2)
   den = divexact(a.data.den, g2)*divexact(b.data.den, g1)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function *(a::Fraction{ZZ}, b::Fraction{ZZ})
   c = Fraction{ZZ}()
   ccall((:fmpq_mul, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(a.data), &(b.data))
   return c
end

function /{T <: Ring}(a::Fraction{T}, b::Fraction{T})
   g1 = gcd(a.data.num, b.data.num)
   g2 = gcd(b.data.den, a.data.den)
   num = divexact(a.data.num, g1)*divexact(b.data.den, g2)
   den = divexact(a.data.den, g2)*divexact(b.data.num, g1)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /(a::Fraction{ZZ}, b::Fraction{ZZ})
   c = Fraction{ZZ}()
   ccall((:fmpq_div, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(a.data), &(b.data))
   return c
end

divexact{T <: Ring}(a::Fraction{T}, b::Fraction{T}) = a/b

function gcd{T <: Ring}(a::Fraction{T}, b::Fraction{T})
   gcd(a.data.num*b.data.den, a.data.den*b.data.num)/(a.data.den*b.data.den)
end

function gcd(a::Fraction{ZZ}, b::Fraction{ZZ})
   c = Fraction{ZZ}()
   ccall((:fmpq_gcd, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(a.data), &(b.data))
   return c
end

###########################################################################################
#
#   Unsafe operators and functions
#
###########################################################################################

function mul!{T <: Ring}(c::Fraction{T}, a::Fraction{T}, b::Fraction{T})
   g1 = gcd(a.data.num, b.data.den)
   g2 = gcd(b.data.num, a.data.den)
   c.data.num = divexact(a.data.num, g1)*divexact(b.data.num, g2)
   c.data.den = divexact(a.data.den, g2)*divexact(b.data.den, g1)
end

function mul!(c::Fraction{ZZ}, a::Fraction{ZZ}, b::Fraction{ZZ})
   ccall((:fmpq_mul, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(a.data), &(b.data))
end

function addeq!{T <: Ring}(c::Fraction{T}, a::Fraction{T})
   num = c.data.num*a.data.den + a.data.num*c.data.den
   den = c.data.den*a.data.den
   g = gcd(num, den)
   c.data.num = divexact(num, g)
   c.data.den = divexact(den, g)
end

function addeq!(c::Fraction{ZZ}, a::Fraction{ZZ})
   ccall((:fmpq_add, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(c.data), &(a.data))
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *{T <: Ring}(a::Fraction{T}, b::Int)
   c = T(b)
   g = gcd(a.data.den, c)
   num = a.data.num*divexact(c, g)
   den = divexact(a.data.den, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function *{T <: Ring}(a::Fraction{T}, b::ZZ)
   c = T(b)
   g = gcd(a.data.den, c)
   num = a.data.num*divexact(c, g)
   den = divexact(a.data.den, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function *(a::Fraction{ZZ}, b::ZZ)
   c = Fraction{ZZ}()
   ccall((:fmpq_mul_fmpz, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{ZZ}), &(c.data), &(a.data), &b)
   return c
end

*(a::Fraction{ZZ}, b::Int) = a*ZZ(b)

function *{T <: Ring}(a::Int, b::Fraction{T})
   c = T(a)
   g = gcd(b.data.den, c)
   num = b.data.num*divexact(c, g)
   den = divexact(b.data.den, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function *{T <: Ring}(a::ZZ, b::Fraction{T})
   c = T(a)
   g = gcd(b.data.den, c)
   num = b.data.num*divexact(c, g)
   den = divexact(b.data.den, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

*(a::ZZ, b::Fraction{ZZ}) = b*a

*(a::Int, b::Fraction{ZZ}) = b*ZZ(a)

function *{T <: Ring}(a::Fraction{T}, b::T)
   g = gcd(a.data.den, b)
   num = a.data.num*divexact(b, g)
   den = divexact(a.data.den, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function *{T <: Ring}(a::T, b::Fraction{T})
   g = gcd(b.data.den, a)
   num = b.data.num*divexact(a, g)
   den = divexact(b.data.den, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /{T <: Ring}(a::Fraction{T}, b::Int)
   c = T(b)
   g = gcd(a.data.num, c)
   num = divexact(a.data.num, g)
   den = a.data.den*divexact(c, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /{T <: Ring}(a::Fraction{T}, b::ZZ)
   c = T(b)
   g = gcd(a.data.num, c)
   num = divexact(a.data.num, g)
   den = a.data.den*divexact(c, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /(a::Fraction{ZZ}, b::ZZ)
   c = Fraction{ZZ}()
   ccall((:fmpq_div_fmpz, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{ZZ}), &(c.data), &(a.data), &b)
   return c
end

/(a::Fraction{ZZ}, b::Int) = a/ZZ(b)

function /{T <: Ring}(a::ZZ, b::Fraction{T})
   c = T(a)
   g = gcd(b.data.num, c)
   num = b.data.den*divexact(c, g)
   den = divexact(b.data.num, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /{T <: Ring}(a::Int, b::Fraction{T})
   c = T(a)
   g = gcd(b.data.num, c)
   num = b.data.den*divexact(c, g)
   den = divexact(b.data.num, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

/(a::ZZ, b::Fraction{ZZ}) = inv(b)*a

/(a::Int, b::Fraction{ZZ}) = inv(b)*ZZ(a)

function /{T <: Ring}(a::Fraction{T}, b::T)
   g = gcd(a.data.num, b)
   num = divexact(a.data.num, g)
   den = a.data.den*divexact(b, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function /{T <: Ring}(a::T, b::Fraction{T})
   g = gcd(b.data.num, a)
   num = b.data.den*divexact(a, g)
   den = divexact(b.data.num, g)
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function +{T <: Ring}(a::Fraction{T}, b::Int)
   (a.data.num + a.data.den*b)/a.data.den
end

function +{T <: Ring}(a::Fraction{T}, b::ZZ)
   (a.data.num + a.data.den*b)/a.data.den
end

function +(a::Fraction{ZZ}, b::Int)
   c = Fraction{ZZ}()
   ccall((:fmpq_add_si, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Int), &(c.data), &(a.data), b)
   return c
end

function +(a::Fraction{ZZ}, b::ZZ)
   c = Fraction{ZZ}()
   ccall((:fmpq_add_fmpz, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{ZZ}), &(c.data), &(a.data), &b)
   return c
end

function +{T <: Ring}(a::Int, b::Fraction{T})
   (a*b.data.den + b.data.num)/b.data.den
end

function +{T <: Ring}(a::ZZ, b::Fraction{T})
   (a*b.data.den + b.data.num)/b.data.den
end

+(a::ZZ, b::Fraction{ZZ}) = b + a

+(a::Int, b::Fraction{ZZ}) = b + a

function +{T <: Ring}(a::Fraction{T}, b::T)
   (a.data.num + a.data.den*b)/a.data.den
end

function +{T <: Ring}(a::T, b::Fraction{T})
   (a*b.data.den + b.data.num)/b.data.den
end

function -{T <: Ring}(a::Fraction{T}, b::Int)
   (a.data.num - a.data.den*b)/a.data.den
end

function -{T <: Ring}(a::Fraction{T}, b::ZZ)
   (a.data.num - a.data.den*b)/a.data.den
end

function -(a::Fraction{ZZ}, b::Int)
   c = Fraction{ZZ}()
   ccall((:fmpq_sub_si, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Int), &(c.data), &(a.data), b)
   return c
end

function -(a::Fraction{ZZ}, b::ZZ)
   c = Fraction{ZZ}()
   ccall((:fmpq_sub_fmpz, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{ZZ}), &(c.data), &(a.data), &b)
   return c
end

function -{T <: Ring}(a::Int, b::Fraction{T})
   (a*b.data.den - b.data.num)/b.data.den
end

function -{T <: Ring}(a::ZZ, b::Fraction{T})
   (a*b.data.den - b.data.num)/b.data.den
end

function -(a::Int, b::Fraction{ZZ})
   c = Fraction{ZZ}()
   ccall((:fmpq_sub_si, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Int), &(c.data), &(b.data), a)
   ccall((:fmpq_neg, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(c.data))
   return c
end

function -(a::ZZ, b::Fraction{ZZ})
   c = Fraction{ZZ}()
   ccall((:fmpq_sub_fmpz, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{ZZ}), &(c.data), &(b.data), &a)
   ccall((:fmpq_neg, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(c.data))
   return c
end

function -{T <: Ring}(a::Fraction{T}, b::T)
   (a.data.num - a.data.den*b)/a.data.den
end

function -{T <: Ring}(a::T, b::Fraction{T})
   (a*b.data.den - b.data.num)/b.data.den
end

function >>(a::Fraction{ZZ}, b::Int)
   c = Fraction{ZZ}()
   ccall((:fmpq_div_2exp, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Int), &(c.data), &(a.data), b)
   return c
end

function <<(a::Fraction{ZZ}, b::Int)
   c = Fraction{ZZ}()
   ccall((:fmpq_mul_2exp, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Int), &(c.data), &(a.data), b)
   return c
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
   Fraction{T}(a.data.num^b, a.data.den^b)
end

function ^(a::Fraction{ZZ}, b::Int)
   c = Fraction{ZZ}()
   ccall((:fmpq_pow_si, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Int), &(c.data), &(a.data), b)
   return c
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

divexact(a::Fraction{ZZ}, b::Int) = a/ZZ(b)

divexact(a::Fraction{ZZ}, b::ZZ) = a/b

divexact{T <: Ring}(a::Fraction{T}, b::Int) = a/ZZ(b)

divexact{T <: Ring}(a::Fraction{T}, b::ZZ) = a/b

divexact{T <: Ring}(a::Fraction{T}, b::T) = a/b

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv{T <: Ring}(a::Fraction{T})
   a.data.num == 0 && throw(DivideError())
   num = a.data.den
   den = a.data.num
   u = canonical_unit(den)
   Fraction{T}(divexact(num, u), divexact(den, u))
end

function inv(a::Fraction{ZZ})
   c = Fraction{ZZ}()
   ccall((:fmpq_inv, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(a.data))
   return c
end

###########################################################################################
#
#   Modular arithmetic
#
###########################################################################################

function mod(a::Fraction{ZZ}, b::ZZ)
   c = ZZ()
   ccall((:fmpq_mod_fmpz, :libflint), Void, (Ptr{ZZ}, Ptr{fmpq}, Ptr{ZZ}), &c, &(a.data), &b)
   return c
end

mod(a::Fraction{ZZ}, b::Int) = mod(a, ZZ(b))

###########################################################################################
#
#   Rational reconstruction
#
###########################################################################################

function reconstruct(a::ZZ, b::ZZ)
   c = Fraction{ZZ}()
   if !bool(ccall((:fmpq_reconstruct_fmpz, :libflint), Cint, (Ptr{fmpq}, Ptr{ZZ}, Ptr{ZZ}), &(c.data), &a, &b))
      error("Impossible rational reconstruction")
   end
   return c
end

reconstruct(a::ZZ, b::Int) =  reconstruct(a, ZZ(b))

reconstruct(a::Int, b::ZZ) =  reconstruct(ZZ(a), b)

reconstruct(a::Int, b::Int) =  reconstruct(ZZ(a), ZZ(b))

###########################################################################################
#
#   Rational enumeration
#
###########################################################################################

function next_minimal(a::Fraction{ZZ})
   a < 0 && throw(DomainError())
   c = Fraction{ZZ}()
   ccall((:fmpq_next_minimal, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(a.data))
   return c
end

function next_signed_minimal(a::Fraction{ZZ})
   c = Fraction{ZZ}()
   ccall((:fmpq_next_signed_minimal, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(a.data))
   return c
end

function next_calkin_wilf(a::Fraction{ZZ})
   a < 0 && throw(DomainError())
   c = Fraction{ZZ}()
   ccall((:fmpq_next_calkin_wilf, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(a.data))
   return c
end

function next_signed_calkin_wilf(a::Fraction{ZZ})
   c = Fraction{ZZ}()
   ccall((:fmpq_next_signed_calkin_wilf, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &(c.data), &(a.data))
   return c
end

###########################################################################################
#
#   Special functions
#
###########################################################################################

function harmonic(n::Int)
   n < 0 && throw(DomainError())
   c = Fraction{ZZ}()
   ccall((:fmpq_harmonic_ui, :libflint), Void, (Ptr{fmpq}, Int), &(c.data), n)
   return c
end

function dedekind_sum(h::ZZ, k::ZZ)
   c = Fraction{ZZ}()
   ccall((:fmpq_dedekind_sum, :libflint), Void, (Ptr{fmpq}, Ptr{ZZ}, Ptr{ZZ}), &(c.data), &h, &k)
   return c
end

dedekind_sum(h::ZZ, k::Int) = dedekind_sum(h, ZZ(k))

dedekind_sum(h::Int, k::ZZ) = dedekind_sum(ZZ(h), k)

dedekind_sum(h::Int, k::Int) = dedekind_sum(ZZ(h), ZZ(k))

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
