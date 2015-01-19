###########################################################################################
#
#   QQ.jl : BigInt rationals
#
###########################################################################################

import Base: gcd

export QQ, FractionField, height, height_bits

###########################################################################################
#
#   Types and memory management
#
###########################################################################################

type RationalField <: Field
end

call(::RationalField) = Rational(0, BigInt(1))

call(::RationalField, a :: Integer) = Rational(a, BigInt(1))

call(::RationalField, a::Rational{BigInt}) = a

QQ = RationalField()

parent(a::Rational{BigInt}) = QQ

elem_type(::RationalField) = Rational{BigInt}

base_ring(a::RationalField) = ZZ

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

zero(a::RationalField) = Rational(0, BigInt(1))

one(a::RationalField) = Rational(1, BigInt(1))

isone(a::Rational{BigInt}) = a == 1

iszero(a::Rational{BigInt}) = a == 0

isunit(a::Rational{BigInt}) = a != 0

function height(a::Rational{BigInt})
   temp = fmpq_readonly(a)
   temp2 = fmpz()
   ccall((:fmpq_height, :libflint), Void, (Ptr{fmpz}, Ptr{fmpq_readonly}), &temp2, &temp)
   return BigInt(temp2)
end

function height_bits(a::Rational{BigInt})
   temp = fmpq_readonly(a)
   return ccall((:fmpq_height_bits, :libflint), Int, (Ptr{fmpq_readonly},), &temp)
end

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io::IO, a::RationalField)
   print(io, "Rational Field")
end

needs_parentheses(x::Rational{BigInt}) = false

is_negative(x::Rational{BigInt}) = x < 0

show_minus_one(::Type{Rational{BigInt}}) = false

###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit(a::Rational{BigInt}) = a

###########################################################################################
#
#   Shifting
#
###########################################################################################

function >>(a::Rational{BigInt}, b::Int)
   temp1 = fmpq_readonly(a)
   temp2 = fmpq()
   ccall((:fmpq_div_2exp, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq_readonly}, Int), &temp2, &temp1, b)
   return Rational(temp2)
end

function <<(a::Rational{BigInt}, b::Int)
   temp1 = fmpq_readonly(a)
   temp2 = fmpq()
   ccall((:fmpq_mul_2exp, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq_readonly}, Int), &temp2, &temp1, b)
   return Rational(temp2)
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

divexact(a::Rational{BigInt}, b::Rational{BigInt}) = a/b

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv(a::Rational{BigInt})
   temp1 = fmpq_readonly(a)
   temp2 = fmpq()
   ccall((:fmpq_inv, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq_readonly}), &temp2, &temp1)
   return Rational(temp2)
end

###########################################################################################
#
#   GCD
#
###########################################################################################

function gcd(a::Rational{BigInt}, b::Rational{BigInt})
   temp1 = fmpq_readonly(a)
   temp2 = fmpq_readonly(b)
   temp3 = fmpq()
   ccall((:fmpq_gcd, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq_readonly}, Ptr{fmpq_readonly}), &temp3, &temp1, &temp2)
   return Rational(temp3)
end

###########################################################################################
#
#   Conversions to/from flint fmpq
#
###########################################################################################

type fmpq_readonly
   num::Int
   den::Int

   function fmpq_readonly(x::Rational{BigInt})
      r = new()
      ccall((:fmpq_init_set_mpz_frac_readonly, :libflint), Void, 
            (Ptr{fmpq_readonly}, Ptr{BigInt}, Ptr{BigInt}), &r, &x.num, &x.den)
      return r
   end
end

type fmpq
   num::Int
   den::Int

   function fmpq()
      a = new(0, 1)
      finalizer(a, _fmpq_clear_fn)
      return a
   end
end

function _fmpq_clear_fn(z::fmpq)
   ccall((:fmpq_clear, :libflint), Void, (Ptr{fmpq},), &z)
end
   
function Rational(z::fmpq)
   r = QQ()
   ccall((:fmpq_get_mpz_frac, :libflint), Void, 
         (Ptr{BigInt}, Ptr{BigInt}, Ptr{fmpq}), &r.num, &r.den, &z)
   return r
end

###########################################################################################
#
#   RationalField constructor
#
###########################################################################################

FractionField(base::IntegerRing) = QQ 

