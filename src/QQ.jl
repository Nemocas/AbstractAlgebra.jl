###########################################################################################
#
#   QQ.jl : Rationals
#
###########################################################################################

import Base: gcd, Rational, isless

export QQ, FractionField, Rational, height, height_bits, isless

###########################################################################################
#
#   Types and memory management
#
###########################################################################################

type RationalField <: Field
end

QQ = RationalField()

type fmpq <: FieldElem
   num::Int
   den::Int

   function fmpq(a::fmpz, b::fmpz)
      z = new()
      ccall((:fmpq_init, :libflint), Void, (Ptr{fmpq},), &z)
      ccall((:fmpq_set_fmpz_frac, :libflint), Void,
            (Ptr{fmpq}, Ptr{fmpz}, Ptr{fmpz}), &z, &a, &b)
      finalizer(z, _fmpq_clear_fn)
      return z
   end

   function fmpq(a::Int, b::Int)
      z = new()
      ccall((:fmpq_init, :libflint), Void, (Ptr{fmpq},), &z)
      ccall((:fmpq_set_si, :libflint), Void,
            (Ptr{fmpq}, Int, Int), &z, a, b)
      finalizer(z, _fmpq_clear_fn)
      return z
   end
end

_fmpq_clear_fn(a::fmpq) = ccall((:fmpq_clear, :libflint), Void, (Ptr{fmpq},), &a)

parent(a::fmpq) = QQ

elem_type(::RationalField) = fmpq

base_ring(a::RationalField) = ZZ

base_ring(a::fmpq) = ZZ

###########################################################################################
#
#   Constructors
#
###########################################################################################

function //(x::fmpz, y::fmpz)
   y == 0 && throw(DivideError())
   g = gcd(x, y)
   return QQ(divexact(x, g), divexact(y, g))
end

//(x::fmpz, y::Integer) = x//ZZ(y)

//(x::Integer, y::fmpz) = ZZ(x)//y

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function hash(a::fmpq)
   h = 0x8a30b0d963237dd5
   return h $ hash(num(a)) $ hash(den(a))
end

zero(a::RationalField) = QQ(0, 1)

one(a::RationalField) = QQ(1, 1)

function num(a::fmpq)
   z = ZZ()
   ccall((:fmpq_numerator, :libflint), Void, (Ptr{fmpz}, Ptr{fmpq}), &z, &a)
   return z
end

function den(a::fmpq)
   z = ZZ()
   ccall((:fmpq_denominator, :libflint), Void, (Ptr{fmpz}, Ptr{fmpq}), &z, &a)
   return z
end

zero(a::RationalField) = QQ(0)

one(a::RationalField) = QQ(1)

isone(a::fmpq) = a == 1

iszero(a::fmpq) = a == 0

isunit(a::fmpq) = a != 0

function height(a::fmpq)
   temp = ZZ()
   ccall((:fmpq_height, :libflint), Void, (Ptr{fmpz}, Ptr{fmpq}), &temp, &a)
   return temp
end

function height_bits(a::fmpq)
   return ccall((:fmpq_height_bits, :libflint), Int, (Ptr{fmpq},), &a)
end

function deepcopy(a::fmpq)
   z = QQ()
   ccall((:fmpq_set, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &z, &a)
   return z
end

###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit(a::fmpq) = a

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io::IO, a::RationalField)
   print(io, "Rational Field")
end

function show(io::IO, a::fmpq)
   print(io, num(a))
   if den(a) != 1
      print("//", den(a))
   end
end

needs_parentheses(x::fmpq) = false

is_negative(x::fmpq) = x < 0

show_minus_one(::Type{fmpq}) = false

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::fmpq)
   z = QQ()
   ccall((:fmpq_neg, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &z, &a)
   return z
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::fmpq, b::fmpq)
   z = QQ()
   ccall((:fmpq_add, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &z, &a, &b)
   return z
end

function -(a::fmpq, b::fmpq)
   z = QQ()
   ccall((:fmpq_sub, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &z, &a, &b)
   return z
end

function *(a::fmpq, b::fmpq)
   z = QQ()
   ccall((:fmpq_mul, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &z, &a, &b)
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function +(a::fmpq, b::Int)
   z = QQ()
   ccall((:fmpq_add_si, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Int), &z, &a, b)
   return z
end

function +(a::fmpq, b::fmpz)
   z = QQ()
   ccall((:fmpq_add_fmpz, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpz}), &z, &a, &b)
   return z
end

+(a::Int, b::fmpq) = b + a

+(a::fmpz, b::fmpq) = b + a

function -(a::fmpq, b::Int)
   z = QQ()
   ccall((:fmpq_sub_si, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Int), &z, &a, b)
   return z
end

function -(a::fmpq, b::fmpz)
   z = QQ()
   ccall((:fmpq_sub_fmpz, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpz}), &z, &a, &b)
   return z
end

function *(a::fmpq, b::fmpz)
   z = QQ()
   ccall((:fmpq_mul_fmpz, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpz}), &z, &a, &b)
   return z
end

*(a::fmpz, b::fmpq) = b + a

###########################################################################################
#
#   Comparison
#
###########################################################################################

function ==(a::fmpq, b::fmpq)
   return ccall((:fmpq_equal, :libflint), Bool, (Ptr{fmpq}, Ptr{fmpq}), &a, &b)
end

function isless(a::fmpq, b::fmpq)
   return ccall((:fmpq_cmp, :libflint), Cint, (Ptr{fmpq}, Ptr{fmpq}), &a, &b) < 0
end

###########################################################################################
#
#   Ad hoc comparison
#
###########################################################################################

function ==(a::fmpq, b::Int)
   return ccall((:fmpq_equal_si, :libflint), Bool, (Ptr{fmpq}, Int), &a, b)
end

==(a::Int, b::fmpq) = b == a

function ==(a::fmpq, b::fmpz)
   return ccall((:fmpq_equal_fmpz, :libflint), Bool, (Ptr{fmpq}, Ptr{fmpz}), &a, &b)
end

==(a::fmpz, b::fmpq) = b == a

function isless(a::fmpq, b::Int)
   z = QQ(b)
   return ccall((:fmpq_cmp, :libflint), Cint, (Ptr{fmpq}, Ptr{fmpq}), &a, &z) < 0
end

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^(a::fmpq, b::Int)
   b < 0 && throw(DomainError())
   temp = QQ()
   ccall((:fmpq_pow_si, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Int), &temp, &a, b)
   return temp
end

###########################################################################################
#
#   Shifting
#
###########################################################################################

function >>(a::fmpq, b::Int)
   temp = QQ()
   ccall((:fmpq_div_2exp, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Int), &temp, &a, b)
   return temp
end

function <<(a::fmpq, b::Int)
   temp = QQ()
   ccall((:fmpq_mul_2exp, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Int), &temp, &a, b)
   return temp2
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact(a::fmpq, b::fmpq)
   z = QQ()
   ccall((:fmpq_div, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &z, &a, &b)
   return z
end

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv(a::fmpq)
   z = QQ()
   ccall((:fmpq_inv, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &z, &a)
   return z
end

###########################################################################################
#
#   GCD
#
###########################################################################################

function gcd(a::fmpq, b::fmpq)
   z = QQ()
   ccall((:fmpq_gcd, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &z, &a, &b)
   return z
end

###############################################################################
#
#   Unsafe operators and functions
#
###############################################################################

function mul!(c::fmpq, a::fmpq, b::fmpq)
   ccall((:fmpq_mul, :libflint), Void,
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &c, &a, &b)
end

function addeq!(c::fmpq, a::fmpq)
   ccall((:fmpq_add, :libflint), Void,
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &c, &c, &a)
end

###########################################################################################
#
#   Conversions to/from flint fmpq
#
###########################################################################################
   
function Rational(z::fmpq)
   r = Rational{BigInt}(0)
   ccall((:fmpq_get_mpz_frac, :libflint), Void, 
         (Ptr{BigInt}, Ptr{BigInt}, Ptr{fmpq}), &r.num, &r.den, &z)
   return r
end

function Rational(z::fmpz)
   return Rational{BigInt}(BigInt(z))
end

###########################################################################################
#
#   Parent object call overloads
#
###########################################################################################

call(a::RationalField) = fmpq(ZZ(0), ZZ(1))

call(a::RationalField, b::Rational{BigInt}) = fmpq(ZZ(b.num), ZZ(b.den)) 

call(a::RationalField, b::Integer) = fmpq(ZZ(b), ZZ(1))

call(a::RationalField, b::Int, c::Int) = fmpq(b, c)

call(a::RationalField, b::fmpz) = fmpq(b, ZZ(1))

call(a::RationalField, b::Integer, c::Integer) = fmpq(ZZ(b), ZZ(c))

call(a::RationalField, b::fmpz, c::Integer) = fmpq(b, ZZ(c))

call(a::RationalField, b::Integer, c::fmpz) = fmpq(ZZ(b), c)

call(a::RationalField, b::fmpz, c::fmpz) = fmpq(b, c)

call(a::RationalField, b::fmpq) = b

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

convert(::Type{fmpq}, a::Integer) = QQ(a)

convert(::Type{fmpq}, a::fmpz) = QQ(a)

Base.promote_rule{T <: Integer}(::Type{fmpq}, ::Type{T}) = fmpq

Base.promote_rule(::Type{fmpq}, ::Type{fmpz}) = fmpq

convert(::Type{Rational{BigInt}}, a::fmpq) = Rational(a)

###########################################################################################
#
#   RationalField constructor
#
###########################################################################################

FractionField(base::IntegerRing) = QQ 

