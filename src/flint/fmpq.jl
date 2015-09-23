###############################################################################
#
#   fmpq.jl : Flint rationals
#
###############################################################################

export fmpq, FlintQQ, FractionField, Rational, FlintRationalField, height,
       height_bits, isless, reconstruct, next_minimal, next_signed_minimal,
       next_calkin_wilf, next_signed_calkin_wilf, dedekind_sum, harmonic,
       bernoulli, bernoulli_cache

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

fmpq(a::Rational{BigInt}) = fmpq(fmpz(a.num), fmpz(a.den))

fmpq(a::Integer) = fmpq(fmpz(a), fmpz(1))

fmpq(a::fmpz) = fmpq(a, fmpz(1))

fmpq(a::Integer, b::Integer) = fmpq(fmpz(a), fmpz(b))

fmpq(a::fmpz, b::Integer) = fmpq(a, fmpz(b))

fmpq(a::Integer, b::fmpz) = fmpq(fmpz(a), b)

parent(a::fmpq) = FlintQQ

elem_type(::FlintRationalField) = fmpq

base_ring(a::FlintRationalField) = FlintZZ

base_ring(a::fmpq) = FlintZZ

###############################################################################
#
#   Constructors
#
###############################################################################

function //(x::fmpz, y::fmpz)
   y == 0 && throw(DivideError())
   g = gcd(x, y)
   return fmpq(divexact(x, g), divexact(y, g))
end

//(x::fmpz, y::Integer) = x//fmpz(y)

//(x::Integer, y::fmpz) = fmpz(x)//y

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function hash(a::fmpq)
   h = 0x8a30b0d963237dd5
   return h $ hash(num(a)) $ hash(den(a))
end

zero(a::FlintRationalField) = fmpq(0, 1)

one(a::FlintRationalField) = fmpq(1, 1)

function num(a::fmpq)
   z = fmpz()
   ccall((:fmpq_numerator, :libflint), Void, (Ptr{fmpz}, Ptr{fmpq}), &z, &a)
   return z
end

function den(a::fmpq)
   z = fmpz()
   ccall((:fmpq_denominator, :libflint), Void, (Ptr{fmpz}, Ptr{fmpq}), &z, &a)
   return z
end

function abs(a::fmpq)
   z = fmpq()
   ccall((:fmpq_abs, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &z, &a)
   return z
end

zero(a::FlintRationalField) = fmpq(0)

one(a::FlintRationalField) = fmpq(1)

isone(a::fmpq) = a == 1

iszero(a::fmpq) = a == 0

isunit(a::fmpq) = a != 0

function height(a::fmpq)
   temp = fmpz()
   ccall((:fmpq_height, :libflint), Void, (Ptr{fmpz}, Ptr{fmpq}), &temp, &a)
   return temp
end

function height_bits(a::fmpq)
   return ccall((:fmpq_height_bits, :libflint), Int, (Ptr{fmpq},), &a)
end

function deepcopy(a::fmpq)
   z = fmpq()
   ccall((:fmpq_set, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &z, &a)
   return z
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::fmpq) = a

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, a::FlintRationalField)
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
   z = fmpq()
   ccall((:fmpq_neg, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &z, &a)
   return z
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::fmpq, b::fmpq)
   z = fmpq()
   ccall((:fmpq_add, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &z, &a, &b)
   return z
end

function -(a::fmpq, b::fmpq)
   z = fmpq()
   ccall((:fmpq_sub, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &z, &a, &b)
   return z
end

function *(a::fmpq, b::fmpq)
   z = fmpq()
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
   z = fmpq()
   ccall((:fmpq_add_si, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Int), &z, &a, b)
   return z
end

function +(a::fmpq, b::fmpz)
   z = fmpq()
   ccall((:fmpq_add_fmpz, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpz}), &z, &a, &b)
   return z
end

+(a::Int, b::fmpq) = b + a

+(a::fmpz, b::fmpq) = b + a

function -(a::fmpq, b::Int)
   z = fmpq()
   ccall((:fmpq_sub_si, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Int), &z, &a, b)
   return z
end

function -(a::fmpq, b::fmpz)
   z = fmpq()
   ccall((:fmpq_sub_fmpz, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpz}), &z, &a, &b)
   return z
end

function *(a::fmpq, b::fmpz)
   z = fmpq()
   ccall((:fmpq_mul_fmpz, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpz}), &z, &a, &b)
   return z
end

*(a::fmpz, b::fmpq) = b*a

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::fmpq, b::fmpq)
   return ccall((:fmpq_equal, :libflint), Bool, 
                (Ptr{fmpq}, Ptr{fmpq}), &a, &b)
end

function isless(a::fmpq, b::fmpq)
   return ccall((:fmpq_cmp, :libflint), Cint, 
                (Ptr{fmpq}, Ptr{fmpq}), &a, &b) < 0
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

function ==(a::fmpq, b::Int)
   return ccall((:fmpq_equal_si, :libflint), Bool, (Ptr{fmpq}, Int), &a, b)
end

==(a::Int, b::fmpq) = b == a

function ==(a::fmpq, b::fmpz)
   return ccall((:fmpq_equal_fmpz, :libflint), Bool, 
                (Ptr{fmpq}, Ptr{fmpz}), &a, &b)
end

==(a::fmpz, b::fmpq) = b == a

function isless(a::fmpq, b::Integer)
   z = fmpq(b)
   return ccall((:fmpq_cmp, :libflint), Cint,
                (Ptr{fmpq}, Ptr{fmpq}), &a, &z) < 0
end

function isless(a::Integer, b::fmpq)
   z = fmpq(a)
   return ccall((:fmpq_cmp, :libflint), Cint,
                (Ptr{fmpq}, Ptr{fmpq}), &z, &b) < 0
end

function isless(a::fmpq, b::fmpz)
   z = fmpq(b)
   return ccall((:fmpq_cmp, :libflint), Cint,
                (Ptr{fmpq}, Ptr{fmpq}), &a, &z) < 0
end

function isless(a::fmpz, b::fmpq)
   z = fmpq(a)
   return ccall((:fmpq_cmp, :libflint), Cint,
                (Ptr{fmpq}, Ptr{fmpq}), &z, &b) < 0
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::fmpq, b::Int)
   temp = fmpq()
   ccall((:fmpq_pow_si, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Int), &temp, &a, b)
   return temp
end

###############################################################################
#
#   Shifting
#
###############################################################################

function >>(a::fmpq, b::Int)
   z = fmpq()
   ccall((:fmpq_div_2exp, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Int), &z, &a, b)
   return z
end

function <<(a::fmpq, b::Int)
   z = fmpq()
   ccall((:fmpq_mul_2exp, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Int), &z, &a, b)
   return z
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(a::fmpq)
   z = fmpq()
   ccall((:fmpq_inv, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &z, &a)
   return z
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::fmpq, b::fmpq)
   z = fmpq()
   ccall((:fmpq_div, :libflint), Void,
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &z, &a, &b)
   return z
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(a::fmpq, b::fmpz)
   z = fmpq()
   ccall((:fmpq_div_fmpz, :libflint), Void,
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpz}), &z, &a, &b)
   return z
end

divexact(a::fmpz, b::fmpq) = inv(b)*a

divexact(a::fmpq, b::Integer) = divexact(a, fmpz(b))

divexact(a::Integer, b::fmpq) = inv(b)*a

###############################################################################
#
#   Modular arithmetic
#
###############################################################################

function mod(a::fmpq, b::fmpz)
   z = fmpz()
   ccall((:fmpq_mod_fmpz, :libflint), Void,
         (Ptr{fmpz}, Ptr{fmpq}, Ptr{fmpz}), &z, &a, &b)
   return z
end

mod(a::fmpq, b::Integer) = mod(a, fmpz(b))

###############################################################################
#
#   GCD
#
###############################################################################

function gcd(a::fmpq, b::fmpq)
   z = fmpq()
   ccall((:fmpq_gcd, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}, Ptr{fmpq}), &z, &a, &b)
   return z
end

###############################################################################
#
#   Rational reconstruction
#
###############################################################################

function reconstruct(a::fmpz, b::fmpz)
   c = fmpq()
   if !Bool(ccall((:fmpq_reconstruct_fmpz, :libflint), Cint, 
                  (Ptr{fmpq}, Ptr{fmpz}, Ptr{fmpz}), &c, &a, &b))
      error("Impossible rational reconstruction")
   end
   return c
end

reconstruct(a::fmpz, b::Integer) =  reconstruct(a, fmpz(b))

reconstruct(a::Integer, b::fmpz) =  reconstruct(fmpz(a), b)

reconstruct(a::Integer, b::Integer) =  reconstruct(fmpz(a), fmpz(b))

###############################################################################
#
#   Rational enumeration
#
###############################################################################

function next_minimal(a::fmpq)
   a < 0 && throw(DomainError())
   c = fmpq()
   ccall((:fmpq_next_minimal, :libflint), Void, (Ptr{fmpq}, Ptr{fmpq}), &c, &a)
   return c
end

function next_signed_minimal(a::fmpq)
   c = fmpq()
   ccall((:fmpq_next_signed_minimal, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}), &c, &a)
   return c
end

function next_calkin_wilf(a::fmpq)
   a < 0 && throw(DomainError())
   c = fmpq()
   ccall((:fmpq_next_calkin_wilf, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}), &c, &a)
   return c
end

function next_signed_calkin_wilf(a::fmpq)
   c = fmpq()
   ccall((:fmpq_next_signed_calkin_wilf, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpq}), &c, &a)
   return c
end

###############################################################################
#
#   Special functions
#
###############################################################################

function harmonic(n::Int)
   n < 0 && throw(DomainError())
   c = fmpq()
   ccall((:fmpq_harmonic_ui, :libflint), Void, (Ptr{fmpq}, Int), &c, n)
   return c
end

function bernoulli(n::Int)
   n < 0 && throw(DomainError())
   c = fmpq()
   ccall((:bernoulli_fmpq_ui, :libarb), Void, (Ptr{fmpq}, Int), &c, n)
   return c
end

function bernoulli_cache(n::Int)
   n = n + 1
   n < 0 && throw(DomainError())
   ccall((:bernoulli_cache_compute, :libarb), Void, (Int,), n)
end

function dedekind_sum(h::fmpz, k::fmpz)
   c = fmpq()
   ccall((:fmpq_dedekind_sum, :libflint), Void, 
         (Ptr{fmpq}, Ptr{fmpz}, Ptr{fmpz}), &c, &h, &k)
   return c
end

dedekind_sum(h::fmpz, k::Integer) = dedekind_sum(h, fmpz(k))

dedekind_sum(h::Integer, k::fmpz) = dedekind_sum(fmpz(h), k)

dedekind_sum(h::Integer, k::Integer) = dedekind_sum(fmpz(h), fmpz(k))

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

###############################################################################
#
#   Conversions to/from flint Julia rationals
#
###############################################################################
   
function Rational(z::fmpq)
   r = Rational{BigInt}(0)
   ccall((:fmpq_get_mpz_frac, :libflint), Void, 
         (Ptr{BigInt}, Ptr{BigInt}, Ptr{fmpq}), &r.num, &r.den, &z)
   return r
end

function Rational(z::fmpz)
   return Rational{BigInt}(BigInt(z))
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

call(a::FlintRationalField) = fmpq(fmpz(0), fmpz(1))

call(a::FlintRationalField, b::Rational{BigInt}) = fmpq(b) 

call(a::FlintRationalField, b::Integer) = fmpq(b)

call(a::FlintRationalField, b::Int, c::Int) = fmpq(b, c)

call(a::FlintRationalField, b::fmpz) = fmpq(b)

call(a::FlintRationalField, b::Integer, c::Integer) = fmpq(b, c)

call(a::FlintRationalField, b::fmpz, c::Integer) = fmpq(b, c)

call(a::FlintRationalField, b::Integer, c::fmpz) = fmpq(b, c)

call(a::FlintRationalField, b::fmpz, c::fmpz) = fmpq(b, c)

call(a::FlintRationalField, b::fmpq) = b

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

convert(::Type{fmpq}, a::Integer) = fmpq(a)

convert(::Type{fmpq}, a::fmpz) = fmpq(a)

Base.promote_rule{T <: Integer}(::Type{fmpq}, ::Type{T}) = fmpq

Base.promote_rule(::Type{fmpq}, ::Type{fmpz}) = fmpq

convert(::Type{Rational{BigInt}}, a::fmpq) = Rational(a)

###############################################################################
#
#   FractionField constructor
#
###############################################################################

FractionField(base::FlintIntegerRing) = FlintQQ

