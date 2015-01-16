###########################################################################################
#
#   ZZ.jl : BigInts
#
###########################################################################################

export ZZ, IntegerRing, parent, show, fmpz, needs_parentheses, is_negative, show_minus_one,
       zero, one, isunit, iszero, isone, invmod, powmod, words

import Base: gcd

type IntegerRing <: Ring
end

call(::IntegerRing) = BigInt()

call(::IntegerRing, a::Integer) = BigInt(a)

call(::IntegerRing, a::String) = BigInt(a)

call(::IntegerRing, a::BigInt) = a

ZZ = IntegerRing()

parent(a::BigInt) = ZZ

elem_type(::IntegerRing) = BigInt

base_ring(a::IntegerRing) = None

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

zero(a::IntegerRing) = BigInt(0)

one(a::IntegerRing) = BigInt(1)

isone(a::BigInt) = a == 1

iszero(a::BigInt) = a == 0

isunit(a::BigInt) = a == 1 || a == -1

words(a::BigInt) = a == 0 ? 0 : div(ndigits(a, 2) + 8*sizeof(Int) - 1, 8*sizeof(Int))

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io::IO, a::IntegerRing)
   print(io, "Integer Ring")
end

needs_parentheses(x::BigInt) = false

is_negative(x::BigInt) = x < 0

show_minus_one(::Type{BigInt}) = false

###########################################################################################
#
#   Modular arithmetic
#
###########################################################################################

function powmod(x::BigInt, p::BigInt, m::BigInt)
    m <= 0 && throw(DomainError())
    if p < 0
       x = invmod(x, m)
       p = -p
    end
    r = ZZ()
    ccall((:__gmpz_powm, :libgmp), Void,
          (Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}), &r, &x, &p, &m)
    return r 
end

function powmod(x::BigInt, p::Int, m::BigInt)
    m <= 0 && throw(DomainError())
    if p < 0
       x = invmod(x, m)
       p = -p
    end
    r = ZZ()
    ccall((:__gmpz_powm_ui, :libgmp), Void,
          (Ptr{BigInt}, Ptr{BigInt}, Int, Ptr{BigInt}), &r, &x, p, &m)
    return r 
end

function invmod(x::BigInt, y::BigInt)
    y <= 0 && throw(DomainError())
    z = ZZ()
    if y == 1
        return zero(parent(x))
    end
    if !ccall((:__gmpz_invert, :libgmp), Bool, 
              (Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}), &z, &x, &y)
       error("Impossible inverse in invmod")
    end
    return z
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function mul!(a::BigInt, b::BigInt, c::BigInt)
   ccall((:__gmpz_mul, :libgmp), Void, (Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}), &a, &b, &c)
end

function addeq!(a::BigInt, b::BigInt)
   ccall((:__gmpz_add, :libgmp), Void, (Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}), &a, &a, &b)
end

###########################################################################################
#
#   Conversions to/from flint fmpz
#
###########################################################################################

type fmpz_readonly
   data :: Int

   function fmpz_readonly(x::BigInt)
      r = new()
      ccall((:fmpz_init_set_readonly, :libflint), Void, 
            (Ptr{fmpz_readonly}, Ptr{BigInt}), &r, &x)
      return r
   end
end

type fmpz
   data :: Int

   function fmpz()
      r = new(0)
      finalizer(r, _fmpz_clear_fn)
      return r
   end
end

function _fmpz_clear_fn(z::fmpz)
   ccall((:fmpz_clear, :libflint), Void, (Ptr{fmpz},), &z)
end

function BigInt(z::fmpz)
   r = BigInt()
   ccall((:fmpz_get_mpz, :libflint), Void, (Ptr{BigInt}, Ptr{fmpz}), &r, &z)
   return r
end

      