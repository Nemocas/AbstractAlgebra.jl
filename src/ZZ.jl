export ZZ, IntegerRing, parent, show, fmpz, needs_parentheses, is_negative, show_minus_one,
       zero, one

type IntegerRing <: Ring
end

call(::IntegerRing, a :: Integer) = BigInt(a)

call(::IntegerRing, a :: String) = BigInt(a)

ZZ = IntegerRing()

function parent(a :: BigInt)
    return ZZ
end

elem_type(::IntegerRing) = BigInt

base(a::IntegerRing) = None

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

zero(a::IntegerRing) = BigInt(0)

one(a::IntegerRing) = BigInt(1)

isone(a::BigInt) = a == 1

iszero(a::BigInt) = a == 0

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io :: IO, a :: IntegerRing)
   print(io, "Integer Ring")
end

needs_parentheses(x::BigInt) = false

is_negative(x::BigInt) = x < 0

show_minus_one(::Type{BigInt}) = false

###########################################################################################
#
#   Conversions to/from flint fmpz
#
###########################################################################################

type fmpz
   data :: Int

   fmpz() = new()
end

function init(z::fmpz)
   z.data = 0
end

function clear(z::fmpz)
   ccall((:fmpz_clear, :libflint), Void, (Ptr{fmpz},), &z)
end

function fmpz_readonly(x::BigInt)
   temp = fmpz()
   ccall((:fmpz_init_set_readonly, :libflint), Void, (Ptr{fmpz}, Ptr{BigInt}), &temp, &x)
   return temp
end
   
function BigInt(z::fmpz)
   r = BigInt()
   ccall((:fmpz_get_mpz, :libflint), Void, (Ptr{BigInt}, Ptr{fmpz}), &r, &z)
   return r
end

      