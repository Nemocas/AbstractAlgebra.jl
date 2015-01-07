export ZZ, IntegerRing, parent, show, fmpz

import Base.show

type IntegerRing <: Ring
end

call(::IntegerRing, a :: Int) = BigInt(a)

call(::IntegerRing, a :: Int128) = BigInt(a)

call(::IntegerRing, a :: String) = BigInt(a)

ZZ = IntegerRing()

function parent(a :: BigInt)
    return ZZ
end

function show(io :: IO, a :: IntegerRing)
   print(io, "Integer Ring")
end

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

      