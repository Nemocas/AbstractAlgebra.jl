###############################################################################
#
#   pari_int.jl : functions for t_INT objects
#
###############################################################################

export PariZZ, pari, fmpz!, factor, PariIntegerRing, pari_int

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

PariZZ = PariIntegerRing()

parent(a::pari_int) = PariZZ

###############################################################################
#
#   Basic manipulation
#
###############################################################################

size(a::pari_int) = (unsafe_load(a.d, 1) & LGBITS) - 2

sign(a::pari_int) = signe(unsafe_load(a.d, 2)::Ptr{Int})

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function show(io :: IO, a :: PariIntegerRing)
   print(io, "Integer Ring")
end

show(io::IO, x::pari_int) = pari_print(io, x.d)

###############################################################################
#
#   Conversions to from BigInt
#
###############################################################################

gensize(a::fmpz) = size(a) + 2

function pari!(x::Ptr{Int}, a::fmpz, s::Int)
   siz = size(a)
   unsafe_store!(x, evaltyp(t_INT) | s, 1) 
   unsafe_store!(x, evalsigne(sign(a)) | s, 2) 
   z = BigInt(a)
   for i in 1:siz
      unsafe_store!(x, reinterpret(Int, unsafe_load(z.d, i)), i + 2)
   end
   return s
end

function pari(a::fmpz)
   s = gensize(a)
   g = pari_int(s)
   pari!(reinterpret(Ptr{Int}, g.d), a, s)
   return g
end

pari(a::Integer) = pari(fmpz(a))

function fmpz!(z::fmpz, g::Ptr{Int})
   const data_offset = 1 + div(2*sizeof(Cint), sizeof(Int))
   s = (unsafe_load(g, 1) & LGBITS) - 2
   sgn = signe(g + sizeof(Int))
   if s == 0
      ccall((:fmpz_zero, :libflint), Void, (Ptr{fmpz},), &z)
   elseif s == 1
      ccall((:fmpz_set_ui, :libflint), Void, 
            (Ptr{fmpz}, Int), &z, unsafe_load(g, 3))
   else
      a = ccall((:_fmpz_promote, :libflint), Ptr{BigInt}, (Ptr{fmpz},), &z)
      ccall((:__gmpz_realloc2, :libgmp), Void, 
            (Ptr{BigInt}, Int), a, s*BITS_IN_WORD)
      for i in 1:s
         unsafe_store!(reinterpret(Ptr{UInt}, 
                       unsafe_load(reinterpret(Ptr{Int}, a), data_offset)), 
                       reinterpret(UInt, unsafe_load(g, i + 2)), i)
      end
      unsafe_store!(reinterpret(Ptr{Cint}, a), Cint(s), 2)
   end
   if sgn < 0
      ccall((:fmpz_neg, :libflint), Void, (Ptr{fmpz}, Ptr{fmpz}), &z, &z)
   end
   return z
end

function (::FlintIntegerRing)(g::pari_int)
   z = fmpz()
   fmpz!(z, g.d)
   return z
end

###############################################################################
#
#   Factorisation
#
###############################################################################

function factor(n::pari_int)
   av = unsafe_load(avma, 1)
   f = ccall((:factorint, :libpari), Ptr{Int}, (Ptr{Int}, Int), n.d, 0)
   fac = PariFactor{pari_int}(f, PariZZ)
   unsafe_store!(avma, av, 1)
   return fac
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (ord::PariIntegerRing)(n::Ptr{Int})
   return fmpz!(fmpz(), n)
end

