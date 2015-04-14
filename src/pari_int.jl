###########################################################################################
#
#   pari_int.jl : functions for t_INT objects
#
###########################################################################################

###########################################################################################
#
#   Constructors
#
###########################################################################################

type PariIntegerRing <: PariRing
end

PariZZ = PariIntegerRing()

type pari_int <: RingElem
   d::Ptr{Int}

   function pari_int(s::Int)
      g = new(ccall((:pari_malloc, :libpari), Ptr{Int}, (Int,), s*BITS_IN_WORD))
      finalizer(g, _pari_int_clear_fn)
      return g
   end
end

_pari_int_clear_fn(g::pari_int) = ccall((:pari_free, :libpari), Void, (Ptr{Uint},), g.d)

parent(a::pari_int) = PariZZ

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

size(a::pari_int) = (unsafe_load(a.d, 1) & LGBITS) - 2

sign(a::pari_int) = signe(unsafe_load(a.d, 2))

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show(io :: IO, a :: PariIntegerRing)
   print(io, "Integer Ring")
end

show(io::IO, x::pari_int) = pari_print(io, x.d)

###########################################################################################
#
#   Conversions to from BigInt
#
###########################################################################################

gensize(a::fmpz) = words(a) + 2

function pari!(x::Ptr{Int}, a::fmpz, s::Int)
   siz = words(a)
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

pari(a::Integer) = pari(ZZ(a))

function ZZ!(z::fmpz, g::Ptr{Int})
   const data_offset = 1 + div(2*sizeof(Cint), sizeof(Int))
   s = (unsafe_load(g, 1) & LGBITS) - 2
   sgn = signe(g + sizeof(Int))
   if s == 0
      ccall((:fmpz_zero, :libflint), Void, (Ptr{fmpz},), &z)
   elseif s == 1
      ccall((:fmpz_set_ui, :libflint), Void, (Ptr{fmpz}, Int), &z, unsafe_load(g, 3))
   else
      a = ccall((:_fmpz_promote, :libflint), Ptr{BigInt}, (Ptr{fmpz},), &z)
      ccall((:__gmpz_realloc2, :libgmp), Void, (Ptr{BigInt}, Int), a, s*BITS_IN_WORD)
      for i in 1:s
         unsafe_store!(reinterpret(Ptr{UInt}, unsafe_load(reinterpret(Ptr{Int}, a), data_offset)), 
                       reinterpret(UInt, unsafe_load(g, i + 2)), i)
      end
      unsafe_store!(reinterpret(Ptr{Cint}, a), Cint(s), 2)
   end
   if sgn < 0
      ccall((:fmpz_neg, :libflint), Void, (Ptr{fmpz}, Ptr{fmpz}), &z, &z)
   end
end

function call(::IntegerRing, g::pari_int)
   z = ZZ()
   ZZ!(z, g.d)
   return z
end

###########################################################################################
#
#   Factorisation
#
###########################################################################################

function factor(n::pari_int)
   av = unsafe_load(avma, 1)
   f = ccall((:factorint, :libpari), Ptr{Int}, (Ptr{Int}, Int), n.d, 0)
   fac = PariFactor{PariIntegerRing}(f, PariZZ)
   unsafe_store!(avma, av, 1)
   return fac
end

###########################################################################################
#
#   Parent object call overloads
#
###########################################################################################

function Base.call(ord::PariIntegerRing, n::Ptr{Int})
   return ZZ!(ZZ(), n)
end

