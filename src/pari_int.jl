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

gensize(a::BigInt) = words(a) + 2

function pari!(x::Ptr{Int}, a::BigInt, s::Int)
   siz = Int(a.size)
   unsafe_store!(x, evaltyp(t_INT) | s, 1) 
   unsafe_store!(x, evalsigne(sign(siz)) | s, 2) 
   for i in 1:abs(siz)
      unsafe_store!(x, reinterpret(Int, unsafe_load(a.d, i)), i + 2)
   end
   return s
end

function pari(a::BigInt)
   s = gensize(a)
   g = pari_int(s)
   pari!(reinterpret(Ptr{Int}, g.d), a, s)
   return g
end

function ZZ!(z::BigInt, g::Ptr{Int})
   s = (unsafe_load(g, 1) & LGBITS) - 2
   sgn = signe(g + sizeof(Int))
   if s > z.alloc
      ccall((:__gmpz_realloc2, :libgmp), Void, (Ptr{BigInt}, Int), &z, s*BITS_IN_WORD)
   end
   for i in 1:s
      unsafe_store!(z.d, reinterpret(UInt, unsafe_load(g, i + 2)), i)
   end
   z.size = sgn < 0 ? -s : s
end

function call(::IntegerRing, g::pari_int)
   z = BigInt()
   ZZ!(z, g.d)
   return z
end

