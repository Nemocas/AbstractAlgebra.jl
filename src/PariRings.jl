###########################################################################################
#
#   PariRings.jl : functions for Pari rings
#
###########################################################################################

import Base: sign

export pari, debug, gclone

abstract PariRing <: Ring

abstract PariField <: PariRing

###########################################################################################
#
#   Pari GEN constants
#
###########################################################################################

const t_INT      =  1
const t_REAL     =  2
const t_INTMOD   =  3
const t_FRAC     =  4
const t_FFELT    =  5
const t_COMPLEX  =  6
const t_PADIC    =  7
const t_QUAD     =  8
const t_POLMOD   =  9
const t_POL      = 10
const t_SER      = 11
const t_RFRAC    = 13
const t_QFR      = 15
const t_QFI      = 16
const t_VEC      = 17
const t_COL      = 18
const t_MAT      = 19
const t_LIST     = 20
const t_STR      = 21
const t_VECSMALL = 22
const t_CLOSURE  = 23
const t_ERROR    = 24
const t_INFINITY = 25

const BITS_IN_WORD = sizeof(Int)*8

const TYPnumBITS  = 7
const SIGNnumBITS = 2
if BITS_IN_WORD == 64
   const VARNnumBITS = 16
else
   const VARNnumBITS = 14
end
const LGnumBITS = (BITS_IN_WORD - 1 - TYPnumBITS)
const VALPnumBITS = (BITS_IN_WORD - SIGNnumBITS - VARNnumBITS)

const TYPSHIFT = (BITS_IN_WORD - TYPnumBITS)
const SIGNSHIFT = (BITS_IN_WORD - SIGNnumBITS)
const VARNSHIFT = VALPnumBITS

const SIGNBITS = ~((1 << SIGNSHIFT) - 1)
const TYPBITS  = ~((1 << TYPSHIFT) - 1)
const LGBITS = (1 << LGnumBITS) - 1

const CLONEBIT = 1<<LGnumBITS

###########################################################################################
#
#   Cloning
#
###########################################################################################

gclone(gen::Ptr{Int}) = ccall((:gclone, :libpari), Ptr{Int}, (Ptr{Int},), gen)

gunclone(gen::Ptr{Int}) = ccall((:gunclone, :libpari), Ptr{Int}, (Ptr{Int},), gen)

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

evaltyp(x::Int) = x << TYPSHIFT

evalsigne(x::Int) = x << SIGNSHIFT

evalvarn(x::Int) = x << VARNSHIFT

setsigne(x::Ptr{Int}, s::Int) = unsafe_store!(x, ((unsafe_load(x, 2)
                                           & (~SIGNBITS)) | evalsigne(s)), 2)

settyp(x::Ptr{Int}, s::Int) = unsafe_store!(x, ((unsafe_load(x, 1)
                                           & (~TYPBITS)) | evaltyp(s)), 1)

varn(x::Ptr{Int}) = (unsafe_load(x, 2) & VARNBITS) >> VARNSHIFT

typ(x::Ptr{Int}) = reinterpret(Int, reinterpret(Uint, unsafe_load(x, 1)) >> TYPSHIFT)

setvarn(x::Ptr{Int}, s::Int) = unsafe_store(x, (unsafe_load(x, 2) & ~VARNBITS) | evalvarn(s), 2)

signe(x::Ptr{Int}) = (unsafe_load(x) >> SIGNSHIFT)

lg(x::Ptr{Int}) = unsafe_load(x) & LGBITS

###########################################################################################
#
#   Debugging
#
###########################################################################################

debug(a::Ptr{Int}) = ccall((:dbgGEN, :libpari), Void, (Ptr{Int}, Int), a, -1)

include("pari_int.jl")

include("pari_poly.jl")

include("pari_polmod.jl")

include("pari_vec.jl")

###########################################################################################
#
#   Printing
#
###########################################################################################

function pari_print(a::Ptr{Int})
   cstr = ccall((:GENtostr, :libpari), Ptr{Uint8}, (Ptr{Int},), a)
   print(bytestring(cstr))
   ccall((:pari_free, :libpari), Void, (Ptr{Uint8},), cstr)
end
   
function pari_print(io::IO, a::Ptr{Int})
   cstr = ccall((:GENtostr, :libpari), Ptr{Uint8}, (Ptr{Int},), a)
   print(io, bytestring(cstr))
   ccall((:pari_free, :libpari), Void, (Ptr{Uint8},), cstr)
end
   
###########################################################################################
#
#   Pari globals
#
###########################################################################################

avma = cglobal((:avma, :libpari), Ptr{Ptr{Int}})

gen_0 = cglobal((:gen_0, :libpari), Ptr{Ptr{Int}})

gen_1 = cglobal((:gen_1, :libpari), Ptr{Ptr{Int}})

###########################################################################################
#
#   Factorization
#
###########################################################################################

type PariFactor{T <: PariRing}
   p::Ptr{Int}
   n::BigInt

   function PariFactor(p::Ptr{Int}, n::BigInt)
      r = new(gclone(p), n)
      finalizer(r, _PariFactor_unclone)
      return r
   end
end

_PariFactor_unclone(f::PariFactor) = gunclone(f.p)
