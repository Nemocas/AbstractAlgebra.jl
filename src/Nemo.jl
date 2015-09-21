VERSION >= v"0.4.0-dev+6521" && __precompile__()

module Nemo

import Base: abs, asin, asinh, atan, atanh, base, bin, call, checkbounds,
             convert, cmp, cos, cosh, dec, deepcopy, den, deserialize, div,
             divrem, exp, factor, gcd, gcdx, getindex, hash, hex, intersect,
             inv, invmod, isequal, isless, isprime, isqrt, lcm, length, log,
             lufact, mod, ndigits, nextpow2, norm, nullspace, num, oct, one, 
             parent, parseint, precision, prevpow2, promote_rule, rank, 
             Rational, rem, reverse, serialize, setindex!, show, sign, sin,
             sinh, size, sqrt, string, sub, tan, tanh, trace, trailing_zeros,
             transpose, transpose!, truncate, var, zero, +, -, *, ==, ^, &, |,
             $, <<, >>, ~, <=, >=, <, >, hcat, vcat, //

export Collection, Ring, Field, CollectionElem, RingElem, FieldElem, Pari,
       Flint, Antic, Generic

export PolyElem, SeriesElem, ResidueElem, FractionElem, MatElem,
       NumberFieldElem

export ZZ, QQ, PadicField, FiniteField, NumberField, CyclotomicField,
       MaximalRealSubfield, MaximalOrder, Ideal

export create_accessors, get_handle, package_handle, allocatemem

export flint_cleanup, flint_set_num_threads

include("AbstractTypes.jl")

###############################################################################
#
#   Set up environment / load libraries
#
###############################################################################

const pkgdir = Pkg.dir("Nemo")
const libdir = Pkg.dir("Nemo", "local", "lib")
const libgmp = Pkg.dir("Nemo", "local", "lib", "libgmp")
const libmpfr = Pkg.dir("Nemo", "local", "lib", "libmpfr")
const libflint = Pkg.dir("Nemo", "local", "lib", "libflint")
const libpari = Pkg.dir("Nemo", "local", "lib", "libpari")
const libarb = Pkg.dir("Nemo", "local", "lib", "libarb")
  
function allocatemem(bytes::Int)
   newsize = pari(fmpz(bytes)).d
   ccall((:gp_allocatemem, :libpari), Void, (Ptr{Int},), newsize)
end

function pari_sigint_handler()
   error("User interrupt")
   return
end

function __init__()

   on_windows = @windows ? true : false
   on_linux = @linux ? true : false

   if "HOSTNAME" in keys(ENV) && ENV["HOSTNAME"] == "juliabox"
       push!(Libdl.DL_LOAD_PATH, "/usr/local/lib")
   elseif on_linux
       push!(Libdl.DL_LOAD_PATH, libdir)
       Libdl.dlopen(libgmp)
       Libdl.dlopen(libmpfr)
       Libdl.dlopen(libflint)
       Libdl.dlopen(libpari)
       Libdl.dlopen(libarb)
   else
      push!(Libdl.DL_LOAD_PATH, libdir)
   end
 
   ccall((:pari_init, libpari), Void, (Int, Int), 300000000, 10000)
  
   global avma = cglobal((:avma, libpari), Ptr{Int})

   global gen_0 = cglobal((:gen_0, libpari), Ptr{Int})

   global gen_1 = cglobal((:gen_1, libpari), Ptr{Int})

   global pari_sigint = cglobal((:cb_pari_sigint, libpari), Ptr{Void})

   unsafe_store!(pari_sigint, cfunction(pari_sigint_handler, Void, ()), 1)

   println("")
   println("Welcome to Nemo version 0.3")
   println("")
   println("Nemo comes with absolutely no warranty whatsoever")
   println("")
end

function flint_set_num_threads(a::Int)
   ccall((:flint_set_num_threads, :libflint), Void, (Int,), a)
end

function flint_cleanup()
   ccall((:flint_cleanup, :libflint), Void, ())
end

###############################################################################
#
#   Load Nemo Rings/Fields/etc
#
###############################################################################

include("generic/GenericTypes.jl")

include("flint/FlintTypes.jl")

include("antic/AnticTypes.jl")

include("pari/PariTypes.jl")

include("Rings.jl")

###########################################################
#
#   Package handle creation
#
###########################################################

const package_handle = [1]

function get_handle()
   package_handle[1] += 1
   return package_handle[1] - 1
end

###############################################################################
#
#   Auxilliary data accessors
#
###############################################################################

function create_accessors(T, S, handle)
   accessor_name = gensym()
   @eval begin
      function $(symbol(:get, accessor_name))(a::$T)
         return a.auxilliary_data[$handle]::$S
      end,
      function $(symbol(:set, accessor_name))(a::$T, b::$S)
         if $handle > length(a.auxilliary_data)
            resize(a.auxilliary_data, $handle)
         end
         a.auxilliary_data[$handle] = b
      end
   end
end

###############################################################################
#
#   Set domain for ZZ, QQ, PadicField, FiniteField to Flint
#
###############################################################################

ZZ = FlintZZ
QQ = FlintQQ
PadicField = FlintPadicField
FiniteField = FlintFiniteField

###############################################################################
#
#   Set domain for NumberField to Antic
#
###############################################################################

NumberField = AnticNumberField
CyclotomicField = AnticCyclotomicField
MaximalRealSubfield = AnticMaximalRealSubfield

###############################################################################
#
#   Set domain for MaximalOrder and Ideal to Pari
#
###############################################################################

MaximalOrder = PariMaximalOrder
Ideal = PariIdeal

###############################################################################
#
#   Test code
#
###############################################################################

include("../test/Rings-test.jl")

end # module
