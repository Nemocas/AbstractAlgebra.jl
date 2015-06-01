module Nemo

import Base: length, call, exp, promote_rule, zero, one, show, divrem, mod, 
             hash, factor, transpose!, rank, inv, lufact, truncate, reverse,
             isqrt, nextpow2, ndigits, invmod, isprime, sub

export Ring, Field, RingElem

export PolyElem, PowerSeriesElem

include("AbstractTypes.jl")

###############################################################################
#
#   Set up environment / load libraries
#
###############################################################################

pkgdir = Pkg.dir("Nemo")

on_windows = @windows ? true : false
on_linux = @linux ? true : false

if on_windows
   push!(Libdl.DL_LOAD_PATH, "$pkgdir\\local\\lib")
else
   try
      if "HOSTNAME" in ENV && ENV["HOSTNAME"] == "juliabox"
         push!(Libdl.DL_LOAD_PATH, "/usr/local/lib")
      elseif on_linux
         push!(Libdl.DL_LOAD_PATH, "$pkgdir/local/lib")
         Libdl.dlopen("$pkgdir/local/lib/libgmp")
         Libdl.dlopen("$pkgdir/local/lib/libmpfr")
         Libdl.dlopen("$pkgdir/local/lib/libflint")
      else
         push!(Libdl.DL_LOAD_PATH, "$pkgdir/local/lib")
      end
   catch
      push!(Libdl.DL_LOAD_PATH, "$pkgdir/local/lib")
   end
end

ccall((:pari_init, :libpari), Void, (Int, Int), 3000000000, 10000)

###############################################################################
#
#   Load Nemo Rings/Fields/etc
#
###############################################################################

include("Rings.jl")

###############################################################################
#
#   Library initialisation message
#
###############################################################################

function __init__()
   println("")
   println("Welcome to Nemo version 0.2")
   println("")
   println("Nemo comes with absolutely no warranty whatsoever")
   println("")
end

###############################################################################
#
#   Test code
#
###############################################################################

include("../test/Rings-test.jl")

end # module
