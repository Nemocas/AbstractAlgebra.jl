module Nemo

pkgdir = Pkg.dir("Nemo")

on_windows = @windows ? true : false

if on_windows
   push!(DL_LOAD_PATH, "$pkgdir\\src\\lib")
else
   try
      if ENV["HOSTNAME"] == "juliabox"
         push!(DL_LOAD_PATH, "/usr/local/lib")
         dlopen("/usr/local/lib/libgmp.so")
         dlopen("/usr/local/lib/libmpfr.so")
         dlopen("/usr/local/lib/libflint.so")
      else
         push!(DL_LOAD_PATH, "$pkgdir/src/lib")
      end
   catch
      push!(DL_LOAD_PATH, "$pkgdir/src/lib")
   end
end

include("Rings.jl")
include("Fields.jl")

function __init__()
   println("")
   println("Welcome to Nemo version 0.1")
   println("")
   println("Nemo comes with absolutely no warranty whatsoever")
   println("")
end

end # module
