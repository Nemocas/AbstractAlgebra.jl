module Nemo

pkgdir = Pkg.dir("Nemo")

on_windows = @windows ? true : false

if on_windows
   push!(DL_LOAD_PATH, "$pkgdir\\src\\lib")
else
   push!(DL_LOAD_PATH, "$pkgdir/src/lib")
end

include("Rings.jl")
include("Fields.jl")

function __init__()
   println("")
   println("Welcome to Nemo version 0.2")
   println("")
   println("Nemo comes with absolutely no warranty whatsoever")
   println("")
end

end # module
