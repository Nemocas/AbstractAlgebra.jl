module Nemo

pkgdir = Pkg.dir("Nemo")

on_windows = @windows ? true : false

if on_windows
   push!(DL_LOAD_PATH, "$pkgdir\\src\\lib")
else
   push!(DL_LOAD_PATH, "$pkgdir/src/lib")
end

ccall((:pari_init, :libpari), Void, (Int, Int), 1000000000, 10000)

include("Rings.jl")

function __init__()
   println("")
   println("Welcome to Nemo version 0.2")
   println("")
   println("Nemo comes with absolutely no warranty whatsoever")
   println("")
end

end # module
