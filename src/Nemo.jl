module Nemo

pkgdir = Pkg.dir("Nemo")

on_windows = @windows ? true : false
on_linux = @linux ? true : false

if on_windows
   push!(Libdl.DL_LOAD_PATH, "$pkgdir\\src\\lib")
else
    if "HOSTNAME" in ENV && ENV["HOSTNAME"] == "juliabox"
       push!(Libdl.DL_LOAD_PATH, "/usr/local/lib")
    elseif on_linux
       push!(Libdl.DL_LOAD_PATH, "$pkgdir/src/lib")
       Libdl.dlopen("$pkgdir/src/lib/libgmp")
       Libdl.dlopen("$pkgdir/src/lib/libmpfr")
       Libdl.dlopen("$pkgdir/src/lib/libflint")
    else
       push!(Libdl.DL_LOAD_PATH, "$pkgdir/src/lib")
    end
end

ccall((:pari_init, :libpari), Void, (Int, Int), 3000000000, 10000)

include("Rings.jl")

function __init__()
   println("")
   println("Welcome to Nemo version 0.2")
   println("")
   println("Nemo comes with absolutely no warranty whatsoever")
   println("")
end

end # module
