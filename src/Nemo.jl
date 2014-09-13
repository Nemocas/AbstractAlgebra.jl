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

end # module
