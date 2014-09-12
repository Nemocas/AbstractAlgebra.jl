module Nemo

pkgdir = Pkg.dir("Nemo")

if on_windows
   push!(DL_LOAD_PATH, "$pkdir\\src\\lib")
else
   push!(DL_LOAD_PATH, "$pkdir/src/lib")
end

include("Rings.jl")
include("Fields.jl")

end # module
