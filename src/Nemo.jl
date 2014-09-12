module Nemo

pkgdir = Pkg.dir("Nemo")
pwd = "$pkgdir/src"
push!(DL_LOAD_PATH, "$pwd/lib")

include("Rings.jl")
include("Fields.jl")

end # module
