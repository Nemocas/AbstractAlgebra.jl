module Nemo

pwd = chomp(readall(`pwd`))
push!(DL_LOAD_PATH, "$pwd/lib")

include("Rings.jl")
include("Fields.jl")

end # module
