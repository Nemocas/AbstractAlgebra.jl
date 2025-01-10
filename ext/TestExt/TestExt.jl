module TestExt

using AbstractAlgebra
isdefined(Base, :get_extension) ? (using Test) : (using ..Test)



end # module
