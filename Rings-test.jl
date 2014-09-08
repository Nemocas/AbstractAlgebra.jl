module Test

using Base.Test, Rings

export test_all

include("ZZ-test.jl")
include("Poly-test.jl")
include("Residue-test.jl")

function test_all()
   test_zz()
   test_poly()
   test_residue()
end

end # module
