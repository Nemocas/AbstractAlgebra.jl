module Test

using Base.Test, Nemo.Rings

export test_all

include("ZZ-test.jl")
include("Poly-test.jl")
include("Residue-test.jl")
include("PowerSeries-test.jl")

function test_all()
   test_zz()
   test_poly()
   test_residue()
   test_series()
end

end # module
