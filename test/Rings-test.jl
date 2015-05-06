module Test

using Base.Test, Nemo

export test_all

include("ZZ-test.jl")
include("Poly-test.jl")
include("Residue-test.jl")
include("PowerSeries-test.jl")
include("nmod_mat-test.jl")

include("Fields-test.jl")

function test_rings()
   test_zz()
   test_poly()
   test_residue()
   test_series()
   test_nmod_mat()
end

function test_all()
   test_rings()
   test_fields()
end

end # module
