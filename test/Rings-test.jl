module Test

using Base.Test, Nemo

export test_all

include("flint/fmpz-test.jl")
include("flint/fmpz_poly-test.jl")
include("flint/fmpz_mod_poly-test.jl")
include("flint/nmod_poly-test.jl")

include("generic/Poly-test.jl")
include("generic/Residue-test.jl")
include("generic/PowerSeries-test.jl")

include("nmod_mat-test.jl")
include("fmpz_series-test.jl")
include("Fields-test.jl")

function test_rings()
   test_fmpz()
   test_fmpz_poly()
   test_fmpz_mod_poly()
   test_nmod_poly()

   test_poly()
   test_residue()
   test_series()

   test_nmod_mat()
   test_fmpz_series()
end

function test_all()
   test_rings()
   test_fields()
end

end # module
