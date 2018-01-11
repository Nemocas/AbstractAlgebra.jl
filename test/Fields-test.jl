include("generic/Fraction-test.jl")

include("julia/Rationals-test.jl")
include("julia/Floats-test.jl")
include("julia/gfelem-test.jl")

function test_fields()
   test_gfelem()

   test_Rationals()
   test_Floats()

   test_gen_frac()
end
