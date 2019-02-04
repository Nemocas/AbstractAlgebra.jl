include("julia/Integers-test.jl")

include("generic/Poly-test.jl")
include("generic/Residue-test.jl")
include("generic/ResidueField-test.jl")
include("generic/RelSeries-test.jl")
include("generic/AbsSeries-test.jl")
include("generic/LaurentSeries-test.jl")
include("generic/PuiseuxSeries-test.jl")
include("generic/Matrix-test.jl")
include("generic/MPoly-test.jl")

function test_gen_rings_broadcast()
   print("Generic.Rings.broadcast...")

   F = GF(3)
   @test F(2) .* [F(1), F(2)] == [F(2), F(1)]
 
   println("PASS")
end

function test_rings()
   test_Integers()

   test_gen_poly()
   test_gen_res()
   test_gen_res_field()
   test_gen_abs_series()
   test_gen_rel_series()
   test_gen_laurent_series()
   test_gen_puiseux_series()
   test_gen_mat()
   test_gen_mpoly()

   test_gen_rings_broadcast()
end
