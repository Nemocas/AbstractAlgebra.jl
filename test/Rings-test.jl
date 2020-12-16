include("julia/Integers-test.jl")

include("generic/Poly-test.jl")
include("generic/LaurentPoly-test.jl")
include("generic/Residue-test.jl")
include("generic/ResidueField-test.jl")
include("generic/RelSeries-test.jl")
include("generic/AbsSeries-test.jl")
include("generic/LaurentSeries-test.jl")
include("generic/PuiseuxSeries-test.jl")
include("generic/Matrix-test.jl")
include("generic/MPoly-test.jl")
include("generic/MPolyFactor-test.jl")

@testset "Generic.Rings.broadcast..." begin
   F = GF(3)
   @test F(2) .* [F(1), F(2)] == [F(2), F(1)]
end

@testset "Generic.Rings.elem/parent_type" begin
   for (R, el, par) in [(GF(3), AbstractAlgebra.GFElem{Int}, AbstractAlgebra.GFField{Int}),
                        (ZZ["x"][1], Generic.Poly{BigInt}, Generic.PolyRing{BigInt})]
      x = R(2)
      @test elem_type(R) == el
      @test elem_type(typeof(R)) == el
      @test parent_type(x) == par
      @test parent_type(typeof(x)) == par
   end

   @test_throws MethodError parent_type('c')
   @test_throws MethodError parent_type(Char)
   @test parent_type(big(1)) == AbstractAlgebra.Integers{BigInt}

   @test_throws MethodError elem_type('c')
   @test_throws MethodError elem_type(Char)
end
