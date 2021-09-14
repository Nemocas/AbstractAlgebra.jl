include("julia/Integers-test.jl")

include("generic/Poly-test.jl")
include("generic/SparsePoly-test.jl")
include("generic/LaurentPoly-test.jl")
include("generic/Residue-test.jl")
include("generic/ResidueField-test.jl")
include("generic/RelSeries-test.jl")
include("generic/AbsSeries-test.jl")
include("generic/LaurentSeries-test.jl")
include("generic/PuiseuxSeries-test.jl")
include("generic/AbsMSeries-test.jl")
include("generic/Matrix-test.jl")
include("generic/MPoly-test.jl")
include("generic/UnivPoly-test.jl")
include("algorithms/MPolyFactor-test.jl")
include("algorithms/DensePoly-test.jl")
include("algorithms/GenericFunctions-test.jl")

@testset "Generic.Rings.broadcast" begin
   F = GF(3)
   @test F(2) .* [F(1), F(2)] == [F(2), F(1)]
end

@testset "Generic.Rings.elem/parent_type" begin
   @test_throws MethodError parent_type('c')
   @test_throws MethodError parent_type(Char)
   @test parent_type(big(1)) == AbstractAlgebra.Integers{BigInt}

   @test_throws MethodError elem_type('c')
   @test_throws MethodError elem_type(Char)
end

include("julia/Integers-test.jl")
