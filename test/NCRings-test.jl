include("generic/MatrixAlgebra-test.jl")
include("generic/NCPoly-test.jl")

@testset "Rings.oftype..." begin
   F = GF(3)
   Fx, x = PolynomialRing(F, "x")
   z = oftype(x, 3)
   @test z == Fx(3)
   @test parent(z) === Fx
end
