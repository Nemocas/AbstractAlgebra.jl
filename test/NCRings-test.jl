include("generic/MatrixAlgebra-test.jl")
include("generic/NCPoly-test.jl")

@testset "NCRings.oftype..." begin
   F = GF(3)
   Fx, x = PolynomialRing(F, "x")
   z = oftype(x, 3)
   @test z == Fx(3)
   @test parent(z) === Fx
end

@testset "NCRings.powers..." begin
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")
   @test_throws DomainError powers(y, 0)
   @test_throws DomainError powers(y, -rand(1:100))
end
