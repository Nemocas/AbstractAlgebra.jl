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
   
   # non-commutative rings
   A = MatrixAlgebra(ZZ, rand(1:9))
   a = rand(A, 1:9)

   B, _ = PolynomialRing(A, "y")
   b = rand(B, 1:9, 1:9)

   # matrices
   n = rand(1:9)
   C = MatrixSpace(ZZ, n, n)
   c = rand(C, 1:9)
   
   # commutative rings
   d = rand() * rand(-9:9)
   e = rand(big.(-9:9))
   f = rand(1:9)//rand(1:9)

   G, _ = PolynomialRing(ZZ, "x")
   g = rand(G, 1:9, 1:9)

   @testset "$T" for (x, T) in (x => string(nameof(typeof(x))) for x in (a, b, c, d, e, f, g))
      @test_throws DomainError powers(x, 0)
      @test_throws DomainError powers(x, -rand(1:100))
      @test_throws DomainError powers(x, 0)
      @test_throws DomainError powers(x, -rand(1:100))
      
      P = powers(x, 1)
      @test length(P) == 2
      @test isone(P[1])
      @test P[2] == x

      n = rand(2:9)
      Q = powers(x, n)
      @test length(Q) == n+1
      @test Q[1:2] == P
      @test Q[3] == x*x
   end

   # powers must error out on non-square matrices
   M = MatrixSpace(ZZ, 2, rand(3:9))
   @test_throws DomainError powers(rand(M, 1:9), rand(1:9))
end
