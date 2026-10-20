using AbstractAlgebra: SparsePolynomialRing

@testset "Generic.SparsePoly.constructors" begin

   R, x = SparsePolynomialRing(ZZ, "x")
   S, y = SparsePolynomialRing(R, "y")

   @test SparsePolynomialRing(R, "y", cached = true)[1] === SparsePolynomialRing(R, "y", cached = true)[1]
   @test SparsePolynomialRing(R, "y", cached = true)[1] !== SparsePolynomialRing(R, "y", cached = false)[1]

   @test S isa Generic.SparsePolyRing

   T, z = SparsePolynomialRing(S, "z")

   @test T isa Generic.SparsePolyRing
end

@testset "Generic.SparsePoly.printing" begin

   R, x = SparsePolynomialRing(ZZ, "x")

   @test string(zero(R)) == "0"
   @test string(one(R)) == "1"
   @test string(x) == "x"
   @test string(5*x^5-3*x^3+2*x^2-x+1) == "5*x^5 - 3*x^3 + 2*x^2 - x + 1"

   S, y = SparsePolynomialRing(R, "y")

   @test string(zero(S)) == "0"
   @test string(one(S)) == "1"
   @test string(y) == "y"
   @test string(x+y+1) == "y + x + 1"
end

@testset "Generic.SparsePoly.adhoc_comparison" begin
   F = GF(7)
   R, x = SparsePolynomialRing(F, :x)

   @test R(F(3)) == F(3)
   @test F(3) == R(F(3))
   @test R(F(3)) == 3
   @test R() == F(0)
   @test x != F(3)
   @test x^2 != F(3)
end

@testset "Generic.SparsePoly.misc" begin

   R, x = SparsePolynomialRing(ZZ, "x")

   @test is_zero(x^15 - 1 - x^15 + 1)
   @test content(24*x + 15) == 3
end

@testset "Generic.SparsePoly.isequal" begin
   # `==` compares series coefficients to the lower precision, `isequal` not
   R, t = power_series_ring(QQ, 10, :t)
   S, x = SparsePolynomialRing(R, :x)
   f = (1 + O(t^2))*x
   g = (1 + O(t^3))*x

   @test f == g
   @test !isequal(f, g)
   @test isequal(f, deepcopy(f))
end
