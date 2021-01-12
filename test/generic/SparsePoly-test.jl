@testset "Generic.SparsePoly.constructors..." begin
   R, x = SparsePolynomialRing(ZZ, "x")
   S, y = SparsePolynomialRing(R, "y")

   @test typeof(S) <: Generic.SparsePolyRing

   T, z = SparsePolynomialRing(S, "z")

   @test typeof(T) <: Generic.SparsePolyRing
end

@testset "Generic.SparsePoly.printing..." begin
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
