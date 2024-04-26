# FIXME/TODO: get these conformance tests to work and pass
#function test_elem(Rx::AbstractAlgebra.Generic.SparsePolyRing)
#   R = base_ring(Rx)
#   x = gen(Rx)
#   return sum(x^(5*i) * test_elem(R) for i in 1:rand(0:6); init=zero(Rx))
#end
#
#@testset "Generic.SparsePoly.conformance" begin
#   R, x = SparsePolynomialRing(ZZ, "x")
#   test_Ring_interface(R)
#end

@testset "Generic.SparsePoly.constructors" begin
   R, x = SparsePolynomialRing(ZZ, "x")
   S, y = SparsePolynomialRing(R, "y")

   @test SparsePolynomialRing(R, "y", cached = true)[1] === SparsePolynomialRing(R, "y", cached = true)[1]
   @test SparsePolynomialRing(R, "y", cached = true)[1] !== SparsePolynomialRing(R, "y", cached = false)[1]

   @test typeof(S) <: Generic.SparsePolyRing

   T, z = SparsePolynomialRing(S, "z")

   @test typeof(T) <: Generic.SparsePolyRing
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
