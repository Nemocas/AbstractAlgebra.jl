using AbstractAlgebra: monomials_degrees

using AbstractAlgebra.Generic: Integers, LaurentPolyWrapRing, LaurentPolyWrap

@testset "Generic.LaurentPoly" begin
   @testset "constructors" begin
      for R in (ZZ, GF(5))
         P, x = PolynomialRing(R, "x")
         L, y = LaurentPolynomialRing(R, "y")

         @test L isa LaurentPolyWrapRing{elem_type(R)}
         @test y isa LaurentPolyWrap{elem_type(R)}

         @test parent_type(y) == typeof(L)
         @test elem_type(L) == typeof(y)

         @test parent(y) == L

         @test base_ring(L) == R
         @test base_ring(L) == base_ring(P)

         @test var(L) == :y
         @test symbols(L) == [:y]
         @test nvars(L) == 1

         @test characteristic(L) == characteristic(R)
      end
   end

   @testset "basic manipulation" begin
      L, y = LaurentPolynomialRing(ZZ, "y")

      @test monomials_degrees(y) == 0:1
      # @test monomials_degrees(y^3) == 0:3

      x = y.poly
      f = LaurentPolyWrap(x, -2)
      @test monomials_degrees(f) == -2:-1

      f = LaurentPolyWrap(1+x^4, -3)
      @test monomials_degrees(f) == -3:1
   end

   @testset "printing" begin
      L, y = LaurentPolynomialRing(ZZ, "y")
      @test string(y) == "y"
      @test string(L) == "Univariate Laurent Polynomial Ring in y over Integers"
   end
end
