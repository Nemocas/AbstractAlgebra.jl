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

#        @test LaurentPolyWrap(x) == y
#        @test LaurentPolyWrap(2x^2 + 3x + 4) == 2y^2 + 3x + 4
#        @test LaurentPolyWrap(2x^2 + 3x + 4, -1) * y == 2y^2 + 3x + 4
      end
   end

   @testset "basic manipulation" begin
      L, y = LaurentPolynomialRing(ZZ, "y")

      @test monomials_degrees(y) == 0:1
      # @test monomials_degrees(y^3) == 0:3
      @test [coeff(y, i) for i=-1:2] == [0, 0, 1, 0]

      x = y.poly
      f = LaurentPolyWrap(x, -2)
      @test monomials_degrees(f) == -2:-1
      @test [coeff(f, i) for i = -3:0] == [0, 0, 1, 0]

      f = LaurentPolyWrap(3 + 2*x^4, -3)
      @test monomials_degrees(f) == -3:1
      @test [coeff(f, i) for i = -4:2] == [0, 3, 0, 0, 0, 2, 0]
   end

   @testset "comparisons" begin
      L, y = LaurentPolynomialRing(ZZ, "y")
      x = y.poly

      @test y == y

      f = LaurentPolyWrap(x^3 + 2x^2 - 1)
      @test f == f
      @test f == LaurentPolyWrap(x^3 + 2x^2 - 1)
      @test f != x
      @test f != LaurentPolyWrap(x^3 + 2x^2 - 1, -2)
   end

   @testset "unary & binary arithmetic operations" begin
      L, y = LaurentPolynomialRing(ZZ, "y")
      x = y.poly

      @test -(-y) == y
#     @test iszero(y + (-y))
      @test y + y - y == y
      @test y*y*y + y*y - y == LaurentPolyWrap(x^3 + x^2 - x)

      f_ = rand(parent(x), 1:9, -9:9)
      f = LaurentPolyWrap(f_)
      @test f*f == LaurentPolyWrap(f_*f_)

      f = LaurentPolyWrap(f_, -3)
      @test y*y*y*y*y*y*f*f == LaurentPolyWrap(f_*f_)
   end

   @testset "printing" begin
      L, y = LaurentPolynomialRing(ZZ, "y")
      @test string(y) == "y"
      @test string(L) == "Univariate Laurent Polynomial Ring in y over Integers"
   end
end
