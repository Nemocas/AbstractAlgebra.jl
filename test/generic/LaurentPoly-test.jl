using AbstractAlgebra: monomials_degrees

using AbstractAlgebra.Generic: Integers, LaurentPolyWrapRing, LaurentPolyWrap

@testset "Generic.LaurentPoly" begin
   @testset "constructors" begin
      for R in (ZZ, GF(5))
         P, _ = PolynomialRing(R, "x0")
         L, y = LaurentPolynomialRing(R, "y")
         x = y.poly

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

         @test LaurentPolyWrap(x) == y
         @test LaurentPolyWrap(2x^2 + 3x + 4) == 2y^2 + 3y + 4
         @test LaurentPolyWrap(2x^2 + 3x + 4, -1) * y == 2y^2 + 3y + 4
      end
   end

   @testset "basic manipulation" begin
      L, y = LaurentPolynomialRing(ZZ, "y")

      @test monomials_degrees(y) == 0:1
      @test monomials_degrees(y^3) == 0:3
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

   @testset "unary & binary & adhoc arithmetic operations" begin
      L, y = LaurentPolynomialRing(ZZ, "y")
      x = y.poly

      @test -(-y) == y
#     @test iszero(y + (-y))
      @test y + y - y == y
      @test 2y*y*y + 3y*y - 5y + 8 == LaurentPolyWrap(2x^3 + 3x^2 - 5x + 8)

      c = rand(-9:9)
      for i = -9:9
         @test coeff(c*y, i) == (i == 1 ? c : 0)
         @test coeff(y*c, i) == (i == 1 ? c : 0)
      end

      fx = rand(parent(x), 1:9, -9:9)
      f = LaurentPolyWrap(fx)
      @test f*f == LaurentPolyWrap(fx*fx)

      f = LaurentPolyWrap(fx, -3)
      @test y*y*y*y*y*y*f*f == LaurentPolyWrap(fx*fx)
      @test y*y*y*y*y*y*f*f == y^6 * f^2
   end

   @testset "powering" begin
      L, y = LaurentPolynomialRing(ZZ, "y")
      x = y.poly

      @test 2y^-2 + 3y^-1 + 4y^0 + 5y + 6y^2 == LaurentPolyWrap(2 + 3x + 4x^2 + 5x^3 + 6x^4, -2)

      fx = rand(parent(x), 1:9, -9:9)
      f = d -> LaurentPolyWrap(fx, d)

      for e in rand(0:9, 2)
         fxe = LaurentPolyWrap(fx^e)
         for i = -3:3
            @test y^(-i*e) * f(i)^e == fxe
         end
      end

      @test_throws DomainError (2y)^-1
      @test_throws DomainError (3y^-1)^-2
      @test_throws DomainError (y + y^2)^-1


      LQ, z = LaurentPolynomialRing(QQ, "z")

      @test (2z)^-1 == 1//2 * z^-1
      @test (3z^-1)^(-2) == 1//9 * z^2

      @test_throws DomainError (z + z^2)^-1
   end

   @testset "printing" begin
      L, y = LaurentPolynomialRing(ZZ, "y")
      @test string(y) == "y"
      @test string(L) == "Univariate Laurent Polynomial Ring in y over Integers"
   end
end
