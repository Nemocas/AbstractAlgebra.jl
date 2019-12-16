using AbstractAlgebra: monomials_degrees

using AbstractAlgebra.Generic: Integers, LaurentPolyWrapRing, LaurentPolyWrap

@testset "Generic.LaurentPoly" begin
   @testset "constructors" begin
      L0, y0 = LaurentPolynomialRing(zz, "y0")

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

         @test L(y) === y
         f = y^2 + y
         @test L(f) === f
         @test L(0) == zero(L)
         @test L(1) == one(L)
         @test L(2) == 2*one(L)
         @test_throws Exception L('x')
         @test_throws Exception L("x")
         @test_throws Exception L(y0)
      end
   end

   @testset "basic manipulation" begin
      L, y = LaurentPolynomialRing(ZZ, "y")
      x = y.poly

      @test monomials_degrees(y) == 0:1
      @test monomials_degrees(y^3) == 0:3
      @test [coeff(y, i) for i=-1:2] == [0, 0, 1, 0]

      @test iszero(zero(L))
      @test iszero(zero(y))
      @test !iszero(one(y))
      @test !iszero(y)

      @test isone(one(L))
      @test isone(one(y))
      @test !isone(zero(y))
      @test !isone(y)

      @test isgen(y)
      @test isgen(gen(L))
      @test !isgen(one(y))
      @test !isgen(zero(y))
      @test y == gen(L)

      @test gens(L)[1] == y
      @test length(gens(L)) == 1

      @test !isunit(zero(L))

      for e = -5:5
         @test isunit(y^e)
      end

      if base_ring(L) isa AbstractAlgebra.Field
         for e = -5:5
            @test isunit(2*y^e)
            @test isunit(3*y^(2e))
         end
      end

      @test lead(zero(y)) == 0
      @test trail(zero(y)) == 0
      @test lead(one(y)) == 1
      @test trail(one(y)) == 1
      @test lead(y) == 1
      @test trail(y) == 1

      @test hash(zero(y)) == hash(zero(y))
      @test hash(one(y)) == hash(one(y))

      f1 = f = LaurentPolyWrap(x, -2)
      @test monomials_degrees(f) == -2:-1
      @test [coeff(f, i) for i = -3:0] == [0, 0, 1, 0]

      @test !isone(f)
      @test !iszero(f)
      @test isone(f^0)
      @test iszero(f-f)

      @test lead(f) == 1
      @test trail(f) == 1

      f2 = f = LaurentPolyWrap(3 + 2*x^4, -3)
      @test monomials_degrees(f) == -3:1
      @test [coeff(f, i) for i = -4:2] == [0, 3, 0, 0, 0, 2, 0]

      @test f == 3y^-3 + 2y

      @test lead(f) == 2
      @test trail(f) == 3

      @test hash(f) != hash(3y^-3 + y)
      @test hash(f) != hash(3y^-2 + 2y)

      setcoeff!(f, -3, big(4))
      @test f == 4y^-3 + 2y

      setcoeff!(f, -3, big(0))
      @test f == 2y

      setcoeff!(f, -50, big(-2))
      @test f == -2y^-50 + 2y

      # TODO: make this work
      @test_broken iszero(setcoeff!(setcoeff!(f, 1, big(0)), -50, big(0)))

      @test !isone(f)
      @test !iszero(f)
      @test isone(f^0)
      @test iszero(f-f)

      for f in (f1, f2, LaurentPolyWrap(rand(parent(x), 0:9, -9:9), rand(-9:9)))
         @test hash(f) == hash(f)
         @test hash(f, rand(UInt)) != hash(f) # very unlikely failure
         @test hash(f-f) == hash(zero(f))
         @test hash(f^1) == hash(f)
         @test hash(f^0) == hash(one(f))
         @test hash(f*f*f) == hash(f^3)

         @test lead(f) == lead(f.poly)
         @test trail(f) == trail(f.poly)
      end
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
      @test iszero(y + (-y))
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

      @test isone(f(rand(-9:9))^0)

      @test_throws DomainError (2y)^-1
      @test_throws DomainError (3y^-1)^-2
      @test_throws DomainError (y + y^2)^-1


      LQ, z = LaurentPolynomialRing(QQ, "z")

      @test (2z)^-1 == 1//2 * z^-1
      @test (3z^-1)^(-2) == 1//9 * z^2

      @test_throws DomainError (z + z^2)^-1
   end

   @testset "change_base_ring & map_coeffs" begin
      Z, z = LaurentPolynomialRing(ZZ, "z")
      Q, q = LaurentPolynomialRing(QQ, "q")

      fz = z^2 - z - 2z^-2

      @test change_base_ring(QQ, z) == q
      @test change_base_ring(QQ, fz) == q^2 - q - 2q^-2
      @test_broken change_base_ring(ZZ, q) == z

      @test map_coeffs(x -> x^2, fz) == z^2 + z + 4z^-2
      @test map_coeffs(one, fz) == z^2 + z + 1 + z^-1 + z^-2
      @test map_coeffs(x -> x^2, q^2 - q - 2q^-2) == q^2 + q + 4q^-2
   end

   @testset "printing" begin
      L, y = LaurentPolynomialRing(ZZ, "y")
      @test string(y) == "y"
      @test string(L) == "Univariate Laurent Polynomial Ring in y over Integers"
   end
end
