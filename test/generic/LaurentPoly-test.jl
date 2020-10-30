using AbstractAlgebra: terms_degrees, LaurentPolyElem

using AbstractAlgebra.Generic: Integers, LaurentPolyWrapRing, LaurentPolyWrap,
                               trail_degree, lead_degree

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
         @test L()  == zero(L)
         c = 3 * one(R)
         @test L(c) == 3
         @test_throws Exception L('x')
         @test_throws Exception L("x")
         @test_throws Exception L(y0)

         @test isdomain_type(typeof(y))
         @test isexact_type(typeof(y))

         R, r = LaurentPolynomialRing(RDF, "r")
         @test !isexact_type(typeof(r))
      end
   end

   @testset "basic manipulation" begin
      L, y = LaurentPolynomialRing(ZZ, "y")
      x = y.poly
      Z, z = LaurentPolynomialRing(L, "z")
      T, t = LaurentPolynomialRing(L.polyring, "t")

      @test terms_degrees(y) == 0:1
      @test terms_degrees(y^3) == 0:3
      @test [coeff(y, i) for i=-1:2] == [0, 0, 1, 0]
      @test trail_degree(y^2) == 2
      @test lead_degree(y^2) == 2

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

      @test ismonomial(y)
      @test ismonomial(y^-3)
      @test !ismonomial(2y)
      @test !ismonomial(y^-1 + y)
      @test ismonomial(z^2)
      # TODO: remove Z constructor below, when ambiguities are fixed
      @test ismonomial(Z(y^-3)*z^4)
      @test ismonomial_recursive(Z(y^-3)*z^4)
      @test !ismonomial_recursive(Z(y+y^2)*z)
      @test ismonomial_recursive(x^2*t^-3)
      @test !ismonomial(x^2*t^-3)
      @test !ismonomial_recursive((x+x^2)*t)

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
      @test terms_degrees(f) == -2:-1
      @test trail_degree(f) == -1
      @test lead_degree(f) == -1
      @test [coeff(f, i) for i = -3:0] == [0, 0, 1, 0]

      @test !isone(f)
      @test !iszero(f)
      @test isone(f^0)
      @test iszero(f-f)

      @test lead(f) == 1
      @test trail(f) == 1

      f2 = f = LaurentPolyWrap(3 + 2*x^4, -3)
      @test terms_degrees(f) == -3:1
      @test trail_degree(f) == -3
      @test lead_degree(f) == 1
      @test [coeff(f, i) for i = -4:2] == [0, 3, 0, 0, 0, 2, 0]

      @test f == 3y^-3 + 2y

      @test lead(f) == 2
      @test trail(f) == 3

      @test canonical_unit(f) == 1

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

      ff = deepcopy(f)
      @test parent(f) === parent(ff)
      @test f == ff && f !== ff
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

      # with polynomials as base ring
      P, x = PolynomialRing(ZZ, "x")
      L, y = LaurentPolynomialRing(P, "y")
      @test parent(x*y) == L
      @test parent(y*x) == L
      @test 3*y == y*3 == LaurentPolyWrap(3*x, 0)

      # with Laurent polynomials as base ring
      P, x = LaurentPolynomialRing(ZZ, "x")
      L, y = LaurentPolynomialRing(P, "y")
      @test parent(x*y) == L
      @test parent(y*x) == L
      @test 3*y == y*3 == LaurentPolyWrap((3*x).poly, 0)

      # as base ring of polynomials
      L, y = LaurentPolynomialRing(ZZ, "y")
      P, x = PolynomialRing(L, "x")
      @test parent(x*y) == P
      @test parent(y*x) == P

      # Inexact field
      R, x = LaurentPolynomialRing(RealField, "x")
      for iter = 1:100
         f = rand(R, 0:10, -1:1)
         g = rand(R, 0:10, -1:1)
         h = rand(R, 0:10, -1:1)
         @test isapprox(f + (g + h), (f + g) + h)
         @test isapprox(f*g, g*f)
         @test isapprox(f*(g + h), f*g + f*h)
         @test isapprox((f - h) + (g + h), f + g)
         @test isapprox((f + g)*(f - g), f*f - g*g)
         @test isapprox(f - g, -(g - f))
      end
      p = 1.2*x
      q = nextfloat(1.2)*x
      @test p != q
      @test isapprox(p, q)
      @test !isapprox(p, q + 1.2x^2)
      @test !isapprox(p, 1.1*x)

      t = 1.2 * x^0
      r = RealField(1.2)
      @test isapprox(r, t) && isapprox(t, r)
      r = RealField(nextfloat(1.2))
      @test isapprox(r, t) && isapprox(t, r)
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
      @test_throws DomainError (y-y)^-1

      LQ, z = LaurentPolynomialRing(QQ, "z")

      @test (2z)^-1 == 1//2 * z^-1
      @test (3z^-1)^(-2) == 1//9 * z^2

      @test_throws DomainError (z + z^2)^-1
   end

   @testset "evaluate" begin
      L, y = LaurentPolynomialRing(ZZ, "y")
      p = 2y+3y^4
      @assert p.mindeg == 0
      for a = Any[-3:3; -10.0:3.3:10;]
         @test evaluate(p, a) == evaluate(p.poly, a)
      end
      q = p - y^-2
      x = y.poly
      t = -x^0 + 2x^3 + 3x^6
      for a = -10.0:3.3:10
         @test evaluate(q, a) == evaluate(t, a) * a^-2
      end
   end

   @testset "unsafe functions" begin
      L, y = LaurentPolynomialRing(ZZ, "y")

      # zero!
      p = y^-2 + 3y + y^3
      q = zero!(p)
      @test q === p
      @test iszero(p)
      p = rand(L, -10:10, -10:10)
      q = zero!(p)
      @test q === p
      @test iszero(p)
      # TODO: add a test for when p.poly is immutable

      # mul!
      p = y+1
      q = y-1
      s = y^0
      t = mul!(s, p, q)
      @test t === s == y^2-1
      p = rand(L, -10:10, -10:10)
      q = rand(L, -10:10, -10:10)
      t = p*q
      mul!(s, p, q)
      @test t == s
      # TODO: add a test for when s.poly is immutable

      # addeq!
      p = rand(L, -10:10, -10:10)
      q = rand(L, -10:10, -10:10)
      t = p + q
      s = addeq!(p, q)
      @test s === p == t
      # TODO: add a test for when p.poly is immutable

      # add!
      p = rand(L, -10:10, -10:10)
      q = rand(L, -10:10, -10:10)
      t = p + q
      s = y^0
      t = add!(s, p, q)
      @test t === s == p + q
      # TODO: add a test for when p.poly is immutable
   end

   @testset "shifting" begin
      L, y = LaurentPolynomialRing(ZZ, "y")

      p = 2y - 3y^-2
      @test shift_left(p, 0) == p
      @test shift_left(p, 3) == 2y^4 - 3y
      @test_throws DomainError shift_left(p, -rand(1:99))

      @test shift_right(p, 0) == p
      @test shift_right(p, 3) == 2y^-2 - 3y^-5
      @test_throws DomainError shift_right(p, -rand(1:99))
   end

   @testset "rand" begin
      L, y = LaurentPolynomialRing(ZZ, "y")

      m = make(L, -5:5, -10:10)
      for f in (rand(m), rand(rng, m),
                rand(L, -5:5, -10:10),
                rand(rng, L, -5:5, -10:10))

         @test f isa LaurentPolyElem{BigInt}
         @test AbstractAlgebra.degrees_range(f) ⊆ -5:5
         for i = -5:5
            @test coeff(f, i) ∈ -10:10
         end
      end
      @test reproducible(m)
      @test reproducible(L, -5:5, -10:10)
   end

   @testset "change_base_ring & map_coeffs" begin
      Z, z = LaurentPolynomialRing(ZZ, "z")
      Q, q = LaurentPolynomialRing(QQ, "q")

      fz = z^2 - z - 2z^-2

      @test change_base_ring(QQ, z) == q
      @test change_base_ring(QQ, fz) == q^2 - q - 2q^-2
      @test_broken change_base_ring(ZZ, q) == z

      @test map_coeffs(x -> x^2, fz) == z^2 + z + 4z^-2
      @test map_coeffs(one, fz) == z^2 + z + z^-2
      @test map_coeffs(x -> x+2, fz) == 3z^2 + z
      @test map_coeffs(x -> x^2, q^2 - q - 2q^-2) == q^2 + q + 4q^-2
   end

   @testset "printing" begin
      L, y = LaurentPolynomialRing(ZZ, "y")
      @test sprint(show, "text/plain", y) == "y"
      @test sprint(show, "text/plain", L) == "Univariate Laurent Polynomial Ring in y over Integers"
      p = y^1; p.mindeg = -3
      @test sprint(show, "text/plain", p) == "y^-2"
      R, z = PolynomialRing(L, "z")
      @test sprint(show, "text/plain", (y^2)*z) == "y^2*z"
      @test sprint(show, "text/plain", 3*(y^0)*z) == "3*z"
      @test sprint(show, "text/plain", -y*z + (-y*z^2)) == "-y*z^2 - y*z"
      @test sprint(show, "text/plain", -y^0*z) == "-z"
   end
end
