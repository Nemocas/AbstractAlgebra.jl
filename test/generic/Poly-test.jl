# In general we want to test over:
#    1) Exact rings, e.g. Z
#    2) Exact fields, e.g. Q or GFp
#    3) Inexact rings, e.g. polynomials over Julia RealField, or power series
#       over Z
#    4) Inexact fields, e.g. Julia RealField
#    5) A field of char p > 0, e.g. GF(p)
#    6) A ring of char p > 0, e.g. Z/pZ
#    7) Commutative ring, not an integral domain, e.g. Z/nZ or Z[x]/(f)
#       with reducible f
# In some cases, we may also wish to test over:
#    8) Polynomial rings, e.g. to test interpolation strategies
#    9) Fraction fields, such as Q, e.g. to test fraction free algorithms,
#       quasidivision, etc.
#   10) Generic towers, e.g. to test ad hoc functions
# Note: only useful to distinguish rings and fields for 1/2, 3/4, 5/6 if the
# algos differ, and 7 can often stand in for 5/6 if the algorithm supports it.

@testset "Generic.Poly.types" begin
   @test dense_poly_type(BigInt) == Generic.Poly{BigInt}
   @test dense_poly_type(Rational{BigInt}) == Generic.Poly{Rational{BigInt}}
end

@testset "Generic.Poly.constructors" begin
   S1 = PolyRing(ZZ)
   S2 = PolyRing(ZZ)

   @test S1 !== S2
   @test isa(S1, Generic.PolyRing)

   R, x = ZZ["x"]
   S1 = R["y"]
   S2 = ZZ["x"]["y"]

   @test PolynomialRing(R, "y", cached = true)[1] === PolynomialRing(R, "y", cached = true)[1]
   @test PolynomialRing(R, "y", cached = true)[1] !== PolynomialRing(R, "y", cached = false)[1]

   for (S, y) in (S1, S2)
      @test base_ring(S) === R
      @test coefficient_ring(S) === R

      @test elem_type(S) == Generic.Poly{elem_type(R)}
      @test elem_type(Generic.PolyRing{elem_type(R)}) == Generic.Poly{elem_type(R)}
      @test parent_type(Generic.Poly{elem_type(R)}) == Generic.PolyRing{elem_type(R)}

      @test typeof(R) <: AbstractAlgebra.Ring
      @test typeof(S) <: Generic.PolyRing

      @test isa(y, PolyElem)
   end

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   @test typeof(S) <: Generic.PolyRing

   @test isa(y, PolyElem)

   T, z = PolynomialRing(S, "z")

   @test typeof(T) <: Generic.PolyRing

   @test isa(z, PolyElem)

   f = x^2 + y^3 + z + 1

   @test isa(f, PolyElem)

   g = S(2)

   @test isa(g, PolyElem)

   h = S(x^2 + 2x + 1)

   @test isa(h, PolyElem)

   j = T(x + 2)

   @test isa(j, PolyElem)

   k = S([x, x + 2, x^2 + 3x + 1])

   @test isa(k, PolyElem)

   l = S(k)

   @test isa(l, PolyElem)

   m = S([1, 2, 3])

   @test isa(m, PolyElem)

   n = S([ZZ(1), ZZ(2), ZZ(3)])

   @test isa(n, PolyElem)

   @test x in [x, y]
   @test x in [y, x]
   @test !(x in [y])

   @test x in keys(Dict(x => 1))
   @test !(y in keys(Dict(x => 1)))
end

@testset "Generic.Poly.iterators" begin
   R, x = PolynomialRing(ZZ, "x")

   C = collect(coefficients(R()))

   @test C == []

   C = collect(coefficients(R(1)))

   @test C == [R(1)]

   C = collect(coefficients(x + 2))

   @test C == [R(2), R(1)]

   C = collect(coefficients(x^2 + 2))

   @test C == [R(2), R(), R(1)]
end

@testset "Generic.Poly.printing" begin
   R, x = PolynomialRing(ZZ, "x")

   @test string(zero(R)) == "0"
   @test string(one(R)) == "1"
   @test string(x) == "x"
   @test string(5*x^5-3*x^3+2*x^2-x+1) == "5*x^5 - 3*x^3 + 2*x^2 - x + 1"

   S, y = PolynomialRing(R, "y")

   @test string(zero(S)) == "0"
   @test string(one(S)) == "1"
   @test string(y) == "y"
   @test string(x+y+1) == "y + x + 1"
end

@testset "Generic.Poly.rand" begin
   R, x = PolynomialRing(ZZ, "x")

   # TODO: test more than just the result type
   test_rand(R, -1:10, -10:10)
   test_rand(R, -1:10, make(ZZ, -10:10))
   for deg in 1:10
      test_rand(R, deg, -10:10)
   end

   S, y = PolynomialRing(R, "y")

   test_rand(S, -1:5, make(R, -1:10, make(ZZ, -10:10)))
   test_rand(S, -1:5, make(R, -1:10, -10:10))
   test_rand(S, -1:5, -1:10, -10:10)

   T, z = PolynomialRing(GF(7), "z")

   test_rand(T, -1:4)
end

@testset "Generic.Poly.similar" begin
   R, x = PolynomialRing(ZZ, "x")

   for iters = 1:10
      f = rand(R, 0:10, -10:10)

      g = similar(f, QQ, "y")
      h = similar(f, "y")
      k = similar(f)

      @test isa(g, PolyElem)
      @test isa(h, PolyElem)
      @test isa(k, PolyElem)

      @test base_ring(g) === QQ

      @test parent(g).S == :y
      @test parent(h).S == :y

      @test length(g) == 0
      @test length(h) == 0
      @test length(k) == 0

      @test parent(g) !== parent(f)
      @test parent(h) !== parent(f)
      @test parent(k) === parent(f)

      p = similar(f, cached=false)
      q = similar(f, "z", cached=false)
      r = similar(f, "z", cached=false)
      s = similar(f)
      t = similar(f)

      @test parent(p) === parent(f)
      @test parent(q) !== parent(r)
      @test parent(s) === parent(t)
   end
end

@testset "Generic.Poly.polynomial" begin
   f = polynomial(ZZ, [1, 2, 3], "y")

   @test isa(f, PolyElem)
   @test base_ring(f) === ZZ
   @test coeff(f, 0) == 1
   @test coeff(f, 2) == 3
   @test parent(f).S == :y

   g = polynomial(ZZ, [1, 2, 3])

   @test isa(g, PolyElem)
   @test base_ring(g) === ZZ
   @test coeff(g, 0) == 1
   @test coeff(g, 2) == 3
   @test parent(g).S == :x

   h = polynomial(ZZ, [1, 2, 3])
   k = polynomial(ZZ, [1, 2, 3], cached=false)
   m = polynomial(ZZ, [1, 2, 3], cached=false)

   @test parent(h) === parent(g)
   @test parent(k) !== parent(m)

   p = polynomial(ZZ, BigInt[])
   q = polynomial(ZZ, [])

   @test isa(p, PolyElem)
   @test isa(q, PolyElem)

   @test length(p) == 0
   @test length(q) == 0
end

@testset "Generic.Poly.zero" begin
   R, x = PolynomialRing(ZZ, "x")

   f = rand(R, 0:10, -10:10)

   g = zero(f, QQ, "y")
   h = zero(f, "y")
   k = zero(f)

   @test isa(g, PolyElem)
   @test isa(h, PolyElem)
   @test isa(k, PolyElem)

   @test length(g) == 0
   @test length(h) == 0
   @test length(k) == 0

   @test base_ring(g) === QQ

   @test parent(g).S == :y
   @test parent(h).S == :y

   @test parent(g) !== parent(f)
   @test parent(h) !== parent(f)
   @test parent(k) === parent(f)

   p = zero(f, cached=false)
   q = zero(f, "z", cached=false)
   r = zero(f, "z", cached=false)
   s = zero(f)
   t = zero(f)

   @test parent(p) === parent(f)
   @test parent(q) !== parent(r)
   @test parent(s) === parent(t)
end

@testset "Generic.Poly.manipulation" begin
   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   @test iszero(zero(S))

   @test isone(one(S))

   @test isgen(gen(S))

   @test isunit(one(S))

   @test ismonic(R(1))
   @test ismonic(x^2 + 3)
   @test ismonic(y^2 + x*y + 1)
   @test !ismonic(x^2*y + 3x + 1)
   @test !ismonic(x^2*y + 3x + 2y + 1)
   @test !ismonic(-R(1))
   @test !ismonic(-2x^2*y + 3x + 2y + 1)

   f = 2x*y + x^2 + 1

   @test leading_coefficient(f) == 2x

   @test trailing_coefficient(2x*y + x^2) == x^2

   @test constant_coefficient(y^2 + 2x) == 2x

   @test tail(2x*y + 2x + y + 1) == 2x + 1
   @test tail(R(3)) == 0
   @test tail(3x*y) == 0
   @test tail(R()) == 0

   @test iszero(constant_coefficient(R()))

   @test degree(f) == 1

   h = x*y^2 + (x + 1)*y + 3

   @test coeff(h, 2) == x

   @test length(h) == 3

   @test canonical_unit(-x*y + x + 1) == -1

   @test deepcopy(h) == h

   @test isterm_recursive(2*x*y^2)
   @test !isterm_recursive(2*(x + 1)*y^2)

   @test !isterm(2*x*y^2 + 1)
   @test isterm(2*x*y^2)

   @test !ismonomial_recursive(2*x*y^2)

   @test ismonomial(y^2)

   @test !ismonomial_recursive(2*x*y^2 + y + 1)
   @test !ismonomial(2*y^2)

   @test characteristic(R) == 0
end

@testset "Generic.Poly.deflation" begin
   R, x = PolynomialRing(ZZ, "x")

   for iter = 1:1000
      f = rand(R, -1:20, -5:5)

      shift, defl = deflation(f)

      @test inflate(deflate(f, shift, defl), shift, defl) == f

      shift = rand(0:10)
      defl = rand(1:10)

      g = inflate(f, shift, defl)
      @test f == deflate(g, shift, defl)

      @test deflate(f, 1) == deflate(f, 0, 1)
      n = rand(1:10)
      @test inflate(f, n) == inflate(f, 0, n)
   end
end

@testset "Generic.Poly.binary_ops" begin
   #  Exact ring
   R, x = PolynomialRing(ZZ, "x")
   for iter = 1:100
      f = rand(R, -1:10, -10:10)
      g = rand(R, -1:10, -10:10)
      h = rand(R, -1:10, -10:10)
      @test f + g == g + f
      @test f + (g + h) == (f + g) + h
      @test f*g == g*f
      @test f*(g + h) == f*g + f*h
      @test (f - h) + (g + h) == f + g
      @test (f + g)*(f - g) == f*f - g*g
      @test f - g == -(g - f)
   end

   # Fake finite field of char 7, degree 2
   S, y = PolynomialRing(GF(7), "y")
   F = ResidueField(S, y^2 + 6y + 3)
   a = F(y)
   R, x = PolynomialRing(F, "x")
   for iter = 1:100
      f = rand(R, -1:10, 0:1)
      g = rand(R, -1:10, 0:1)
      h = rand(R, -1:10, 0:1)
      @test f + g == g + f
      @test f + (g + h) == (f + g) + h
      @test f*g == g*f
      @test f*(g + h) == f*g + f*h
      @test (f - h) + (g + h) == f + g
      @test (f + g)*(f - g) == f*f - g*g
      @test f - g == -(g - f)
   end

   #  Inexact field
   R, x = PolynomialRing(RealField, "x")
   for iter = 1:100
      f = rand(R, -1:10, -1:1)
      g = rand(R, -1:10, -1:1)
      h = rand(R, -1:10, -1:1)
      @test isapprox(f + (g + h), (f + g) + h)
      @test isapprox(f*g, g*f)
      @test isapprox(f*(g + h), f*g + f*h)
      @test isapprox((f - h) + (g + h), f + g)
      @test isapprox((f + g)*(f - g), f*f - g*g)
      @test isapprox(f - g, -(g - f))
   end

   # Non-integral domain
   T = ResidueRing(ZZ, 6)
   R, x = T["x"]
   for iter = 1:100
      f = rand(R, -1:10, 0:5)
      g = rand(R, -1:10, 0:5)
      h = rand(R, -1:10, 0:5)
      @test f + (g + h) == (f + g) + h
      @test f*g == g*f
      @test f*(g + h) == f*g + f*h
      @test (f - h) + (g + h) == f + g
      @test (f + g)*(f - g) == f*f - g*g
      @test f - g == -(g - f)
   end
end

@testset "Generic.Poly.adhoc_binary" begin
   # Exact ring
   R, x = ZZ["x"]
   for iter = 1:500
      f = rand(R, -1:10, -10:10)
      c1 = rand(ZZ, -10:10)
      c2 = rand(ZZ, -10:10)
      d1 = rand(zz, -10:10)
      d2 = rand(zz, -10:10)

      @test c1*f - c2*f == (c1 - c2)*f
      @test c1*f + c2*f == (c1 + c2)*f
      @test d1*f - d2*f == (d1 - d2)*f
      @test d1*f + d2*f == (d1 + d2)*f

      @test f*c1 - f*c2 == f*(c1 - c2)
      @test f*c1 + f*c2 == f*(c1 + c2)
      @test f*d1 - f*d2 == f*(d1 - d2)
      @test f*d1 + f*d2 == f*(d1 + d2)
   end

   # Fake finite field of char 7, degree 2
   S, y = PolynomialRing(GF(7), "y")
   F = ResidueField(S, y^2 + 6y + 3)
   a = F(y)
   R, x = PolynomialRing(F, "x")
   for iter = 1:500
      f = rand(R, -1:10, 0:1)
      c1 = rand(ZZ, -10:10)
      c2 = rand(ZZ, -10:10)
      d1 = rand(zz, -10:10)
      d2 = rand(zz, -10:10)

      @test c1*f - c2*f == (c1 - c2)*f
      @test c1*f + c2*f == (c1 + c2)*f
      @test d1*f - d2*f == (d1 - d2)*f
      @test d1*f + d2*f == (d1 + d2)*f

      @test f*c1 - f*c2 == f*(c1 - c2)
      @test f*c1 + f*c2 == f*(c1 + c2)
      @test f*d1 - f*d2 == f*(d1 - d2)
      @test f*d1 + f*d2 == f*(d1 + d2)
   end

   # Inexact field
   R, x = RealField["x"]
   for iter = 1:500
      f = rand(R, -1:10, -1:1)
      c1 = rand(ZZ, -10:10)
      c2 = rand(ZZ, -10:10)
      d1 = rand(RealField, -1:1)
      d2 = rand(RealField, -1:1)

      @test isapprox(c1*f - c2*f, (c1 - c2)*f)
      @test isapprox(c1*f + c2*f, (c1 + c2)*f)
      @test isapprox(d1*f - d2*f, (d1 - d2)*f)
      @test isapprox(d1*f + d2*f, (d1 + d2)*f)

      @test isapprox(f*c1 - f*c2, f*(c1 - c2))
      @test isapprox(f*c1 + f*c2, f*(c1 + c2))
      @test isapprox(f*d1 - f*d2, f*(d1 - d2))
      @test isapprox(f*d1 + f*d2, f*(d1 + d2))
   end

   # Non-integral domain
   R = ResidueRing(ZZ, 6)
   S, x = R["x"]
   for iter = 1:500
      f = rand(S, -1:10, 0:5)
      c1 = rand(ZZ, -10:10)
      c2 = rand(ZZ, -10:10)
      d1 = rand(zz, -10:10)
      d2 = rand(zz, -10:10)
      a1 = rand(R, 0:5)
      a2 = rand(R, 0:5)

      @test a1*f - a2*f == (a1 - a2)*f
      @test a1*f + a2*f == (a1 + a2)*f
      @test c1*f - c2*f == (c1 - c2)*f
      @test c1*f + c2*f == (c1 + c2)*f
      @test d1*f - d2*f == (d1 - d2)*f
      @test d1*f + d2*f == (d1 + d2)*f

      @test f*a1 - f*a2 == f*(a1 - a2)
      @test f*a1 + f*a2 == f*(a1 + a2)
      @test f*c1 - f*c2 == f*(c1 - c2)
      @test f*c1 + f*c2 == f*(c1 + c2)
      @test f*d1 - f*d2 == f*(d1 - d2)
      @test f*d1 + f*d2 == f*(d1 + d2)
   end

   # Generic tower
   R, x = ZZ["x"]
   S, y = R["y"]
   for iter = 1:100
      f = rand(S, -1:10, 0:5, -10:10)
      c1 = rand(ZZ, -10:10)
      c2 = rand(ZZ, -10:10)
      d1 = rand(R, -1:5, -10:10)
      d2 = rand(R, -1:5, -10:10)

      @test c1*f - c2*f == (c1 - c2)*f
      @test c1*f + c2*f == (c1 + c2)*f
      @test d1*f - d2*f == (d1 - d2)*f
      @test d1*f + d2*f == (d1 + d2)*f

      @test f*c1 - f*c2 == f*(c1 - c2)
      @test f*c1 + f*c2 == f*(c1 + c2)
      @test f*d1 - f*d2 == f*(d1 - d2)
      @test f*d1 + f*d2 == f*(d1 + d2)
   end
end

@testset "Generic.Poly.comparison" begin
   # Exact ring
   R, x = ZZ["x"]
   for iter = 1:500
      f = rand(R, -1:10, -10:10)
      g = deepcopy(f)
      h = rand(R, 0:10, -10:10)

      @test f == g
      @test isequal(f, g)
      @test f != g + h
   end

   # Fake finite field of char 7, degree 2
   S, y = PolynomialRing(GF(7), "y")
   F = ResidueField(S, y^2 + 6y + 3)
   a = F(y)
   R, x = PolynomialRing(F, "x")
   for iter = 1:500
      f = rand(R, -1:10, 0:1)
      g = deepcopy(f)
      h = rand(R, 0:10, 0:1)

      @test f == g
      @test isequal(f, g)
      @test f != g + h
   end

   # Inexact field
   R, x = RealField["x"]
   for iter = 1:500
      f = rand(R, 0:10, -1:1)
      g = deepcopy(f)
      h = rand(R, 0:10, -1:1)

      @test f == g
      @test isequal(f, g)
      @test f != g + h
   end

   # Non-integral domain
   R = ResidueRing(ZZ, 6)
   S, x = R["x"]
   for iter = 1:500
      f = rand(S, -1:10, 0:5)
      g = deepcopy(f)
      h = rand(S, 0:10, 0:5)

      @test f == g
      @test isequal(f, g)
      @test f != g + h
   end
end

@testset "Generic.Poly.adhoc_comparison" begin
   # Exact ring
   R, x = ZZ["x"]
   for iter = 1:500
      f = rand(R, 0:10, -10:10)
      c1 = rand(ZZ, -10:10)
      d1 = rand(zz, -10:10)

      @test R(c1) == c1
      @test c1 == R(c1)
      @test R(d1) == d1
      @test d1 == R(d1)

      @test R(c1) != c1 + f
      @test c1 != R(c1) + f
      @test R(d1) != d1 + f
      @test d1 != R(d1) + f
   end

   # Fake finite field of char 7, degree 2
   S, y = PolynomialRing(GF(7), "y")
   F = ResidueField(S, y^2 + 6y + 3)
   a = F(y)
   R, x = PolynomialRing(F, "x")
   for iter = 1:500
      f = rand(R, 0:10, 0:1)
      c1 = rand(ZZ, -10:10)
      d1 = rand(zz, -10:10)

      @test R(c1) == c1
      @test c1 == R(c1)
      @test R(d1) == d1
      @test d1 == R(d1)

      @test R(c1) != c1 + f
      @test c1 != R(c1) + f
      @test R(d1) != d1 + f
      @test d1 != R(d1) + f
   end

   # Inexact field
   R, x = RealField["x"]
   for iter = 1:500
      f = rand(R, 0:10, -1:1)
      c1 = rand(ZZ, -10:10)
      d1 = rand(RealField, -1:1)

      @test R(c1) == c1
      @test c1 == R(c1)
      @test R(d1) == d1
      @test d1 == R(d1)

      @test R(c1) != c1 + f
      @test c1 != R(c1) + f
      @test R(d1) != d1 + f
      @test d1 != R(d1) + f
   end

   # Non-integral domain
   R = ResidueRing(ZZ, 6)
   S, x = R["x"]
   for iter = 1:500
      f = rand(S, 0:10, 0:5)
      c1 = rand(ZZ, -10:10)
      d1 = rand(zz, -10:10)
      a1 = rand(R, 0:5)

      @test S(a1) == a1
      @test a1 == S(a1)
      @test S(c1) == c1
      @test c1 == S(c1)
      @test S(d1) == d1
      @test d1 == S(d1)

      @test S(a1) != a1 + f
      @test a1 != S(a1) + f
      @test S(c1) != c1 + f
      @test c1 != S(c1) + f
      @test S(d1) != d1 + f
      @test d1 != S(d1) + f
   end

   # Generic tower
   R, x = ZZ["x"]
   S, y = R["y"]
   for iter = 1:100
      f = rand(S, 0:10, 0:5, -10:10)
      c1 = rand(ZZ, -10:10)
      d1 = rand(R, -1:5, -10:10)

      @test S(c1) == c1
      @test c1 == S(c1)
      @test S(d1) == d1
      @test d1 == S(d1)

      @test S(c1) != c1 + f
      @test c1 != S(c1) + f
      @test S(d1) != d1 + f
      @test d1 != S(d1) + f
   end
end

@testset "Generic.Poly.unary_ops" begin
   #  Exact ring
   R, x = PolynomialRing(ZZ, "x")
   for iter = 1:300
      f = rand(R, -1:10, -10:10)

      @test -(-f) == f
      @test iszero(f + (-f))
   end

   # Fake finite field of char 7, degree 2
   S, y = PolynomialRing(GF(7), "y")
   F = ResidueField(S, y^2 + 6y + 3)
   a = F(y)
   R, x = PolynomialRing(F, "x")
   for iter = 1:300
      f = rand(R, -1:10, 0:1)

      @test -(-f) == f
      @test iszero(f + (-f))
   end

   #  Inexact field
   R, x = PolynomialRing(RealField, "x")
   for iter = 1:300
      f = rand(R, -1:10, -1:1)

      @test -(-f) == f
      @test iszero(f + (-f))
   end

   # Non-integral domain
   T = ResidueRing(ZZ, 6)
   R, x = T["x"]
   for iter = 1:300
      f = rand(R, -1:10, 0:5)

      @test -(-f) == f
      @test iszero(f + (-f))
   end
end

@testset "Generic.Poly.truncation" begin
   #  Exact ring
   R, x = PolynomialRing(ZZ, "x")
   for iter = 1:300
      f = rand(R, -1:10, -10:10)
      g = rand(R, -1:10, -10:10)
      n = rand(0:20)

      @test truncate(f*g, n) == mullow(f, g, n)
      @test truncate(f*g, n) == mullow_karatsuba(f, g, n, rand(1:5))
   end

   for iter = 1:100
      lena = rand(0:10)
      lenb = rand(0:10)
      a = rand(R, lena:lena, -10:10)
      b = rand(R, lenb:lenb, -10:10)
      c = a*b
      for i = 0:length(c) - 1
         d = mulhigh_n(a, b, i)
         f = c - d
         @test length(f) < length(c) - i
         for j = length(c):length(c) - i
            @test coeff(d, j) == coeff(c, j)
         end
      end
   end

   for iter = 1:100
      lena = rand(-1:10)
      lenb = rand(1:10)
      a = rand(R, lena:lena, -10:10)
      b = R()
      while iszero(b)
         b = rand(R, lenb:lenb, -10:10)
      end
      c = a*b
      for i = 0:length(b)
         d = divhigh(c, b, i)
         f = a - d
         @test degree(f) < i
         for j = 0:i - 1
            @test coeff(f, j) == coeff(a, j)
         end
      end
   end

   for iter = 1:100
      lena = rand(-1:10)
      lenb = rand(1:10)
      a = rand(R, lena:lena, -10:10)
      b = R()
      while iszero(b)
         b = rand(R, lenb:lenb, -10:10)
      end
      c = a*b
      for i = 0:length(a)
         d = divexact_low(c, b, i)
         f = truncate(a, i)
         @test f == d
      end
   end

   # Fake finite field of char 7, degree 2
   S, y = PolynomialRing(GF(7), "y")
   F = ResidueField(S, y^2 + 6y + 3)
   a = F(y)
   R, x = PolynomialRing(F, "x")
   for iter = 1:300
      f = rand(R, -1:10, 0:1)
      g = rand(R, -1:10, 0:1)
      n = rand(0:20)

      @test truncate(f*g, n) == mullow(f, g, n)
      @test truncate(f*g, n) == mullow_karatsuba(f, g, n, rand(1:5))
   end

   #  Inexact field
   R, x = PolynomialRing(RealField, "x")
   for iter = 1:300
      f = rand(R, -1:10, -1:1)
      g = rand(R, -1:10, -1:1)
      n = rand(0:20)

      @test isapprox(truncate(f*g, n), mullow(f, g, n))
   end

   # Non-integral domain
   T = ResidueRing(ZZ, 6)
   R, x = T["x"]
   for iter = 1:300
      f = rand(R, -1:10, 0:5)
      g = rand(R, -1:10, 0:5)
      n = rand(0:20)

      r = mullow(f, g, n)

      @test truncate(f*g, n) == r
      @test truncate(f*g, n) == mullow_karatsuba(f, g, n, rand(1:10))
      @test r == 0 || !iszero(leading_coefficient(r))
   end
end

@testset "Generic.Poly.reverse" begin
   #  Exact ring
   R, x = ZZ["x"]
   for iter = 1:300
      f = rand(R, -1:10, -10:10)
      len = rand(length(f):12)
      frev = reverse(f, len)

      shift = 0
      for i = 1:len
         if coeff(f, i - 1) != 0
            break
         end
         shift += 1
      end

      @test length(frev) == len - shift
      @test f == reverse(frev, len)
   end

   f = rand(R, -1:10, -10:10)
   @test_throws DomainError reverse(f, -1)
   @test_throws DomainError reverse(f, -rand(2:100))

   # Fake finite field of char 7, degree 2
   S, y = PolynomialRing(GF(7), "y")
   F = ResidueField(S, y^2 + 6y + 3)
   a = F(y)
   R, x = PolynomialRing(F, "x")
   for iter = 1:300
      f = rand(R, -1:10, 0:1)
      len = rand(length(f):12)
      frev = reverse(f, len)

      shift = 0
      for i = 1:len
         if coeff(f, i - 1) != 0
            break
         end
         shift += 1
      end

      @test length(frev) == len - shift
      @test f == reverse(frev, len)
   end

   f = rand(R, -1:10, 0:1)
   @test_throws DomainError reverse(f, -1)
   @test_throws DomainError reverse(f, -rand(2:100))

   #  Inexact field
   R, x = PolynomialRing(RealField, "x")
   for iter = 1:300
      f = rand(R, -1:10, -1:1)
      len = rand(length(f):12)
      frev = reverse(f, len)

      shift = 0
      for i = 1:len
         if coeff(f, i - 1) != 0
            break
         end
         shift += 1
      end

      @test length(frev) == len - shift
      @test f == reverse(frev, len)
   end

   f = rand(R, -1:10, -1:1)
   @test_throws DomainError reverse(f, -1)
   @test_throws DomainError reverse(f, -rand(2:100))

   #  Non-integral domain
   T = ResidueRing(ZZ, 6)
   R, x = T["x"]
   for iter = 1:300
      f = rand(R, -1:10, 0:5)
      len = rand(length(f):12)
      frev = reverse(f, len)

      shift = 0
      for i = 1:len
         if coeff(f, i - 1) != 0
            break
         end
         shift += 1
      end

      @test length(frev) == len - shift
      @test f == reverse(frev, len)
   end

   f = rand(R, -1:10, 0:5)
   @test_throws DomainError reverse(f, -1)
   @test_throws DomainError reverse(f, -rand(2:100))
end

@testset "Generic.Poly.shift" begin
   # Exact ring
   R, x = ZZ["x"]
   for iter = 1:300
      f = rand(R, -1:10, -10:10)
      s = rand(0:10)
      g = s == 0 ? R() : rand(R, 0:s - 1, -10:10)

      @test shift_right(shift_left(f, s) + g, s) == f
      @test shift_left(f, s) == x^s*f
      @test length(shift_right(f, s)) == max(0, length(f) - s)
   end

   f = rand(R, -1:10, -10:10)
   @test_throws DomainError shift_right(f, -1)
   @test_throws DomainError shift_right(f, -rand(2:100))
   @test_throws DomainError shift_left(f, -1)
   @test_throws DomainError shift_left(f, -rand(2:100))

   # Fake finite field of char 7, degree 2
   S, y = PolynomialRing(GF(7), "y")
   F = ResidueField(S, y^2 + 6y + 3)
   a = F(y)
   R, x = PolynomialRing(F, "x")
   for iter = 1:300
      f = rand(R, -1:10, 0:1)
      s = rand(0:10)
      g = s == 0 ? R() : rand(R, 0:s - 1, 0:1)

      @test shift_right(shift_left(f, s) + g, s) == f
      @test shift_left(f, s) == x^s*f
      @test length(shift_right(f, s)) == max(0, length(f) - s)
   end

   f = rand(R, -1:10, 0:1)
   @test_throws DomainError shift_right(f, -1)
   @test_throws DomainError shift_right(f, -rand(2:100))
   @test_throws DomainError shift_left(f, -1)
   @test_throws DomainError shift_left(f, -rand(2:100))

   # Inexact field
   R, x = PolynomialRing(RealField, "x")
   for iter = 1:300
      f = rand(R, -1:10, -1:1)
      s = rand(0:10)
      g = s == 0 ? R() : rand(R, 0:s - 1, -1:1)

      @test shift_right(shift_left(f, s) + g, s) == f
      @test shift_left(f, s) == x^s*f
      @test length(shift_right(f, s)) == max(0, length(f) - s)
   end

   f = rand(R, -1:10, -1:1)
   @test_throws DomainError shift_right(f, -1)
   @test_throws DomainError shift_right(f, -rand(2:100))
   @test_throws DomainError shift_left(f, -1)
   @test_throws DomainError shift_left(f, -rand(2:100))

   # Non-integral domain
   T = ResidueRing(ZZ, 6)
   R, x = T["x"]
   for iter = 1:300
      f = rand(R, -1:10, 0:5)
      s = rand(0:10)
      g = s == 0 ? R() : rand(R, 0:s - 1, 0:5)

      @test shift_right(shift_left(f, s) + g, s) == f
      @test shift_left(f, s) == x^s*f
      @test length(shift_right(f, s)) == max(0, length(f) - s)
   end

   f = rand(R, -1:10, 0:5)
   @test_throws DomainError shift_right(f, -1)
   @test_throws DomainError shift_right(f, -rand(2:100))
   @test_throws DomainError shift_left(f, -1)
   @test_throws DomainError shift_left(f, -rand(2:100))
end

@testset "Generic.Poly.powering" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter = 1:10
      f = rand(R, -1:10, -10:10)
      r2 = one(R)

      for expn = 0:10
         r1 = f^expn

         @test r1 == r2

         r2 *= f
      end
   end

   f = rand(R, -1:10, -10:10)
   @test_throws DomainError f^-1
   @test_throws DomainError f^-rand(2:100)
   @test_throws DomainError pow_multinomial(f, -1)
   @test_throws DomainError pow_multinomial(f, -rand(2:100))

   # Fake finite field of char 7, degree 2
   S, y = PolynomialRing(GF(7), "y")
   F = ResidueField(S, y^2 + 6y + 3)
   a = F(y)
   R, x = PolynomialRing(F, "x")

   for iter = 1:10
      f = rand(R, -1:10, 0:1)
      r2 = one(R)

      for expn = 0:10
         r1 = f^expn

         @test r1 == r2

         r2 *= f
      end
   end

   f = rand(R, -1:10, 0:1)
   @test_throws DomainError f^-1
   @test_throws DomainError f^-rand(2:100)
   @test_throws DomainError pow_multinomial(f, -1)
   @test_throws DomainError pow_multinomial(f, -rand(2:100))

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter = 1:10
      f = rand(R, -1:10, -1:1)
      r2 = one(R)

      for expn = 0:4 # cannot set high power here
         r1 = f^expn

         @test isapprox(r1, r2)

         r2 *= f
      end
   end

   f = rand(R, -1:10, -1:1)
   @test_throws DomainError f^-1
   @test_throws DomainError f^-rand(2:100)
   @test_throws DomainError pow_multinomial(f, -1)
   @test_throws DomainError pow_multinomial(f, -rand(2:100))

   # Non-integral domain
   for iter = 1:10
      n = rand(2:26)

      Zn = ResidueRing(ZZ, n)
      R, x = PolynomialRing(Zn, "x")

      f = rand(R, -1:10, 0:n - 1)
      r2 = one(R)

      for expn = 0:10
         r1 = f^expn

         @test r1 == r2

         r2 *= f
      end
   end

   f = rand(R, -1:10, 0:rand(1:25))
   @test_throws DomainError f^-1
   @test_throws DomainError f^-rand(2:100)
   @test_throws DomainError pow_multinomial(f, -1)
   @test_throws DomainError pow_multinomial(f, -rand(2:100))
end

if false
   @testset "Generic.Poly.modular_arithmetic" begin
      # Exact ring
      R = ResidueRing(ZZ, 23)
      S, x = PolynomialRing(R, "x")

      for iter = 1:100
         f = rand(S, -1:5, 0:22)
         g = rand(S, -1:5, 0:22)
         h = rand(S, -1:5, 0:22)
         k = rand(S, 0:5, 0:22)

         @test mulmod(mulmod(f, g, k), h, k) == mulmod(f, mulmod(g, h, k), k)
      end

      for iter = 1:100
         f = rand(S, 0:5, 0:22)
         g = rand(S, 0:5, 0:22)
         while gcd(f, g) != 1
            f = rand(S, 0:5, 0:22)
            g = rand(S, 0:5, 0:22)
         end

         @test mulmod(invmod(f, g), f, g) == mod(S(1), g)
      end

      for iter = 1:100
         f = rand(S, -1:5, 0:22)
         g = rand(S, 0:5, 0:22)
         p = mod(one(S), g)

         for expn = 0:5
            r = powermod(f, expn, g)

            @test r == p

            p = mulmod(p, f, g)
         end
      end

      # Fake finite field of char 7, degree 2
      R, y = PolynomialRing(GF(7), "y")
      F = ResidueField(R, y^2 + 6y + 3)
      a = F(y)
      S, x = PolynomialRing(F, "x")

      for iter = 1:100
         f = rand(S, -1:5, 0:1)
         g = rand(S, -1:5, 0:1)
         h = rand(S, -1:5, 0:1)
         k = rand(S, 0:5, 0:1)

         @test mulmod(mulmod(f, g, k), h, k) == mulmod(f, mulmod(g, h, k), k)
      end

      for iter = 1:100
         f = rand(S, 0:5, 0:1)
         g = rand(S, 0:5, 0:1)
         while gcd(f, g) != 1
            f = rand(S, 0:5, 0:1)
            g = rand(S, 0:5, 0:1)
         end

         @test mulmod(invmod(f, g), f, g) == mod(S(1), g)
      end

      for iter = 1:100
         f = rand(S, -1:5, 0:1)
         g = rand(S, 0:5, 0:1)
         p = mod(one(S), g)

         for expn = 0:5
            r = powermod(f, expn, g)

            @test r == p

            p = mulmod(p, f, g)
         end
      end

      # Inexact field
      S, x = PolynomialRing(RealField, "x")

      for iter = 1:100
         f = rand(S, -1:5, -1:1)
         g = rand(S, -1:5, -1:1)
         h = rand(S, -1:5, -1:1)
         k = rand(S, 0:5, -1:1)

         @test isapprox(mulmod(mulmod(f, g, k), h, k), mulmod(f, mulmod(g, h, k), k))
      end

      for iter = 1:100
         f = rand(S, 0:5, -1:1)
         g = rand(S, 0:5, -1:1)
         while gcd(f, g) != 1
            f = rand(S, 0:5, -1:1)
            g = rand(S, 0:5, -1:1)
         end

         @test isapprox(mulmod(invmod(f, g), f, g), mod(S(1), g))
      end

      for iter = 1:100
         f = rand(S, -1:5, -1:1)
         g = rand(S, 0:5, -1:1)
         p = mod(one(S), g)

         for expn = 0:5
            r = powermod(f, expn, g)

            @test isapprox(r, p)

            p = mulmod(p, f, g)
         end
      end

      # Exact field
      R, x = PolynomialRing(QQ, "y")

      for iter = 1:10
         f = rand(R, -1:5, -10:10)
         g = rand(R, -1:5, -10:10)
         h = rand(R, -1:5, -10:10)
         k = rand(R, 0:5, -10:10)

         @test mulmod(mulmod(f, g, k), h, k) == mulmod(f, mulmod(g, h, k), k)
      end

      for iter = 1:10
         f = rand(R, 0:5, -10:10)
         g = rand(R, 0:5, -10:10)

         while gcd(f, g) != 1
            f = rand(R, 0:5, -10:10)
            g = rand(R, 0:5, -10:10)
         end

         @test mulmod(invmod(f, g), f, g) == mod(R(1), g)
      end

      for iter = 1:10
         f = rand(R, -1:5, -10:10)
         g = rand(R, 0:5, -10:10)
         p = mod(one(R), g)

         for expn = 0:5
            r = powermod(f, expn, g)

            @test r == p

            p = mulmod(p, f, g)
         end
      end
   end
end

@testset "Generic.Poly.exact_division" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter = 1:100
      f = rand(R, -1:10, -100:100)
      g = rand(R, 0:10, -100:100)

      @test divexact(f*g, g) == f
   end

   @test_throws ArgumentError divexact(x^2, x - 1)

   # Fake finite field of char 7, degree 2
   S, y = PolynomialRing(GF(7), "y")
   F = ResidueField(S, y^2 + 6y + 3)
   a = F(y)
   R, x = PolynomialRing(F, "x")

   for iter = 1:100
      f = rand(R, -1:10, 0:1)
      g = rand(R, 0:10, 0:1)

      @test divexact(f*g, g) == f
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter = 1:100
      f = rand(R, -1:10, -1:1)
      g = rand(R, 0:10, -1:1)

      @test isapprox(divexact(f*g, g), f)
   end

   # Characteristic p ring
   n = 23
   Zn = ResidueRing(ZZ, n)
   R, x = PolynomialRing(Zn, "x")

   for iter = 1:100
      f = rand(R, -1:10, 0:n - 1)
      g = rand(R, 0:10, 0:n - 1)

      @test divexact(f*g, g) == f
   end
end

@testset "Generic.Poly.adhoc_exact_division" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter = 1:100
      f = rand(R, -1:10, -100:100)
      g = ZZ()
      while g == 0
         g = rand(ZZ, -10:10)
      end

      @test divexact(f*g, g) == f
      @test divexact_left(f*g, g) == f
      @test divexact_right(f*g, g) == f

      h = 0
      while h == 0
         h = rand(-10:10)
      end

      @test divexact(f*h, h) == f
   end

   # Fake finite field of char 7, degree 2
   S, y = PolynomialRing(GF(7), "y")
   F = ResidueField(S, y^2 + 6y + 3)
   a = F(y)
   R, x = PolynomialRing(F, "x")

   for iter = 1:100
      f = rand(R, -1:10, 0:1)
      g = ZZ()
      while g == 0
         g = rand(ZZ, 1:6)
      end

      @test divexact(f*g, g) == f

      h = 0
      while h == 0
         h = rand(1:6)
      end

      @test divexact(f*h, h) == f
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter = 1:100
      f = rand(R, -1:10, -1:1)
      g = ZZ()
      while g == 0
         g = rand(RealField, -1:1)
      end

      @test isapprox(divexact(f*g, g), f)

      h = 0
      while h == 0
         h = rand(-10:10)
      end

      @test isapprox(divexact(f*h, h), f)
   end

   # Characteristic p ring
   n = 23
   Zn = ResidueRing(ZZ, n)
   R, x = PolynomialRing(Zn, "x")

   for iter = 1:100
      f = rand(R, -1:10, 0:22)
      g = rand(Zn, 1:22)

      @test divexact(f*g, g) == f

      h = 0
      while (h % n) == 0
         h = rand(-100:100)
      end

      @test divexact(f*h, h) == f
   end

   # Generic tower
   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   for iter = 1:100
      f = rand(S, -1:10, 0:10, -100:100)
      g = rand(R, 0:10, -100:100)

      @test divexact(f*g, g) == f

      h = ZZ()
      while h == 0
         h = rand(ZZ, -10:10)
      end

      @test divexact(f*h, h) == f
   end
end

@testset "Generic.Poly.euclidean_division" begin
   # Exact ring
   R = ResidueRing(ZZ, 23)
   S, x = PolynomialRing(R, "x")

   for iter = 1:100
      f = rand(S, -1:5, 0:22)
      g = rand(S, -1:5, 0:22)
      h = rand(S, 0:5, 0:22)

      @test mod(f + g, h) == mod(f, h) + mod(g, h)
   end

   for iter = 1:10
      f = rand(S, -1:5, 0:22)
      g = rand(S, 0:5, 0:22)

      q, r = divrem(f, g)
      @test q*g + r == f

      @test mod(f, g) == r
   end

   # Fake finite field of char 7, degree 2
   R, y = PolynomialRing(GF(7), "y")
   F = ResidueField(R, y^2 + 6y + 3)
   a = F(y)
   S, x = PolynomialRing(F, "x")

   for iter = 1:100
      f = rand(S, -1:5, 0:1)
      g = rand(S, -1:5, 0:1)
      h = rand(S, 0:5, 0:1)

      @test mod(f + g, h) == mod(f, h) + mod(g, h)
   end

   for iter = 1:10
      f = rand(S, -1:5, 0:1)
      g = rand(S, 0:5, 0:1)

      q, r = divrem(f, g)
      @test q*g + r == f

      @test mod(f, g) == r
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter = 1:100
      f = rand(R, -1:5, -1:1)
      g = rand(R, -1:5, -1:1)
      h = rand(R, 0:5, -1:1)

      @test isapprox(mod(f + g, h), mod(f, h) + mod(g, h))
   end

   for iter = 1:10
      f = rand(R, -1:5, -1:1)
      g = rand(R, 0:5, -1:1)

      q, r = divrem(f, g)
      @test isapprox(q*g + r, f)

      @test isapprox(mod(f, g), r)
   end

   # Exact field
   R, x = PolynomialRing(QQ, "x")

   for iter = 1:100
      f = rand(R, -1:5, -10:10)
      g = rand(R, -1:5, -10:10)
      h = rand(R, 0:5, -10:10)

      @test mod(f + g, h) == mod(f, h) + mod(g, h)
   end

   for iter = 1:10
      f = rand(R, -1:5, -10:10)
      g = rand(R, 0:5, -10:10)

      q, r = divrem(f, g)
      @test q*g + r == f

      @test mod(f, g) == r
   end
end

@testset "Generic.Poly.ad_hoc_euclidean_division" begin
   # Exact field
   R, x = PolynomialRing(QQ, "x")

   for iter = 1:100
      f = rand(R, -1:5, -10:10)
      c = QQ()
      while iszero(c)
         c = rand(QQ, -10:10)
      end

      @test div(f, c)*c == f
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter = 1:100
      f = rand(R, -1:5, -1:1)
      c = RealField()
      while isapprox(c, RealField())
         c = rand(RealField, -1:1)
      end

      @test isapprox(div(f, c)*c, f)
   end

   # Residue ring
   S = ResidueRing(ZZ, 23)
   R, x = PolynomialRing(S, "x")

   for iter = 1:100
      f = rand(R, -1:5, 0:22)
      c = R()
      while iszero(c)
         c = rand(S, 0:22)
      end

      @test div(f, c)*c == f
   end
end

@testset "Generic.Poly.pseudodivision" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter = 1:100
      f = rand(R, -1:5, -10:10)
      g = rand(R, 0:5, -10:10)

      q, r = pseudodivrem(f, g)

      if length(f) < length(g)
         @test f == r && q == 0
      else
         @test q*g + r == f*leading_coefficient(g)^(length(f) - length(g) + 1)
      end

      @test pseudorem(f, g) == r
   end

   # Characteristic p ring
   Zn = ResidueRing(ZZ, 23)
   R, x = PolynomialRing(Zn, "x")

   for iter = 1:100
      f = rand(R, -1:5, 0:22)
      g = rand(R, 0:5, 0:22)

      q, r = pseudodivrem(f, g)

      if length(f) < length(g)
         @test f == r && q == 0
      else
         @test q*g + r == f*leading_coefficient(g)^(length(f) - length(g) + 1)
      end

      @test pseudorem(f, g) == r
   end
end

@testset "Generic.Poly.content_primpart_gcd" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter = 1:100
      f = rand(R, -1:10, -10:10)

      g = R()
      while g == 0
         g = rand(ZZ, -10:10)
      end

      @test content(f*g) == divexact(g, canonical_unit(g))*content(f)

      @test primpart(f*g) == canonical_unit(g)*primpart(f)
   end

   for iter = 1:20
      f = rand(R, -1:10, -10:10)
      g = rand(R, -1:10, -10:10)
      h = rand(R, 0:10, -10:10)

      @test gcd(f*h, g*h) == divexact(h, canonical_unit(leading_coefficient(h)))*gcd(f, g)

      @test lcm(f, h) == divexact(f*h, gcd(f, h))
   end

   # Exact field
   R, x = PolynomialRing(QQ, "x")

   for iter = 1:100
      f = rand(R, -1:5, -10:10)

      g = QQ()
      while g == 0
         g = rand(QQ, -10:10)
      end

      @test content(f*g) == content(f)*gcd(g, QQ()) # must normalise g correctly

      @test primpart(f*g) == primpart(f)*divexact(g, gcd(g, QQ()))
   end

   for iter = 1:20
      f = rand(R, -1:5, -10:10)
      g = rand(R, -1:5, -10:10)
      h = rand(R, 0:5, -10:10)

      @test gcd(f*h, g*h) == inv(leading_coefficient(h))*h*gcd(f, g)
   end

   for iter = 1:10
      f = rand(R, 0:5, -10:10)
      g = rand(R, 0:5, -10:10)
      while gcd(f, g) != 1
         f = rand(R, 0:5, -10:10)
         g = rand(R, 0:5, -10:10)
      end

      d, inv = gcdinv(f, g)

      @test d == gcd(f, g)

      @test mod(f*inv, g) == mod(R(1), g)
   end

   # Characteristic p ring
   R = ResidueRing(ZZ, 23)
   S, x = PolynomialRing(R, "x")

   for iter = 1:100
      f = rand(S, -1:10, 0:22)
      g = rand(R, 1:22)

      @test content(f*g) == divexact(g, canonical_unit(g))*content(f)

      @test primpart(f*g) == canonical_unit(g)*primpart(f)
   end

   for iter = 1:100
      f = rand(S, 0:5, 0:22)
      g = rand(S, 0:5, 0:22)
      while gcd(f, g) != 1
         f = rand(S, 0:5, 0:22)
         g = rand(S, 0:5, 0:22)
      end

      d, inv = gcdinv(f, g)

      @test d == gcd(f, g)

      @test mod(f*inv, g) == mod(S(1), g)
   end

   # Characteristic p field
   R = GF(23)
   S, x = PolynomialRing(R, "x")

   for iter = 1:100
      f = rand(S, -1:10)
      g = R()
      while g == 0
         g = rand(R)
      end

      @test content(f*g) == divexact(g, canonical_unit(g))*content(f)

      @test primpart(f*g) == canonical_unit(g)*primpart(f)
   end

   for iter = 1:100
      f = rand(S, 0:5)
      g = rand(S, 0:5)
      while gcd(f, g) != 1
         f = rand(S, 0:5)
         g = rand(S, 0:5)
      end

      d, inv = gcdinv(f, g)

      @test d == gcd(f, g)

      @test mod(f*inv, g) == mod(S(1), g)
   end
end

@testset "Generic.Poly.evaluation" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter in 1:10
      f = rand(R, -1:4, -10:10)
      g = rand(R, -1:4, -10:10)

      d = rand(ZZ, -10:10)

      @test evaluate(g, evaluate(f, d)) == evaluate(subst(g, f), d)
   end

   for iter in 1:10
      f = rand(R, -1:4, -10:10)
      g = rand(R, -1:4, -10:10)

      d = rand(-10:10)

      @test evaluate(g, evaluate(f, d)) == evaluate(subst(g, f), d)
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter in 1:10
      f = rand(R, -1:4, 0:1)
      g = rand(R, -1:4, 0:1)

      d = rand(RealField, 0:1)

      @test isapprox(evaluate(g, evaluate(f, d)), evaluate(subst(g, f), d))
   end

   for iter in 1:10
      f = rand(R, -1:4, 0:1)
      g = rand(R, -1:4, 0:1)

      d = rand(-10:10)

      @test isapprox(evaluate(g, evaluate(f, d)), evaluate(subst(g, f), d))
   end

   # Non-integral domain
   Zn = ResidueRing(ZZ, 23)
   R, x = PolynomialRing(Zn, "x")

   for iter in 1:10
      f = rand(R, -1:4, 0:22)
      g = rand(R, -1:4, 0:22)

      d = rand(Zn, 0:22)

      @test evaluate(g, evaluate(f, d)) == evaluate(subst(g, f), d)
   end

   for iter in 1:10
      f = rand(R, -1:4, 0:22)
      g = rand(R, -1:4, 0:22)

      d = rand(-100:100)

      @test evaluate(g, evaluate(f, d)) == evaluate(subst(g, f), d)
   end
end

@testset "Generic.Poly.composition" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter in 1:10
      f = rand(R, -1:5, -10:10)
      g = rand(R, -1:5, -10:10)
      h = rand(R, -1:5, -10:10)

      @test compose(f, compose(g, h)) == compose(compose(f, g), h)
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter in 1:10
      f = rand(R, -1:5, 0:1)
      g = rand(R, -1:5, 0:1)
      h = rand(R, -1:5, 0:1)

      @test isapprox(compose(f, compose(g, h)), compose(compose(f, g), h))
   end

   # Non-integral domain
   Zn = ResidueRing(ZZ, 6)
   R, x = PolynomialRing(Zn, "x")

   for iter in 1:10
      f = rand(R, -1:5, 0:5)
      g = rand(R, -1:5, 0:5)
      h = rand(R, -1:5, 0:5)

      @test compose(f, compose(g, h)) == compose(compose(f, g), h)
   end
end

@testset "Generic.Poly.derivative" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter in 1:10
      f = rand(R, -1:4, -100:100)
      g = rand(R, -1:4, -100:100)

      @test derivative(f + g) == derivative(g) + derivative(f)

      @test derivative(g*f) == derivative(g)*f + derivative(f)*g
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter in 1:10
      f = rand(R, -1:4, 0:1)
      g = rand(R, -1:4, 0:1)

      @test isapprox(derivative(f + g), derivative(g) + derivative(f))

      @test isapprox(derivative(g*f), derivative(g)*f + derivative(f)*g)
   end

   # Non-integral domain
   Zn = ResidueRing(ZZ, 6)
   R, x = PolynomialRing(Zn, "x")

   for iter in 1:10
      f = rand(R, -1:4, 0:5)
      g = rand(R, -1:4, 0:5)

      @test derivative(f + g) == derivative(g) + derivative(f)

      @test derivative(g*f) == derivative(g)*f + derivative(f)*g
   end
end

@testset "Generic.Poly.integral" begin
   # Exact field
   R, x = PolynomialRing(QQ, "x")

   for iter in 1:10
      f = rand(R, -1:10, -100:100)

      @test derivative(integral(f)) == f

      g = rand(R, -1:2, -100:100)

      @test integral(f + g) == integral(g) + integral(f)
      @test integral(f)*integral(g) == integral(integral(f)*g + integral(g)*f)
   end

   # Characteristic p ring
   Zn = ResidueRing(ZZ, 23)
   R, x = PolynomialRing(Zn, "x")

   for iter in 1:10
      f = rand(R, -1:10, 0:22)

      @test derivative(integral(f)) == f

      g = rand(R, -1:10, 0:22)

      @test integral(f + g) == integral(g) + integral(f)
      @test integral(f)*integral(g) == integral(integral(f)*g + integral(g)*f)
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter in 1:10
      f = rand(R, -1:10, 0:1)

      @test isapprox(derivative(integral(f)), f)

      g = rand(R, -1:10, 0:1)

      @test isapprox(integral(f + g), integral(g) + integral(f))
      @test isapprox(integral(f)*integral(g), integral(integral(f)*g + integral(g)*f))
   end
end

@testset "Generic.Poly.sylvester_matrix" begin
   R, x = PolynomialRing(ZZ, "x")

   for iter in 1:10
      f = rand(R, 1:5, -10:10)
      g = rand(R, 1:5, -10:10)

      d1 = degree(f)
      d2 = degree(g)

      f1 = rand(R, 0:d2-1, -10:10)
      g1 = rand(R, 0:d1-1, -10:10)

      w = matrix(ZZ, 1, d2, [coeff(f1, d2 - i) for i in 1:d2])
      w = hcat(w, matrix(ZZ, 1, d1, [coeff(g1, d1 - i) for i in 1:d1]))

      h = f1 * f + g1 * g

      v = matrix(ZZ, 1, d1 + d2, [coeff(h, d1 + d2 - i) for i in 1:d1 + d2])
      M = sylvester_matrix(f, g)
      @test v == w * M
   end
end

@testset "Generic.Poly.resultant" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter in 1:10
      f = rand(R, -1:5, -10:10)
      g = rand(R, -1:5, -10:10)
      h = rand(R, -1:5, -10:10)

      @test resultant(f*g, h) == resultant(f, h) * resultant(g, h)
      @test resultant(f, g*h) == resultant(f, g) * resultant(f, h)

      @test resultant(f, g) == resultant_subresultant(f, g)
      @test resultant_ducos(f, g) == resultant_subresultant(f, g)
   end

   # Exact field
   R, x = PolynomialRing(QQ, "x")

   for iter in 1:10
      f = rand(R, -1:5, -10:10)
      g = rand(R, -1:5, -10:10)
      h = rand(R, -1:5, -10:10)

      @test resultant(f*g, h) == resultant(f, h) * resultant(g, h)
      @test resultant(f, g*h) == resultant(f, g) * resultant(f, h)

      @test resultant(f, g) == resultant_subresultant(f, g)
      @test resultant_ducos(f, g) == resultant_subresultant(f, g)
      @test resultant(f, g) == Generic.resultant_lehmer(f, g)
   end

   # Characteristic p ring
   Zn = ResidueRing(ZZ, 23)
   R, x = PolynomialRing(Zn, "x")

   for iter in 1:10
      f = rand(R, -1:5, 0:22)
      g = rand(R, -1:5, 0:22)
      h = rand(R, -1:5, 0:22)

      @test resultant(f*g, h) == resultant(f, h)*resultant(g, h)
      @test resultant(f, g*h) == resultant(f, g)*resultant(f, h)

      @test resultant(f, g) == resultant_subresultant(f, g)
      @test resultant_ducos(f, g) == resultant_subresultant(f, g)
   end

   # Characteristic p field
   R, x = PolynomialRing(GF(23), "x")

   for iter in 1:10
      f = rand(R, -1:5)
      g = rand(R, -1:5)
      h = rand(R, -1:5)

      @test resultant(f*g, h) == resultant(f, h)*resultant(g, h)
      @test resultant(f, g*h) == resultant(f, g)*resultant(f, h)

      @test resultant(f, g) == resultant_subresultant(f, g)
      @test resultant_ducos(f, g) == resultant_subresultant(f, g)
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter in 1:10
      f = rand(R, -1:5, 0:1)
      g = rand(R, -1:5, 0:1)
      h = rand(R, -1:5, 0:1)

      @test isapprox(resultant(f*g, h), resultant(f, h)*resultant(g, h))
      @test isapprox(resultant(f, g*h), resultant(f, g)*resultant(f, h))

      @test isapprox(resultant(f, g), resultant_subresultant(f, g))
   end

   # Non-integral domain
   Zn = ResidueRing(ZZ, 6)
   R, x = PolynomialRing(Zn, "x")

   for iter in 1:10
      f = rand(R, -1:5, 0:5)
      g = rand(R, -1:5, 0:5)
      h = rand(R, -1:5, 0:5)

      @test leading_coefficient(f)*leading_coefficient(g) == 0 ||
            resultant(f*g, h) == resultant(f, h)*resultant(g, h)
      @test leading_coefficient(g)*leading_coefficient(h) == 0 ||
            resultant(f, g*h) == resultant(f, g)*resultant(f, h)
   end
end

@testset "Generic.Poly.discriminant" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter = 1:100
      f = rand(R, 1:5, -10:10)
      g = rand(R, 1:5, -10:10)

      # See http://www2.math.uu.se/~svante/papers/sjN5.pdf 3.10
      # The identity on Wikipedia is incorrect as of 07.10.2017
      @test discriminant(f*g) == discriminant(f)*discriminant(g)*resultant(g, f)^2
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter = 1:100
      f = rand(R, 1:10, 0:1)
      g = rand(R, 1:10, 0:1)

      # See http://www2.math.uu.se/~svante/papers/sjN5.pdf 3.10
      # The identity on Wikipedia is incorrect as of 07.10.2017
      @test isapprox(discriminant(f*g), discriminant(f)*discriminant(g)*resultant(g, f)^2)
   end

#   TODO: Fix issue #291
#   # Non-integral domain
#   Zn = ResidueRing(ZZ, 6)
#   R, x = PolynomialRing(Zn, "x")
#
#   for iter = 1:100
#      f = R()
#      g = R()
#      while length(f) < 2 || length(g) < 2
#         f = rand(R, 1:5, 0:5)
#         g = rand(R, 1:5, 0:5)
#      end
#
#      # See http://www2.math.uu.se/~svante/papers/sjN5.pdf 3.10
#      # The identity on Wikipedia is incorrect as of 07.10.2017
#      @test discriminant(f*g) == discriminant(f)*discriminant(g)*resultant(g, f)^2
#   end
end

@testset "Generic.Poly.resx" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter in 1:100
      f = rand(R, 1:5, -10:10)
      g = rand(R, 1:5, -10:10)
      r, u, v = resx(f, g)

      @test u*f + v*g == r
      @test r == resultant(f, g)

      h = rand(R, -1:5, -10:10)
      r, u, v = resx(f*h, g*h)

      @test (u*f + v*g)*h == r
   end

   # Exact field
   R, x = PolynomialRing(QQ, "x")

   for iter in 1:100
      f = rand(R, 1:5, -10:10)
      g = rand(R, 1:5, -10:10)
      r, u, v = resx(f, g)

      @test u*f + v*g == r
      @test r == resultant(f, g)

      h = rand(R, -1:5, -10:10)
      r, u, v = resx(f*h, g*h)

      @test (u*f + v*g)*h == r
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter in 1:100
      f = rand(R, 1:5, 0:1)
      g = rand(R, 1:5, 0:1)
      r, u, v = resx(f, g)

      @test isapprox(u*f + v*g, r)
      @test isapprox(r, resultant(f, g))

      h = rand(R, -1:5, 0:1)
      r, u, v = resx(f*h, g*h)

      @test isapprox((u*f + v*g)*h, r)
   end

   # Characteristic p ring
   Zn = ResidueRing(ZZ, 23)
   R, x = PolynomialRing(Zn, "x")

   for iter in 1:100
      f = rand(R, 1:5, 0:22)
      g = rand(R, 1:5, 0:22)
      r, u, v = resx(f, g)

      @test u*f + v*g == r
      @test r == resultant(f, g)

      h = rand(R, 0:5, 0:22)
      r, u, v = resx(f*h, g*h)

      @test (u*f + v*g)*h == r
   end

   # Characteristic p field
   R, x = PolynomialRing(GF(23), "x")

   for iter in 1:100
      f = rand(R, 1:5)
      g = rand(R, 1:5)
      r, u, v = resx(f, g)

      @test u*f + v*g == r
      @test r == resultant(f, g)

      h = rand(R, -1:5)
      r, u, v = resx(f*h, g*h)

      @test (u*f + v*g)*h == r
   end

#   TODO: Fix issue #293
#   Test will cause impossible inverse in the mean time
#
#   # Non-integral domain
#   Zn = ResidueRing(ZZ, 6)
#   R, x = PolynomialRing(Zn, "x")
#
#   for iter in 1:100
#      f = R()
#      g = R()
#      while length(f) <= 1 && length(g) <= 1
#         f = rand(R, 0:5, 0:5)
#         g = rand(R, 0:5, 0:5)
#      end
#      r, u, v = resx(f, g)
#
#      @test u*f + v*g == r
#      @test r == resultant(f, g)
#
#      h = R()
#      h = rand(R, 0:5, 0:5)
#      r, u, v = resx(f*h, g*h)
#
#      @test (u*f + v*g)*h == r
#   end
end

@testset "Generic.Poly.gcdx" begin
   # Exact field
   R, x = PolynomialRing(QQ, "x")

   for iter in 1:100
      f = rand(R, 1:5, -10:10)
      g = rand(R, 1:5, -10:10)
      r, u, v = gcdx(f, g)

      @test u*f + v*g == r
      @test r == gcd(f, g)

      h = R()
      h = rand(R, -1:5, -10:10)
      r, u, v = gcdx(f*h, g*h)

      @test (u*f + v*g)*h == r
   end

   # Characteristic p ring
   Zn = ResidueRing(ZZ, 23)
   R, x = PolynomialRing(Zn, "x")

   for iter in 1:100
      f = rand(R, 1:5, 0:22)
      g = rand(R, 1:5, 0:22)
      r, u, v = gcdx(f, g)

      @test u*f + v*g == r
      @test r == gcd(f, g)

      h = rand(R, -1:5, 0:22)
      r, u, v = gcdx(f*h, g*h)

      @test (u*f + v*g)*h == r
   end

   # Characteristic p field
   R, x = PolynomialRing(GF(23), "x")

   for iter in 1:100
      f = rand(R, 1:5)
      g = rand(R, 1:5)
      r, u, v = gcdx(f, g)

      @test u*f + v*g == r
      @test r == gcd(f, g)

      h = rand(R, -1:5)
      r, u, v = gcdx(f*h, g*h)

      @test (u*f + v*g)*h == r
   end

   # Fake finite field of char 7, degree 2
   S, y = PolynomialRing(GF(7), "y")
   F = ResidueField(S, y^2 + 6y + 3)
   a = F(y)
   R, x = PolynomialRing(F, "x")

   for iter in 1:100
      f = rand(R, 1:5, 0:1)
      g = rand(R, 1:5, 0:1)
      r, u, v = gcdx(f, g)

      @test u*f + v*g == r
      @test r == gcd(f, g)

      h = rand(R, -1:5, 0:1)
      r, u, v = gcdx(f*h, g*h)

      @test (u*f + v*g)*h == r
   end
end

@testset "Generic.Poly.newton_representation" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter in 1:10
      f = rand(R, -1:10, -100:100)

      g = deepcopy(f)
      roots = BigInt[rand(ZZ, -10:10) for i in 1:length(f)]
      monomial_to_newton!(g.coeffs, roots)
      newton_to_monomial!(g.coeffs, roots)

      @test f == g
   end

   # Exact field
   R, x = PolynomialRing(QQ, "x")

   for iter in 1:10
      f = rand(R, -1:10, -100:100)

      g = deepcopy(f)
      roots = Rational{BigInt}[rand(QQ, -10:10) for i in 1:length(f)]
      monomial_to_newton!(g.coeffs, roots)
      newton_to_monomial!(g.coeffs, roots)

      @test f == g
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter in 1:10
      f = rand(R, -1:10, 0:1)

      g = deepcopy(f)
      roots = BigFloat[rand(RealField, 0:1) for i in 1:length(f)]
      monomial_to_newton!(g.coeffs, roots)
      newton_to_monomial!(g.coeffs, roots)

      @test isapprox(f, g)
   end

   # Characteristic p ring
   Zn = ResidueRing(ZZ, 23)
   R, x = PolynomialRing(Zn, "x")

   for iter in 1:10
      f = rand(R, -1:10, 0:22)

      g = deepcopy(f)
      roots = elem_type(Zn)[rand(Zn, 0:22) for i in 1:length(f)]
      monomial_to_newton!(g.coeffs, roots)
      newton_to_monomial!(g.coeffs, roots)

      @test f == g
   end

   # Characteristic p ring
   K = GF(23)
   R, x = PolynomialRing(K, "x")

   for iter in 1:10
      f = rand(R, -1:10)

      g = deepcopy(f)
      roots = elem_type(K)[rand(K) for i in 1:length(f)]
      monomial_to_newton!(g.coeffs, roots)
      newton_to_monomial!(g.coeffs, roots)

      @test f == g
   end

   # Non-integral domain
   Zn = ResidueRing(ZZ, 6)
   R, x = PolynomialRing(Zn, "x")

   for iter in 1:10
      f = rand(R, -1:10, 0:5)

      g = deepcopy(f)
      roots = elem_type(Zn)[rand(Zn, 0:5) for i in 1:length(f)]
      monomial_to_newton!(g.coeffs, roots)
      newton_to_monomial!(g.coeffs, roots)

      @test f == g
   end
end

@testset "Generic.Poly.interpolation" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for iter in 1:10
      p = rand(R, 0:10, -10:10)

      xs = BigInt[i for i in 1:length(p)]
      ys = [p(i) for i in 1:length(p)]

      f = interpolate(R, xs, ys)

      @test f == p
   end

   # Exact field
   R, x = PolynomialRing(QQ, "x")

   for iter in 1:10
      p = rand(R, 0:10, -10:10)

      xs = Rational{BigInt}[i for i in 1:length(p)]
      ys = [p(i) for i in 1:length(p)]

      f = interpolate(R, xs, ys)

      @test f == p
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for iter in 1:10
      p = rand(R, 0:10, 0:1)

      xs = BigFloat[i for i in 1:length(p)]
      ys = [p(i) for i in 1:length(p)]

      f = interpolate(R, xs, ys)

      @test isapprox(f, p)
   end

   # Characteristic p ring
   Zn = ResidueRing(ZZ, 23)
   R, x = PolynomialRing(Zn, "x")

   for iter in 1:10
      p = rand(R, 0:10, 0:22)

      xs = elem_type(Zn)[Zn(i) for i in 1:length(p)]
      ys = [p(i) for i in 1:length(p)]

      f = interpolate(R, xs, ys)

      @test f == p
   end

   # Characteristic p field
   K = GF(23)
   R, x = PolynomialRing(K, "x")

   for iter in 1:10
      p = rand(R, 0:10)

      xs = elem_type(K)[K(i) for i in 1:length(p)]
      ys = [p(i) for i in 1:length(p)]

      f = interpolate(R, xs, ys)

      @test f == p
   end

#   TODO: Fix issue #294 (if possible)
#   # Non-integral domain
#   Zn = ResidueRing(ZZ, 6)
#   R, x = PolynomialRing(Zn, "x")
#
#   for iter in 1:10
#      p = R()
#      while p == 0
#         p = rand(R, 0:10, 0:5)
#      end
#
#      xs = elem_type(Zn)[Zn(i) for i in 1:length(p)]
#      ys = [p(i) for i in 1:length(p)]
#
#      f = interpolate(R, xs, ys)
#
#      @test f == p
#   end
end

@testset "Generic.Poly.special" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   for n in 10:20
      T = chebyshev_t(n, x)
      dT = derivative(T)
      ddT = derivative(dT)

      @test (1 - x^2)*ddT + n^2*T == x*dT

      U = chebyshev_u(n - 1, x)
      dU = derivative(U)
      ddU = derivative(dU)

      @test (1 - x^2)*ddU + (n-1)*(n+1)*U == 3*x*dU

      @test T^2 == 1 + (x^2 - 1)*U^2
   end

   # Exact field
   R, x = PolynomialRing(QQ, "x")

   for n in 10:20
      T = chebyshev_t(n, x)
      dT = derivative(T)
      ddT = derivative(dT)

      @test (1 - x^2)*ddT + n^2*T == x*dT

      U = chebyshev_u(n - 1, x)
      dU = derivative(U)
      ddU = derivative(dU)

      @test (1 - x^2)*ddU + (n-1)*(n+1)*U == 3*x*dU

      @test T^2 == 1 + (x^2 - 1)*U^2
   end

   # Inexact field
   R, x = PolynomialRing(RealField, "x")

   for n in 10:20
      T = chebyshev_t(n, x)
      dT = derivative(T)
      ddT = derivative(dT)

      @test (1 - x^2)*ddT + n^2*T == x*dT

      U = chebyshev_u(n - 1, x)
      dU = derivative(U)
      ddU = derivative(dU)

      @test (1 - x^2)*ddU + (n-1)*(n+1)*U == 3*x*dU

      @test T^2 == 1 + (x^2 - 1)*U^2
   end

   # Characteristic p ring
   Zn = ResidueRing(ZZ, 23)
   R, x = PolynomialRing(Zn, "x")

   for n in 10:20
      T = chebyshev_t(n, x)
      dT = derivative(T)
      ddT = derivative(dT)

      @test (1 - x^2)*ddT + n^2*T == x*dT

      U = chebyshev_u(n - 1, x)
      dU = derivative(U)
      ddU = derivative(dU)

      @test (1 - x^2)*ddU + (n-1)*(n+1)*U == 3*x*dU

      @test T^2 == 1 + (x^2 - 1)*U^2
   end

   # Characteristic p field
   R, x = PolynomialRing(GF(23), "x")

   for n in 10:20
      T = chebyshev_t(n, x)
      dT = derivative(T)
      ddT = derivative(dT)

      @test (1 - x^2)*ddT + n^2*T == x*dT

      U = chebyshev_u(n - 1, x)
      dU = derivative(U)
      ddU = derivative(dU)

      @test (1 - x^2)*ddU + (n-1)*(n+1)*U == 3*x*dU

      @test T^2 == 1 + (x^2 - 1)*U^2
   end

   # Non-integral domain
   Zn = ResidueRing(ZZ, 6)
   R, x = PolynomialRing(Zn, "x")

   for n in 10:20
      T = chebyshev_t(n, x)
      dT = derivative(T)
      ddT = derivative(dT)

      @test (1 - x^2)*ddT + n^2*T == x*dT

      U = chebyshev_u(n - 1, x)
      dU = derivative(U)
      ddU = derivative(dU)

      @test (1 - x^2)*ddU + (n-1)*(n+1)*U == 3*x*dU

      @test T^2 == 1 + (x^2 - 1)*U^2
   end
end

@testset "Generic.Poly.mul_karatsuba" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")

   f = x + y + 2z^2 + 1

   @test mul_karatsuba(f^10, f^10) == mul_classical(f^10, f^10)
   @test mul_karatsuba(f^10, f^30) == mul_classical(f^10, f^30)

   @test mul_karatsuba(f^10, f^30, 5) == mul_classical(f^10, f^30)
   @test mul_karatsuba(f^10, f^30, 10) == mul_classical(f^10, f^30)

end

@testset "Generic.Poly.mul_ks" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")

   f = x + y + 2z^2 + 1

   @test mul_ks(f^10, f^10) == mul_classical(f^10, f^10)
   @test mul_ks(f^10, f^30) == mul_classical(f^10, f^30)
end

@testset "Generic.Poly.remove_valuation" begin
   # Exact ring
   R, x = PolynomialRing(ZZ, "x")

   @test_throws ErrorException remove(R(1), R(0))
   @test_throws ErrorException remove(R(1), R(-1))
   @test_throws ErrorException remove(R(0), R(1))
   @test_throws ErrorException remove(R(0), R(2))

   for iter = 1:10
      d = true
      f = R()
      g = R()
      while d
         f = rand(R, 0:10, -10:10)
         g = rand(R, 0:10, -10:10)
         while isunit(g)
           g = rand(R, 0:10, -10:10)
         end

         d, q = divides(f, g)
         if d
           @test g * q == f
         end
      end

      s = rand(0:10)

      v, q = remove(f*g^s, g)

      @test valuation(f*g^s, g) == s
      @test q == f
      @test v == s

      v, q = divides(f*g, f)

      @test v
      @test q == g

      if length(f) > 1
         v, q = divides(f*g + 1, f)

         @test !v
      end
   end

   # Exact field
   R, x = PolynomialRing(QQ, "x")

   @test_throws ErrorException remove(R(1), R(0))
   @test_throws ErrorException remove(R(1), R(1))
   @test_throws ErrorException remove(R(0), R(x))
   @test_throws ErrorException remove(R(1), R(2))

   for iter = 1:10
      d = true
      f = R()
      g = R()
      while d
         f = rand(R, 0:10, -10:10)
         g = rand(R, 0:10, -10:10)
         while isunit(g)
            g = rand(R, 0:10, -10:10)
         end

         d, q = divides(f, g)
      end

      s = rand(0:10)

      v, q = remove(f*g^s, g)

      @test valuation(f*g^s, g) == s
      @test q == f
      @test v == s

      v, q = divides(f*g, f)

      @test v
      @test q == g

      if length(f) > 1
         v, q = divides(f*g + 1, f)

         @test !v
      end
   end

   # Characteristic p ring
   Zn = ResidueRing(ZZ, 23)
   R, x = PolynomialRing(Zn, "x")

   for iter = 1:10
      d = true
      f = R()
      g = R()
      while d
         f = rand(R, 0:10, 0:22)
         g = rand(R, 0:10, 0:22)
         while isunit(g)
            f = rand(R, 0:10, 0:22)
            g = rand(R, 0:10, 0:22)
         end

         d, q = divides(f, g)
      end

      s = rand(0:10)

      v, q = remove(f*g^s, g)

      @test valuation(f*g^s, g) == s
      @test q == f
      @test v == s

      v, q = divides(f*g, f)

      @test v
      @test q == g

      if length(f) > 1
         v, q = divides(f*g + 1, f)

         @test !v
      end
   end

   # Characteristic p field
   R, x = PolynomialRing(GF(23), "x")

   for iter = 1:10
      d = true
      f = R()
      g = R()
      while d
         f = rand(R, 0:10)
         g = rand(R, 0:10)
         while isunit(g)
            g = rand(R, 0:10)
         end

         d, q = divides(f, g)
      end

      s = rand(0:10)

      v, q = remove(f*g^s, g)

      @test valuation(f*g^s, g) == s
      @test q == f
      @test v == s

      v, q = divides(f*g, f)

      @test v
      @test q == g

      if length(f) > 1
         v, q = divides(f*g + 1, f)

         @test !v
      end
   end
end

@testset "Generic.Poly.square_root" begin
   # Exact ring
   S, x = PolynomialRing(ZZ, "x")
   for iter = 1:10
      f = rand(S, -1:20, -20:20)

      p = f^2

      @test issquare(p)

      q = sqrt(f^2)

      @test q^2 == f^2

      q = sqrt(f^2, false)

      @test q^2 == f^2

      if f != 0
         @test_throws ErrorException sqrt(f^2*x)
      end
   end

   # Exact field
   S, x = PolynomialRing(QQ, "x")
   for iter = 1:10
      f = rand(S, -1:20, -20:20)

      p = f^2

      @test issquare(p)

      q = sqrt(f^2)

      @test q^2 == f^2

      q = sqrt(f^2, false)

      @test q^2 == f^2

      if f != 0
         @test_throws ErrorException sqrt(f^2*x)
      end
   end

   # Characteristic p field
   for p in [2, 7, 19, 65537, ZZ(2), ZZ(7), ZZ(19), ZZ(65537)]
      R = ResidueField(ZZ, p)

      S, x = PolynomialRing(R, "x")

      for iter = 1:10
         f = rand(S, -1:20, 0:Int(p))

         s = f^2

         @test issquare(s)

         q = sqrt(f^2)

         @test q^2 == f^2

         q = sqrt(f^2, false)

         @test q^2 == f^2

         if f != 0
            @test_throws ErrorException sqrt(f^2*x)
         end
      end
   end
end

@testset "Generic.Poly.generic_eval" begin
   R, x = PolynomialRing(ZZ, "x")

   for iter in 1:10
      f = rand(R, -1:2, -100:100)
      g = rand(R, -1:2, -100:100)
      h = rand(R, -1:2, -100:100)

      @test f(g(h)) == f(g)(h)
   end

   R, x = PolynomialRing(ZZ, "x")

   f = x
   b = a = QQ(13)
   for i in 1:5
      g = x^2 + rand(R, -1:1, -1:1)
      f = g(f)
      b = g(b)

      @test b == f(a)
   end
end

@testset "Generic.Poly.change_base_ring" begin
   Zx, x = PolynomialRing(ZZ,'x')
   @test 1 == map_coefficients(sqrt, x^0)
   p = Zx([i for i in 1:10])
   q = Zx([i for i in 10:-1:1])
   pq = p * q
   for R in [QQ,GF(2),GF(13),ZZ]
      pR = change_base_ring(R, p)
      qR = change_base_ring(R, q, parent = parent(pR))
      @test parent(qR) === parent(pR)
      pqR = change_base_ring(R, pq, parent = parent(pR))
      @test pR * qR == pqR
   end

   ps = map_coefficients(z -> z^2, p)
   @test ps == Zx([i^2 for i in 1:10])

   f = x^2 + 3x^3 + 2x^6
   @test map_coefficients(one, f) == x^2 + x^3 + x^6
   f2 = map_coefficients(t -> t+2, f)
   @test f2 == 3x^2 + 5x^3 + 4x^6
   for i in [0, 1, 4, 5]
      @test coeff(f2, i) !== coeff(f, i)
   end

   F = GF(11)
   P, y = PolynomialRing(F, 'x')
   @test map_coefficients(t -> F(t) + 2, f) == 3y^2 + 5y^3 + 4y^6
end

@testset "Generic.Poly.printing" begin
   M = MatrixAlgebra(ZZ, 3)
   _, x = M['x']
   @test string(M(-1)*x) isa String
end

@testset "Generic.Poly.polynomial_to_power_sums" begin
   R, x = PolynomialRing(QQ, "x")

   for iters = 1:100
      f = rand(R, 0:10, -10:10)
      f += x^length(f) # make monic
      while iszero(constant_coefficient(f))
         f += rand(-10:10)
      end

      V = polynomial_to_power_sums(f)

      @test power_sums_to_polynomial(V, R) == f
      @test isa(power_sums_to_polynomial(V), PolyElem)

      num = rand(0:2*degree(f))
      W = polynomial_to_power_sums(f, num)
      d = min(num, degree(f))
      @test collect(Base.Iterators.take(V, d)) == collect(Base.Iterators.take(W, d))
   end

   R, x = PolynomialRing(ZZ, "x")

   for iters = 1:100
      f = rand(R, 0:10, -10:10)
      f += x^length(f) # make monic
      while iszero(constant_coefficient(f))
         f += rand(-10:10)
      end

      V = polynomial_to_power_sums(f)

      @test power_sums_to_polynomial(V, R) == f
      @test isa(power_sums_to_polynomial(V), PolyElem)

      num = rand(0:2*degree(f))
      W = polynomial_to_power_sums(f, num)
      d = min(num, degree(f))
      @test collect(Base.Iterators.take(V, d)) == collect(Base.Iterators.take(W, d))
   end
end 
