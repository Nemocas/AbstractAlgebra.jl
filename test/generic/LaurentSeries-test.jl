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

@testset "Generic.LaurentSeries.constructors" begin
   R, x = laurent_series_ring(ZZ, 30, "x")

   @test laurent_series_ring(ZZ, 30, "x", cached = true)[1] === laurent_series_ring(ZZ, 30, "x", cached = true)[1]
   @test laurent_series_ring(ZZ, 30, "x", cached = false)[1] !== laurent_series_ring(ZZ, 30, "x", cached = true)[1]

   S, t = polynomial_ring(QQ, "t")
   T, y = laurent_series_ring(S, 30, "y")

   U, z = laurent_series_field(QQ, 30, "z")

   @test elem_type(R) == Generic.LaurentSeriesRingElem{BigInt}
   @test elem_type(Generic.LaurentSeriesRing{BigInt}) == Generic.LaurentSeriesRingElem{BigInt}
   @test parent_type(Generic.LaurentSeriesRingElem{BigInt}) == Generic.LaurentSeriesRing{BigInt}

   @test elem_type(T) == Generic.LaurentSeriesRingElem{elem_type(S)}
   @test elem_type(Generic.LaurentSeriesRing{elem_type(S)}) == Generic.LaurentSeriesRingElem{elem_type(S)}
   @test parent_type(Generic.LaurentSeriesRingElem{elem_type(S)}) == Generic.LaurentSeriesRing{elem_type(S)}

   @test elem_type(U) == Generic.LaurentSeriesFieldElem{Rational{BigInt}}
   @test elem_type(Generic.LaurentSeriesField{Rational{BigInt}}) == Generic.LaurentSeriesFieldElem{Rational{BigInt}}
   @test parent_type(Generic.LaurentSeriesFieldElem{Rational{BigInt}}) == Generic.LaurentSeriesField{Rational{BigInt}}

   @test isa(R, Generic.LaurentSeriesRing)
   @test isa(T, Generic.LaurentSeriesRing)
   @test isa(U, Generic.LaurentSeriesField)

   a1 = x^3 + 2x + 1
   a2 = (t^2 + 1)*y^2 + (t + 3)*y + O(y^4)
   a3 = z^3 + 2z + 1

   @test isa(a1, Generic.LaurentSeriesElem)
   @test isa(a2, Generic.LaurentSeriesElem)
   @test isa(a3, Generic.LaurentSeriesElem)

   b1 = R(a1)
   b2 = T(a2)
   b3 = U(a3)

   @test isa(b1, Generic.LaurentSeriesElem)
   @test isa(b2, Generic.LaurentSeriesElem)
   @test isa(b3, Generic.LaurentSeriesElem)

   c1 = R(BigInt[1, 3, 5], 3, 5, 0, 1)
   c2 = T([t + 1, t, S(1)], 3, 5, 0, 1)
   c3 = U(Rational{BigInt}[1, 3, 5], 3, 5, 0, 1)

   @test isa(c1, Generic.LaurentSeriesElem)
   @test isa(c2, Generic.LaurentSeriesElem)
   @test isa(c3, Generic.LaurentSeriesElem)

   g1 = R(1)
   h1 = R(ZZ(2))
   k1 = R()

   g2 = T(1)
   h2 = T(ZZ(2))
   h2 = T(QQ(2, 3))
   k2 = T()

   g3 = U(1)
   h3 = U(ZZ(2))
   k3 = U()

   @test isa(g1, Generic.LaurentSeriesElem)
   @test isa(g2, Generic.LaurentSeriesElem)
   @test isa(g3, Generic.LaurentSeriesElem)
   @test isa(h1, Generic.LaurentSeriesElem)
   @test isa(h2, Generic.LaurentSeriesElem)
   @test isa(h3, Generic.LaurentSeriesElem)
   @test isa(k1, Generic.LaurentSeriesElem)
   @test isa(k2, Generic.LaurentSeriesElem)
   @test isa(k3, Generic.LaurentSeriesElem)

   l = T(t)

   @test isa(l, Generic.LaurentSeriesElem)

   @test x in [x, y]
   @test x in [y, x]
   @test !(x in [y])

   @test x in keys(Dict(x => 1))
   @test !(y in keys(Dict(x => 1)))
   
   R, x = laurent_series_ring(ZZ, 30, "x")
   RR, x = laurent_series_ring(ZZ, 30, "x")
   @test R === RR
   RR, x = laurent_series_ring(ZZ, 30, "x", cached = false)
   @test R !== RR
end

@testset "Generic.LaurentSeries.rand" begin
   R, x = laurent_series_ring(ZZ, 10, "x")

   test_rand(R, -12:12, -10:10)
   test_rand(R, -12:12, make(ZZ, -10:10))

   R, x = laurent_series_field(RealField, 10, "x")

   test_rand(R, -12:12, -1:1)
   test_rand(R, -12:12, make(RealField, -1:1))
end

@testset "Generic.LaurentSeries.manipulation" begin
   R, t = polynomial_ring(QQ, "t")
   S, x = laurent_series_ring(R, 30, "x")

   @test max_precision(S) == 30

   a = 2x + x^3
   b = O(x^4)

   @test pol_length(a) == 2
   @test pol_length(b) == 0

   @test valuation(a) == 1
   @test valuation(b) == 4

   @test precision(a) == 31
   @test precision(b) == 4

   @test is_gen(gen(S))

   @test iszero(zero(S))

   @test isone(one(S))

   @test is_unit(-1 + x + 2x^2)

   @test isequal(deepcopy(a), a)
   @test isequal(deepcopy(b), b)

   @test normalise(a, 2) == 2

   @test coeff(a, 1) == 2
   @test coeff(b, 7) == 0

   @test_throws DomainError polcoeff(a, -1)
   @test_throws DomainError polcoeff(a, -rand(2:100))

   @test_throws DomainError upscale(a, 0)
   @test_throws DomainError upscale(a, -rand(1:100))
   @test_throws DomainError downscale(a, 0)
   @test_throws DomainError downscale(a, -rand(1:100))

   @test characteristic(S) == 0

   T, = residue_ring(ZZ, 7)
   U, y = laurent_series_ring(T, 10, "y")

   @test modulus(T) == 7
end

@testset "Generic.LaurentSeries.similar" begin
   R, x = laurent_series_ring(ZZ, 10, "x")

   for iters = 1:10
      f = rand(R, 0:10, -10:10)

      g = similar(f, QQ, "y")
      h = similar(f, "y")
      k = similar(f)
      m = similar(f, QQ, 5)
      n = similar(f, 5)

      @test isa(g, Generic.LaurentSeriesFieldElem)
      @test isa(h, Generic.LaurentSeriesRingElem)
      @test isa(k, Generic.LaurentSeriesRingElem)
      @test isa(m, Generic.LaurentSeriesFieldElem)
      @test isa(n, Generic.LaurentSeriesRingElem)

      @test base_ring(g) === QQ
      @test base_ring(m) === QQ

      @test parent(g).S == :y
      @test parent(h).S == :y

      @test iszero(g)
      @test iszero(h)
      @test iszero(k)
      @test iszero(m)
      @test iszero(n)

      @test parent(g) !== parent(f)
      @test parent(h) !== parent(f)
      @test parent(k) === parent(f)
      @test parent(m) !== parent(f)
      @test parent(n) !== parent(f)

      p = similar(f, cached=false)
      q = similar(f, "z", cached=false)
      r = similar(f, "z", cached=false)
      s = similar(f)
      t = similar(f)

      @test parent(p) === parent(f)
      @test parent(q) !== parent(r)
      @test parent(s) === parent(t)

      g = zero(f, QQ, "y")
      h = zero(f, "y")
      k = zero(f)
      m = zero(f, QQ, 5)
      n = zero(f, 5)

      @test isa(g, Generic.LaurentSeriesFieldElem)
      @test isa(h, Generic.LaurentSeriesRingElem)
      @test isa(k, Generic.LaurentSeriesRingElem)
      @test isa(m, Generic.LaurentSeriesFieldElem)
      @test isa(n, Generic.LaurentSeriesRingElem)
   end
end

@testset "Generic.LaurentSeries.change_base_ring" begin
   Zx, x = laurent_series_ring(ZZ, 10, "x")
   @test 1 == map_coefficients(sqrt, x^0)
   p = Zx(BigInt[i for i in 1:10], 10, 11, 5, 1)
   q = Zx(BigInt[i for i in 10:-1:1], 10, 11, 5, 1)
   pq = p * q
   for R in [QQ, GF(2), GF(13), ZZ]
      pR = change_base_ring(R, p)
      qR = change_base_ring(R, q, parent = parent(pR))
      @test parent(qR) === parent(pR)
      pqR = change_base_ring(R, pq, parent = parent(pR))
      @test pR * qR == pqR
   end

   ps = map_coefficients(z -> z^2, p)
   @test ps == Zx(BigInt[i^2 for i in 1:10], 10, 11, 5, 1)

   f = x^2 + 3x^3 + 2x^6
   @test map_coefficients(one, f) == x^2 + x^3 + x^6
   f2 = map_coefficients(t -> t + 2, f)
   @test f2 == 3x^2 + 5x^3 + 4x^6
   for i in [0, 1, 4, 5]
      @test coeff(f2, i) !== coeff(f, i)
   end

   F = GF(11)
   P, y = laurent_series_ring(F, 10, "x")
   @test map_coefficients(t -> F(t) + 2, f) == 3y^2 + 5y^3 + 4y^6
end

@testset "Generic.LaurentSeries.laurent_series" begin
   f = laurent_series(ZZ, [1, 2, 3], 3, 8, 2, 2, "y")

   @test isa(f, Generic.LaurentSeriesRingElem)
   @test base_ring(f) === ZZ
   @test coeff(f, 2) == 1
   @test coeff(f, 4) == 2
   @test f.scale == 2
   @test parent(f).S == :y

   g = laurent_series(ZZ, [1, 2, 3], 3, 7, 4, 1)

   @test isa(g, Generic.LaurentSeriesRingElem)
   @test base_ring(g) === ZZ
   @test coeff(g, 4) == 1
   @test coeff(g, 6) == 3
   @test g.scale == 1
   @test parent(g).S == :x

   h = laurent_series(ZZ, [1, 2, 3], 2, 7, 1, 1)
   k = laurent_series(ZZ, [1, 2, 3], 1, 6, 0, 1, cached=false)
   m = laurent_series(ZZ, [1, 2, 3], 3, 9, 5, 1, cached=false)

   @test parent(h) === parent(g)
   @test parent(k) !== parent(m)

   p = laurent_series(ZZ, BigInt[], 0, 3, 1, 1)
   q = laurent_series(ZZ, [], 0, 3, 2, 1)

   @test isa(p, Generic.LaurentSeriesRingElem)
   @test isa(q, Generic.LaurentSeriesRingElem)

   @test pol_length(p) == 0
   @test pol_length(q) == 0

   r = laurent_series(QQ, BigInt[1, 2, 3], 3, 11, 8, 1)

   @test isa(r, Generic.LaurentSeriesFieldElem)

   s = laurent_series(ZZ, [1, 2, 3], 3, 5, 0, 1; max_precision=10)
   
   @test max_precision(parent(s)) == 10
end

@testset "Generic.LaurentSeries.unary_ops" begin
   #  Exact ring
   R, x = laurent_series_ring(ZZ, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, -10:10)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end

   #  Inexact field
   R, x = laurent_series_field(RealField, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, -1:1)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = laurent_series_ring(T, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 0:5)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end
end

@testset "Generic.LaurentSeries.binary_ops" begin
   #  Exact ring
   R, x = laurent_series_ring(ZZ, 10, "x")
   for iter = 1:100
      f = rand(R, -12:12, -10:10)
      g = rand(R, -12:12, -10:10)
      h = rand(R, -12:12, -10:10)
      @test isequal(f + g, g + f)
      @test isequal(f + (g + h), (f + g) + h)
      @test isequal(f*g, g*f)
      @test isequal(f*(g*h), (f*g)*h)
      @test isequal(f - g, -(g - f))
      @test (f - h) + h == f
      @test f*(g + h) == f*g + f*h
      @test f*(g - h) == f*g - f*h
   end

   #  Inexact field
   R, x = laurent_series_field(RealField, 10, "x")
   for iter = 1:100
      f = rand(R, -12:12, -1:1)
      g = rand(R, -12:12, -1:1)
      h = rand(R, -12:12, -1:1)
      @test isapprox(f + (g + h), (f + g) + h)
      @test isapprox(f*g, g*f)
      @test isapprox(f*(g*h), (f*g)*h)
      @test isapprox(f - g, -(g - f))
      @test isapprox((f - h) + h, f)
      @test isapprox(f*(g + h), f*g + f*h)
      @test isapprox(f*(g - h), f*g - f*h)
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = laurent_series_ring(T, 10, "x")
   for iter = 1:100
      f = rand(R, -12:12, 0:5)
      g = rand(R, -12:12, 0:5)
      h = rand(R, -12:12, 0:5)
      @test isequal(f + (g + h), (f + g) + h)
      @test isequal(f*g, g*f)
      @test f*(g*h) == (f*g)*h
      @test isequal(f - g, -(g - f))
      @test (f - h) + h == f
      @test f*(g + h) == f*g + f*h
      @test f*(g - h) == f*g - f*h
   end

   # Regression test for bug #484
   R, x = laurent_series_ring(QQ, 20,"x")
   a = sum(n*x^n for n in 1:17)
   b = 83 + 43*x^10 + O(x^20)
   c = a + b
   ccp = deepcopy(c)
   addeq!(a, b)
   @test c - a == 0
   @test c - ccp == 0
end

@testset "Generic.LaurentSeries.inplace_binary_ops" begin
   #  Exact ring
   R, x = laurent_series_ring(ZZ, 10, "x")
   for iter = 1:100
      f = rand(R, -12:12, -10:10)
      g = rand(R, -12:12, -10:10)
      h = rand(R, -12:12, -10:10)
      r = R()
      add!(r, f, g)
      @test isequal(r, f + g)
      mul!(r, f, g)
      @test isequal(r, f*g)
   end

   #  Inexact field
   R, x = laurent_series_field(RealField, 10, "x")
   for iter = 1:100
      f = rand(R, -12:12, -1:1)
      g = rand(R, -12:12, -1:1)
      h = rand(R, -12:12, -1:1)
      r = R()
      add!(r, f, g)
      @test isapprox(r, f + g)
      mul!(r, f, g)
      @test isapprox(r, f*g)
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = laurent_series_ring(T, 10, "x")
   for iter = 1:100
      f = rand(R, -12:12, 0:5)
      g = rand(R, -12:12, 0:5)
      h = rand(R, -12:12, 0:5)
      r = R()
      add!(r, f, g)
      @test isequal(r, f + g)
      mul!(r, f, g)
      @test isequal(r, f*g)
   end
end

@testset "Generic.LaurentSeries.adhoc_binary_ops" begin
   # Exact ring
   R, x = laurent_series_ring(ZZ, 10, "x")
   for iter = 1:500
      f = rand(R, -12:12, -10:10)
      c1 = rand(ZZ, -10:10)
      c2 = rand(ZZ, -10:10)
      d1 = rand(zz, -10:10)
      d2 = rand(zz, -10:10)

      @test isequal(c1*f - c2*f, (c1 - c2)*f)
      @test isequal(c1*f + c2*f, (c1 + c2)*f)
      @test isequal(d1*f - d2*f, (d1 - d2)*f)
      @test isequal(d1*f + d2*f, (d1 + d2)*f)

      @test isequal(f*c1 - f*c2, f*(c1 - c2))
      @test isequal(f*c1 + f*c2, f*(c1 + c2))
      @test isequal(f*d1 - f*d2, f*(d1 - d2))
      @test isequal(f*d1 + f*d2, f*(d1 + d2))
   end

   # Inexact field
   R, x = laurent_series_field(RealField, 10, "x")
   for iter = 1:500
      f = rand(R, -12:12, -1:1)
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
   R, = residue_ring(ZZ, 6)
   S, x = laurent_series_ring(R, 10, "x")
   for iter = 1:500
      f = rand(S, -12:12, 0:5)
      c1 = rand(ZZ, -10:10)
      c2 = rand(ZZ, -10:10)
      d1 = rand(zz, -10:10)
      d2 = rand(zz, -10:10)
      a1 = rand(R, 0:5)
      a2 = rand(R, 0:5)

      @test isequal(a1*f - a2*f, (a1 - a2)*f)
      @test isequal(a1*f + a2*f, (a1 + a2)*f)
      @test isequal(c1*f - c2*f, (c1 - c2)*f)
      @test isequal(c1*f + c2*f, (c1 + c2)*f)
      @test isequal(d1*f - d2*f, (d1 - d2)*f)
      @test isequal(d1*f + d2*f, (d1 + d2)*f)

      @test isequal(f*a1 - f*a2, f*(a1 - a2))
      @test isequal(f*a1 + f*a2, f*(a1 + a2))
      @test isequal(f*c1 - f*c2, f*(c1 - c2))
      @test isequal(f*c1 + f*c2, f*(c1 + c2))
      @test isequal(f*d1 - f*d2, f*(d1 - d2))
      @test isequal(f*d1 + f*d2, f*(d1 + d2))
   end

   # Generic tower
   R, x = ZZ["x"]
   S, y = laurent_series_ring(R, 10, "y")
   for iter = 1:100
      f = rand(S, -12:12, 0:5, -10:10)
      c1 = rand(ZZ, -10:10)
      c2 = rand(ZZ, -10:10)
      d1 = rand(R, 0:5, -10:10)
      d2 = rand(R, 0:5, -10:10)

      @test isequal(c1*f - c2*f, (c1 - c2)*f)
      @test isequal(c1*f + c2*f, (c1 + c2)*f)
      @test isequal(d1*f - d2*f, (d1 - d2)*f)
      @test isequal(d1*f + d2*f, (d1 + d2)*f)

      @test isequal(f*c1 - f*c2, f*(c1 - c2))
      @test isequal(f*c1 + f*c2, f*(c1 + c2))
      @test isequal(f*d1 - f*d2, f*(d1 - d2))
      @test isequal(f*d1 + f*d2, f*(d1 + d2))
   end
end

@testset "Generic.LaurentSeries.comparison" begin
   # Exact ring
   R, x = laurent_series_ring(ZZ, 10, "x")
   for iter = 1:500
      f = rand(R, -12:12, -10:10)
      g = deepcopy(f)
      h = R()
      while iszero(h)
         h = rand(R, -12:12, -10:10)
      end

      @test f == g
      @test isequal(f, g)
      @test (precision(h) > min(precision(f), precision(g)) || f != g + h)
      @test (precision(h) > min(precision(f), precision(g)) || !isequal(f, g + h))
   end

   # Inexact field
   R, x = laurent_series_field(RealField, 10, "x")
   for iter = 1:500
      f = rand(R, -12:12, -1:1)
      g = deepcopy(f)
      h = R()
      while iszero(h)
         h = rand(R, -12:12, -1:1)
      end

      @test f == g
      @test isequal(f,  g)
      @test (precision(h) > min(precision(f), precision(g)) || f != g + h)
      @test (precision(h) > min(precision(f), precision(g)) || !isequal(f, g + h))
   end

   # Non-integral domain
   R, = residue_ring(ZZ, 6)
   S, x = laurent_series_ring(R, 10, "x")
   for iter = 1:500
      f = rand(S, -12:12, 0:5)
      g = deepcopy(f)
      h = R()
      while iszero(h)
         h = rand(S, -12:12, 0:5)
      end

      @test f == g
      @test isequal(f,  g)
      @test (precision(h) > min(precision(f), precision(g)) || f != g + h)
      @test (precision(h) > min(precision(f), precision(g)) || !isequal(f, g + h))
   end
end

@testset "Generic.LaurentSeries.adhoc_comparison" begin
   # Exact ring
   R, x = laurent_series_ring(ZZ, 10, "x")
   for iter = 1:500
      f = R()
      while f == 0
         f = rand(R, 0:0, -10:10)
      end
      f += rand(R, 1:12, -10:10)
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
   R, x = laurent_series_field(RealField, 10, "x")
   for iter = 1:500
      f = R()
      while isapprox(f, R())
         f = rand(R, 0:0, -1:1)
      end
      f += rand(R, 1:12, -1:1)
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
   R, = residue_ring(ZZ, 6)
   S, x = laurent_series_ring(R, 10, "x")
   for iter = 1:500
      f = S()
      while f == 0
         f = rand(S, 0:0, 0:5)
      end
      f += rand(S, 1:12, 0:5)
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
   S, y = laurent_series_ring(R, 10, "y")
   for iter = 1:100
      f = S()
      while f == 0
         f = rand(S, 0:0, 0:5, -10:10)
      end
      f += rand(S, 1:12, 0:5, -10:10)
      c1 = rand(ZZ, -10:10)
      d1 = rand(R, 0:5, -10:10)

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

@testset "Generic.LaurentSeries.powering" begin
   # Exact ring
   R, x = laurent_series_ring(ZZ, 10, "x")

   for iter = 1:100
      f = rand(R, -12:12, -10:10)
      r2 = R(1)

      for expn = 0:10
         r1 = f^expn

         @test isequal(r1, r2)

         r2 *= f
      end
   end

   # Inexact field
   R, x = laurent_series_field(RealField, 10, "x")

   for iter = 1:100
      f = rand(R, -12:12, -1:1)
      r2 = R(1)

      for expn = 0:4 # cannot set high power here
         r1 = f^expn

         @test isapprox(r1, r2)

         r2 *= f
      end
   end

   # Non-integral domain
   for iter = 1:100
      n = rand(2:26)

      Zn, = residue_ring(ZZ, n)
      R, x = laurent_series_ring(Zn, 10, "x")

      f = rand(R, -12:12, 0:n - 1)
      r2 = R(1)

      for expn = 0:10
         r1 = f^expn

         @test isequal(r1, r2)

         r2 *= f
      end
   end

   # regression test (see #967)
   Zn, = residue_ring(ZZ, 4)
   R, x = laurent_series_ring(Zn, 5, "x")
   f = 2*x^6 + O(x^11)

   @test isequal(f*f, f^2)
   @test isequal(f^0, one(R))
end

@testset "Generic.LaurentSeries.shift" begin
   # Exact ring
   R, x = laurent_series_ring(ZZ, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, -10:10)
      s = rand(0:12)

      @test isequal(shift_right(shift_left(f, s), s), f)
      @test isequal(shift_left(f, s), x^s*f)
      @test precision(shift_right(f, s)) == precision(f) - s
   end

   # Inexact field
   R, x = laurent_series_field(RealField, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, -1:1)
      s = rand(0:12)

      @test isapprox(shift_right(shift_left(f, s), s), f)
      @test isapprox(shift_left(f, s), x^s*f)
      @test precision(shift_right(f, s)) == precision(f) - s
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = laurent_series_ring(T, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 0:5)
      s = rand(0:12)

      @test isequal(shift_right(shift_left(f, s), s), f)
      @test isequal(shift_left(f, s), x^s*f)
      @test precision(shift_right(f, s)) == precision(f) - s
   end
end

@testset "Generic.LaurentSeries.truncation" begin
   # Exact ring
   R, x = laurent_series_ring(ZZ, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, -10:10)
      s = rand(-12:12)

      @test truncate(f, s) == f
      @test isequal(truncate(f, s), f + O(x^s))
      @test precision(truncate(f, s)) == min(precision(f), s)
   end

   # Inexact field
   R, x = laurent_series_field(RealField, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, -1:1)
      s = rand(-12:12)

      @test truncate(f, s) == f
      @test isapprox(truncate(f, s), f + O(x^s))
      @test precision(truncate(f, s)) == min(precision(f), s)
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = laurent_series_ring(T, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 0:5)
      s = rand(-12:12)

      @test truncate(f, s) == f
      @test isequal(truncate(f, s), f + O(x^s))
      @test precision(truncate(f, s)) == min(precision(f), s)
   end
end

@testset "Generic.LaurentSeries.inversion" begin
   # Exact ring
   R, x = laurent_series_ring(ZZ, 10, "x")
   for iter = 1:300
      f = R()
      while iszero(f) || !is_unit(polcoeff(f, 0))
         f = rand(R, -12:12, -10:10)
      end

      @test f*inv(f) == 1
   end

   # Inexact field
   R, x = laurent_series_field(RealField, 10, "x")
   for iter = 1:300
      f = R()
      while iszero(f)
         f = rand(R, -12:12, -1:1)
      end

      @test isapprox(f*inv(f), R(1))
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = laurent_series_ring(T, 10, "x")
   for iter = 1:300
      f = R()
      while iszero(f) || !is_unit(polcoeff(f, 0))
         f = rand(R, -12:12, 0:5)
      end

      @test f*inv(f) == 1
   end
end

@testset "Generic.LaurentSeries.square_root" begin
   # Exact ring
   R, x = laurent_series_ring(ZZ, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, -10:10)
      g = f^2

      @test isequal(sqrt(g)^2, g)

      @test is_square(g)

      f1, s1 = is_square_with_sqrt(g)

      @test f1 && s1^2 == g

      if g*x != 0
          @test_throws ErrorException sqrt(g*x)
      end
   end

   f = 4 + 12*x + 17*x^2 + 12*x^3 + 25*x^4 + 30*x^5 + 20*x^6 + 25*x^8 + O(x^10)

   @test_throws ErrorException sqrt(f)

   f = 7 + 12*x + 17*x^2 + 12*x^3 + 25*x^4 + 30*x^5 + 20*x^6 + 25*x^8 + O(x^10)

   @test_throws ErrorException sqrt(f)

   # Inexact field
   R, x = laurent_series_field(RealField, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, -1:1)
      g = f^2

      @test isapprox(sqrt(g)^2, g)
   end

   f = 0.4684049219189188040246563105029053986072540283203125*x +
        0.7468877608885249319570220905006863176822662353515625*x^2 +
        0.70502205920556537677157393773086369037628173828125*x^3 +
        0.0338681672304976810750076765543781220912933349609375*x^4 +
        0.520844397208659604103786477935500442981719970703125*x^5 +
        0.3137863950850567906769583714776672422885894775390625*x^6 +
        0.593894760292587253758256338187493383884429931640625*x^7 +
        0.62582529337203851582671632058918476104736328125*x^8 +
        0.8365754395386633124331865474232472479343414306640625*x^9 +
        O(x^10)

   flag, s = is_square_with_sqrt(f^2)

   @test flag && isapprox(s, f)

   # Characteristic p field
   for p in [2, 7, 19, 65537, ZZ(7), ZZ(19), ZZ(65537)]
      R, = residue_field(ZZ, p)

      S, x = laurent_series_field(R, 10, "x")

      for iter = 1:10
          f = rand(S, 0:10, 0:Int(p))

          s = f^2

          @test is_square(s)

          q = sqrt(s)

          @test q^2 == s

          q = sqrt(s; check=false)

          @test q^2 == s

          f1, s1 = is_square_with_sqrt(s)

          @test f1 && s1^2 == s

          if s*x != 0
              @test_throws ErrorException sqrt(s*x)
          end
      end
   end

    R, = residue_field(ZZ, 2)
    T, y = polynomial_ring(R, "x")

    S, x = laurent_series_ring(T, 10, "x")

    f = 1 + y^2*x^2 + (y^2 + y + 1)*x^4 + O(x^10)

    @test_throws ErrorException sqrt(f)

    f = x^2 + x^3

    @test_throws ErrorException sqrt(f)

    f = x + y^4*x^4 + (y^4 + y^2 + 1)*x^8 + O(x^11)

    @test_throws ErrorException sqrt(f)

    f = zero(S)

    @test iszero(sqrt(f))
end

@testset "Generic.LaurentSeries.exact_division" begin
   # Exact ring
   R, x = laurent_series_ring(ZZ, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, -10:10)
      g = rand(R, -12:12, -10:10)
      while iszero(g) || !is_unit(polcoeff(g, 0))
         g = rand(R, -12:12, -10:10)
      end

      @test divexact(f, g)*g == f
   end

   # Inexact field
   R, x = laurent_series_field(RealField, 10, "x")
   for iter = 1:300
      s = rand(0:12)
      f = rand(R, -12:12, -1:1)
      g = rand(R, -12:12, -1:1)
      while iszero(g) || iszero(polcoeff(g, 0))
         g = rand(R, -12:12, -1:1)
      end

      @test isapprox(divexact(f, g)*g, f)
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = laurent_series_ring(T, 3, "x")
   for iter = 1:3000
      f = rand(R, -12:12, 0:5)
      g = rand(R, -12:12, 0:5)
      while iszero(g) || !is_unit(polcoeff(g, 0))
         g = rand(R, -12:12, 0:5)
      end

      @test divexact(f, g)*g == f
   end
end

@testset "Generic.LaurentSeries.adhoc_exact_division" begin
   # Exact field
   R, x = laurent_series_ring(ZZ, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, -10:10)
      c = ZZ()
      while c == 0
         c = rand(ZZ, -10:10)
      end

      @test isequal(divexact(f*c, c), f)
   end

   # Inexact field
   R, x = laurent_series_field(RealField, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, -1:1)
      c = RealField()
      while isapprox(c, 0)
         c = rand(RealField, -1:1)
      end

      @test isapprox(divexact(f*c, c), f)
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = laurent_series_ring(T, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 0:5)
      c = T()
      while !is_unit(c)
         c = rand(T, 0:5)
      end

      @test isequal(divexact(f*c, c), f)
   end
end

@testset "Generic.LaurentSeries.derivative_integral" begin
   # Exact field
   S, x = laurent_series_ring(QQ, 10, "x")

   for iter = 1:100
      f = rand(S, -10:10, -10:10)

      const_coeff = S(coeff(f, 0))
      const_coeff = set_precision!(const_coeff, precision(f))

      @test isequal(integral(derivative(f)) + const_coeff, f)
   end

   # Inexact field
   S, x = laurent_series_ring(RealField, 10, "x")

   for iter = 1:100
      f = rand(S, -10:10, -10:10)

      const_coeff = S(coeff(f, 0))
      set_precision!(const_coeff, precision(f))

      @test isapprox(integral(derivative(f)) + const_coeff, f)
   end
end

@testset "Generic.LaurentSeries.special_functions" begin
   # Exact field
   S, x = laurent_series_ring(QQ, 10, "x")

   for iter = 1:100
      @test exp(x + O(x^iter)) == exp(x + O(x^(iter - 1)))

      f = S()
      while !is_unit(f)
         f = rand(S, 0:0, -10:10)
      end
      g = S()
      while !is_unit(g) || !is_unit(f + g)
         g = rand(S, 0:0, -10:10)
      end
      f *= x
      g *= x

      @test isequal(exp(f)*exp(g), exp(f + g))

      @test isequal(log(exp(f)), f)

      while !isone(coeff(f, 0))
         f = rand(S, 0:0, -10:10)
      end

      @test isequal(exp(log(f)), f)
   end

   # Inexact field
   S, x = laurent_series_field(RealField, 10, "x")

   for iter = 1:100
      @test isapprox(exp(x + O(x^iter)), exp(x + O(x^(iter - 1))))

      f = S()
      while coeff(f, 0) == 0
         f = rand(S, 0:0, -1:1)
      end
      g = S()
      while coeff(g, 0) == 0 || coeff(f + g, 0) == 0
         g = rand(S, 0:0, -1:1)
      end
      f *= x
      g *= x

      while coeff(f, 0) <= 0
         f = rand(S, 0:0, -1:1)
      end

      @test isapprox(exp(f)*exp(g), exp(f + g))

      @test isapprox(log(exp(f)), f)

      while coeff(f, 0) <= 0
         f = rand(S, 0:0, -10:10)
      end

      @test isapprox(exp(log(f)), f)
   end

   # Non-integral domain
   R, = residue_ring(ZZ, 143)
   S, x = laurent_series_ring(R, 10, "x")

   for iter = 1:10
      @test exp(x + O(x^iter)) == exp(x + O(x^(iter - 1)))

      f = S()
      while !is_unit(coeff(f, 0))
         f = rand(S, 0:0, -10:10)
      end
      g = S()
      while !is_unit(coeff(g, 0)) || !is_unit(coeff(f + g, 0))
         g = rand(S, 0:0, -10:10)
      end
      f *= x
      g *= x

      @test isequal(exp(f)*exp(g), exp(f + g))
   end
end
