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

@testset "Generic.PuiseuxSeries.constructors" begin
   R, x = puiseux_series_ring(ZZ, 30, "x")

   S, t = polynomial_ring(QQ, "t")
   T, y = puiseux_series_ring(S, 30, "y")

   @test puiseux_series_ring(S, 30, "y", cached = true)[1] === puiseux_series_ring(S, 30, "y", cached = true)[1]
   @test puiseux_series_ring(S, 30, "y", cached = false)[1] !== puiseux_series_ring(S, 30, "y", cached = true)[1]

   U, z = puiseux_series_field(QQ, 30, "z")

   @test elem_type(R) == Generic.PuiseuxSeriesRingElem{BigInt}
   @test elem_type(Generic.PuiseuxSeriesRing{BigInt}) == Generic.PuiseuxSeriesRingElem{BigInt}
   @test parent_type(Generic.PuiseuxSeriesRingElem{BigInt}) == Generic.PuiseuxSeriesRing{BigInt}

   @test elem_type(T) == Generic.PuiseuxSeriesRingElem{elem_type(S)}
   @test elem_type(Generic.PuiseuxSeriesRing{elem_type(S)}) == Generic.PuiseuxSeriesRingElem{elem_type(S)}
   @test parent_type(Generic.PuiseuxSeriesRingElem{elem_type(S)}) == Generic.PuiseuxSeriesRing{elem_type(S)}

   @test elem_type(U) == Generic.PuiseuxSeriesFieldElem{Rational{BigInt}}
   @test elem_type(Generic.PuiseuxSeriesField{Rational{BigInt}}) == Generic.PuiseuxSeriesFieldElem{Rational{BigInt}}
   @test parent_type(Generic.PuiseuxSeriesFieldElem{Rational{BigInt}}) == Generic.PuiseuxSeriesField{Rational{BigInt}}

   @test isa(R, Generic.PuiseuxSeriesRing)
   @test isa(T, Generic.PuiseuxSeriesRing)
   @test isa(U, Generic.PuiseuxSeriesField)

   a1 = x^3 + 2x + 1
   a2 = (t^2 + 1)*y^2 + (t + 3)*y + O(y^4)
   a3 = z^3 + 2z + 1

   @test isa(a1, Generic.PuiseuxSeriesElem)
   @test isa(a2, Generic.PuiseuxSeriesElem)
   @test isa(a3, Generic.PuiseuxSeriesElem)

   b1 = R(a1)
   b2 = T(a2)
   b3 = U(a3)

   @test isa(b1, Generic.PuiseuxSeriesElem)
   @test isa(b2, Generic.PuiseuxSeriesElem)
   @test isa(b3, Generic.PuiseuxSeriesElem)

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

   @test isa(g1, Generic.PuiseuxSeriesElem)
   @test isa(g2, Generic.PuiseuxSeriesElem)
   @test isa(g3, Generic.PuiseuxSeriesElem)
   @test isa(h1, Generic.PuiseuxSeriesElem)
   @test isa(h2, Generic.PuiseuxSeriesElem)
   @test isa(h3, Generic.PuiseuxSeriesElem)
   @test isa(k1, Generic.PuiseuxSeriesElem)
   @test isa(k2, Generic.PuiseuxSeriesElem)
   @test isa(k3, Generic.PuiseuxSeriesElem)

   l = T(t)

   @test isa(l, Generic.PuiseuxSeriesElem)

   @test x in [x, y]
   @test x in [y, x]
   @test !(x in [y])

   @test x in keys(Dict(x => 1))
   @test !(y in keys(Dict(x => 1)))
end

@testset "Generic.PuiseuxSeries.printing" begin
   R, x = puiseux_series_ring(ZZ, 10, "x")

   @test occursin("O", string(x^(-1//2) + 1 - x + x^2 + x^5))
end

@testset "Generic.PuiseuxSeries.rand" begin
   R, x = puiseux_series_ring(ZZ, 10, "x")

   test_rand(R, -12:12, 1:6, -10:10)
   test_rand(R, -12:12, 1:6, make(ZZ, -10:10))

   R, x = puiseux_series_field(RealField, 10, "x")

   test_rand(R, -12:12, 1:6, -1:1)
   test_rand(R, -12:12, 1:6, make(RealField, -1:1))
end

@testset "Generic.PuiseuxSeries.manipulation" begin
   R, t = polynomial_ring(QQ, "t")
   S, x = puiseux_series_ring(R, 30, "x")

   @test max_precision(S) == 30

   a = 2x + x^3
   b = O(x^4)

   @test valuation(a) == 1
   @test valuation(b) == 4

   @test precision(a) == 31
   @test precision(b) == 4

   a = 2x^(1//3) + x^(2//3) + 3x

   @test coeff(a, 1) == 3
   @test coeff(a, 1//3) == 2

   @test is_gen(gen(S))

   @test iszero(zero(S))

   @test isone(one(S))

   @test is_unit(-1 + x + 2x^2)

   @test isequal(deepcopy(a), a)
   @test isequal(deepcopy(b), b)

   @test characteristic(S) == 0

   T, = residue_ring(ZZ, 7)
   U, y = puiseux_series_ring(T, 10, "y")

   @test modulus(T) == 7

   R, x = puiseux_series_ring(QQ, 10, "x")

   for iter = 1:100
      f = rand(R, -10:10, 1:6, -10:10)

      prec = rand(0:10)//rand(1:10) + 10

      f = set_precision!(f, prec)

      @test precision(f) == prec
   end

   for iter = 1:100
      f = rand(R, -10:10, 1:6, -10:10)

      val = rand(-10:10)//rand(1:10)

      f = set_valuation!(f, val)

      @test valuation(f) == val
   end
end

@testset "Generic.PuiseuxSeries.change_base_ring" begin
   Zx, x = laurent_series_ring(ZZ, 10, "x")
   @test 1 == map_coefficients(sqrt, x^0)
   lp = Zx(BigInt[i for i in 1:10], 10, 11, 5, 1)
   lq = Zx(BigInt[i for i in 10:-1:1], 10, 11, 5, 1)
   Zxx, xx = puiseux_series_ring(ZZ, 10, "x")
   p = Zxx(lp, 1)
   q = Zxx(lp, 1)
   pq = p * q

   for R in [QQ, GF(2), GF(13), ZZ]
      pR = change_base_ring(R, p)
      qR = change_base_ring(R, q, parent = parent(pR))
      @test parent(qR) === parent(pR)
      pqR = change_base_ring(R, pq, parent = parent(pR))
      @test pR * qR == pqR
   end

   ps = map_coefficients(z -> z^2, p)
   @test ps == Zxx(Zx(BigInt[i^2 for i in 1:10], 10, 11, 5, 1), 1)

   f = xx^2 + 3xx^3 + 2xx^6
   @test map_coefficients(one, f) == xx^2 + xx^3 + xx^6
   f2 = map_coefficients(t -> t + 2, f)
   @test f2 == 3xx^2 + 5xx^3 + 4xx^6
   for i in [0, 1, 4, 5]
      @test coeff(f2, i) !== coeff(f, i)
   end

   F = GF(11)
   P, y = puiseux_series_ring(F, 10, "x")
   @test map_coefficients(t -> F(t) + 2, f) == 3y^2 + 5y^3 + 4y^6
end

@testset "Generic.PuiseuxSeries.unary_ops" begin
   #  Exact ring
   R, x = puiseux_series_ring(ZZ, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 1:6, -10:10)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end

   #  Inexact field
   R, x = puiseux_series_field(RealField, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 1:6,-1:1)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = puiseux_series_ring(T, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 1:6, 0:5)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end
end

@testset "Generic.PuiseuxSeries.binary_ops" begin
   #  Exact ring
   R, x = puiseux_series_ring(ZZ, 10, "x")
   for iter = 1:100
      f = rand(R, -12:12, 1:6, -10:10)
      g = rand(R, -12:12, 1:6, -10:10)
      h = rand(R, -12:12, 1:6, -10:10)
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
   R, x = puiseux_series_field(RealField, 10, "x")
   for iter = 1:100
      f = rand(R, -12:12, 1:6, -1:1)
      g = rand(R, -12:12, 1:6, -1:1)
      h = rand(R, -12:12, 1:6, -1:1)
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
   R, x = puiseux_series_ring(T, 10, "x")
   for iter = 1:100
      f = rand(R, -12:12, 1:6, 0:5)
      g = rand(R, -12:12, 1:6, 0:5)
      h = rand(R, -12:12, 1:6, 0:5)
      @test isequal(f + (g + h), (f + g) + h)
      @test isequal(f*g, g*f)
      @test f*(g*h) == (f*g)*h
      @test isequal(f - g, -(g - f))
      @test (f - h) + h == f
      @test f*(g + h) == f*g + f*h
      @test f*(g - h) == f*g - f*h
   end

   # Regression test for Nemo issue 319
   R, x = puiseux_series_ring(ZZ, 4, "x")
   cs = [one(ZZ), zero(ZZ), one(ZZ), zero(ZZ)]
   f = R(R.laurent_ring(cs, 4, 4, 0, 1, false), 2)
   f = deepcopy(f) # triggers rescale
   g = deepcopy(f)

   ff = f*f;

   @test isequal(g, f)
end

@testset "Generic.PuiseuxSeries.adhoc_binary_ops" begin
   # Exact ring
   R, x = puiseux_series_ring(ZZ, 10, "x")
   for iter = 1:500
      f = rand(R, -12:12, 1:6, -10:10)
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
   R, x = puiseux_series_field(RealField, 10, "x")
   for iter = 1:500
      f = rand(R, -12:12, 1:6, -1:1)
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
   S, x = puiseux_series_ring(R, 10, "x")
   for iter = 1:500
      f = rand(S, -12:12, 1:6, 0:5)
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
   S, y = puiseux_series_ring(R, 10, "y")
   for iter = 1:100
      f = rand(S, -12:12, 1:6, 0:5, -10:10)
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

@testset "Generic.PuiseuxSeries.comparison" begin
   # Exact ring
   R, x = puiseux_series_ring(ZZ, 10, "x")
   for iter = 1:500
      f = rand(R, -12:12, 1:6, -10:10)
      g = deepcopy(f)
      h = R()
      while iszero(h)
         h = rand(R, -12:12, 1:6, -10:10)
      end

      @test f == g
      @test isequal(f, g)
      @test (precision(h) > min(precision(f), precision(g)) || f != g + h)
      @test (precision(h) > min(precision(f), precision(g)) || !isequal(f, g + h))
   end

   # Inexact field
   R, x = puiseux_series_field(RealField, 10, "x")
   for iter = 1:500
      f = rand(R, -12:12, 1:6, -1:1)
      g = deepcopy(f)
      h = R()
      while iszero(h)
         h = rand(R, -12:12, 1:6, -1:1)
      end

      @test f == g
      @test isequal(f,  g)
      @test (precision(h) > min(precision(f), precision(g)) || f != g + h)
      @test (precision(h) > min(precision(f), precision(g)) || !isequal(f, g + h))
   end

   # Non-integral domain
   R, = residue_ring(ZZ, 6)
   S, x = puiseux_series_ring(R, 10, "x")
   for iter = 1:500
      f = rand(S, -12:12, 1:6, 0:5)
      g = deepcopy(f)
      h = R()
      while iszero(h)
         h = rand(S, -12:12, 1:6, 0:5)
      end

      @test f == g
      @test isequal(f,  g)
      @test (precision(h) > min(precision(f), precision(g)) || f != g + h)
      @test (precision(h) > min(precision(f), precision(g)) || !isequal(f, g + h))
   end
end

@testset "Generic.PuiseuxSeries.adhoc_comparison" begin
   # Exact ring
   R, x = puiseux_series_ring(ZZ, 10, "x")
   for iter = 1:500
      f = R()
      while f == 0
         f = rand(R, 0:0, 1:6, -10:10)
      end
      f += rand(R, 1:12, 1:6, -10:10)
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
   R, x = puiseux_series_field(RealField, 10, "x")
   for iter = 1:500
      f = R()
      while isapprox(f, R())
         f = rand(R, 0:0, 1:6, -1:1)
      end
      f += rand(R, 1:12, 1:6, -1:1)
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
   S, x = puiseux_series_ring(R, 10, "x")
   for iter = 1:500
      f = S()
      while f == 0
         f = rand(S, 0:0, 1:6, 0:5)
      end
      f += rand(S, 1:12, 1:6, 0:5)
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
   S, y = puiseux_series_ring(R, 10, "y")
   for iter = 1:100
      f = S()
      while f == 0
         f = rand(S, 0:0, 1:6, 0:5, -10:10)
      end
      f += rand(S, 1:12, 1:6, 0:5, -10:10)
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

@testset "Generic.PuiseuxSeries.powering" begin
   # Exact ring
   R, x = puiseux_series_ring(ZZ, 10, "x")

   for iter = 1:100
      f = rand(R, -12:12, 1:6, -10:10)
      r2 = R(1)

      for expn = 0:10
         r1 = f^expn

         @test isequal(r1, r2)

         r2 *= f
      end
   end

   # Inexact field
   R, x = puiseux_series_field(RealField, 10, "x")

   for iter = 1:100
      f = rand(R, -12:12, 1:6, -1:1)
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
      R, x = puiseux_series_ring(Zn, 10, "x")

      f = rand(R, -12:12, 1:6, 0:n - 1)
      r2 = R(1)

      for expn = 0:10
         r1 = f^expn

         @test isequal(r1, r2)

         r2 *= f
      end
   end
end

@testset "Generic.PuiseuxSeries.inversion" begin
   # Exact ring
   R, x = puiseux_series_ring(ZZ, 10, "x")
   for iter = 1:300
      f = R()
      while iszero(f) || !is_unit(coeff(f, valuation(f)))
         f = rand(R, -12:12, 1:6, -10:10)
      end

      @test f*inv(f) == 1
   end

   # Inexact field
   R, x = puiseux_series_field(RealField, 10, "x")
   for iter = 1:300
      f = R()
      while iszero(f)
         f = rand(R, -12:12, 1:6, -1:1)
      end

      @test isapprox(f*inv(f), R(1))
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = puiseux_series_ring(T, 10, "x")
   for iter = 1:300
      f = R()
      while iszero(f) || !is_unit(coeff(f, valuation(f)))
         f = rand(R, -12:12, 1:6, 0:5)
      end

      @test f*inv(f) == 1
   end
end

@testset "Generic.PuiseuxSeries.square_root" begin
   # Exact ring
   R, x = puiseux_series_ring(ZZ, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 1:6, -10:10)
      g = f^2

      @test isequal(sqrt(g)^2, g)

      @test is_square(g)

      f1, s1 = is_square_with_sqrt(g)

      @test f1 && s1^2 == g
   end

   # Inexact field
   R, x = puiseux_series_field(RealField, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 1:6, -1:1)
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

      S, x = puiseux_series_field(R, 10, "x")

      for iter = 1:10
          f = rand(S, -12:12, 1:6, 0:Int(p))

          s = f^2

          @test is_square(s)

          q = sqrt(s)

          @test q^2 == s

          q = sqrt(s; check=false)

          @test q^2 == s

          f1, s1 = is_square_with_sqrt(s)

          @test f1 && s1^2 == s
      end
   end

   R, = residue_field(ZZ, 2)
   T, y = polynomial_ring(R, "x")

   S, x = puiseux_series_ring(T, 10, "x")

   f = x + x^(3//2) + x^2 + x^(5//2) +O(x^6)

   @test is_square(f^2)

   f = x + x^(3//2) + y*x^2 + 7*x^(5//2) +O(x^6)

   @test_throws ErrorException sqrt(f)
end

@testset "Generic.PuiseuxSeries.exact_division" begin
   # Exact ring
   R, x = puiseux_series_ring(ZZ, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 1:6, -10:10)
      g = rand(R, -12:12, 1:6, -10:10)
      while iszero(g) || !is_unit(coeff(g, valuation(g)))
         g = rand(R, -12:12, 1:6, -10:10)
      end

      @test divexact(f, g)*g == f
   end

   # Inexact field
   R, x = puiseux_series_field(RealField, 10, "x")
   for iter = 1:300
      s = rand(0:12)
      f = rand(R, -12:12, 1:6, -1:1)
      g = rand(R, -12:12, 1:6, -1:1)
      while iszero(g) || iszero(coeff(g, valuation(g)))
         g = rand(R, -12:12, 1:6, -1:1)
      end

      @test isapprox(divexact(f, g)*g, f)
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = puiseux_series_ring(T, 10, "x")
   for iter = 1:300
      s = rand(0:12)
      f = rand(R, -12:12, 1:6, 0:5)
      g = rand(R, -12:12, 1:6, 0:5)
      while iszero(g) || !is_unit(coeff(g, valuation(g)))
         g = rand(R, -12:12, 1:6, 0:5)
      end

      @test divexact(f, g)*g == f
   end
end

@testset "Generic.PuiseuxSeries.adhoc_exact_division" begin
   # Exact field
   R, x = puiseux_series_ring(ZZ, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 1:6, -10:10)
      c = ZZ()
      while c == 0
         c = rand(ZZ, -10:10)
      end

      @test isequal(divexact(f*c, c), f)
   end

   # Inexact field
   R, x = puiseux_series_field(RealField, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 1:6, -1:1)
      c = RealField()
      while isapprox(c, 0)
         c = rand(RealField, -1:1)
      end

      @test isapprox(divexact(f*c, c), f)
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = puiseux_series_ring(T, 10, "x")
   for iter = 1:300
      f = rand(R, -12:12, 1:6, 0:5)
      c = T()
      while !is_unit(c)
         c = rand(T, 0:5)
      end

      @test isequal(divexact(f*c, c), f)
   end
end

@testset "Generic.PuiseuxSeries.derivative_integral" begin
   # Exact field
   S, x = puiseux_series_ring(QQ, 10, "x")

   for iter = 1:100
      f = rand(S, -10:10, 1:10, -10:10)

      const_coeff = S(coeff(f, 0))
      const_coeff = set_precision!(const_coeff, precision(f))

      @test isequal(integral(derivative(f)) + const_coeff, f)
   end
 
   # Inexact field
   S, x = puiseux_series_field(RealField, 10, "x")

   for iter = 1:100
      f = rand(S, -10:10, 1:10, -10:10)

      const_coeff = S(coeff(f, 0))
      const_coeff = set_precision!(const_coeff, precision(f))

      @test isapprox(integral(derivative(f)) + const_coeff, f)
   end

   # Non-integral domain
   R, = residue_ring(ZZ, 143)
   S, x = puiseux_series_ring(R, 5, "x")

   for iter = 1:100
      f = rand(S, -5:5, 1:10, -10:10)

      const_coeff = S(coeff(f, 0))
      const_coeff = set_precision!(const_coeff, precision(f))

      @test isequal(integral(derivative(f)) + const_coeff, f)
   end
end

@testset "Generic.PuiseuxSeries.special_functions" begin
   # Exact field
   S, x = puiseux_series_ring(QQ, 10, "x")

   for iter = 1:100
      @test exp(x + O(x^iter)) == exp(x + O(x^(iter - 1)))

      f = S()
      while !is_unit(f)
         f = rand(S, 0:0, 1:6, -10:10)
      end
      g = S()
      while !is_unit(g) || !is_unit(f + g)
         g = rand(S, 0:0, 1:6, -10:10)
      end
      f *= x
      g *= x

      @test isequal(exp(f)*exp(g), exp(f + g))

      @test isequal(log(exp(f)), f)

      while !isone(coeff(f, 0))
         f = rand(S, 0:0, 1:6, -10:10)
      end

      @test isequal(exp(log(f)), f)
   end

   # Inexact field
   S, x = puiseux_series_field(RealField, 10, "x")

   for iter = 1:100
      @test isapprox(exp(x + O(x^iter)), exp(x + O(x^(iter - 1))))

      f = S()
      while coeff(f, 0) == 0
         f = rand(S, 0:0, 1:6, -1:1)
      end
      g = S()
      while coeff(g, 0) == 0 || coeff(f + g, 0) == 0
         g = rand(S, 0:0, 1:6, -1:1)
      end
      f *= x
      g *= x

      @test isapprox(exp(f)*exp(g), exp(f + g))

      @test isapprox(log(exp(f)), f)

      while coeff(f, 0) <= 0
         f = rand(S, 0:0, 1:6, -10:10)
      end

      @test isapprox(exp(log(f)), f)
   end

   # Non-integral domain
   R, = residue_ring(ZZ, 143)
   S, x = puiseux_series_ring(R, 5, "x")

   for iter = 1:10
      @test exp(x + O(x^iter)) == exp(x + O(x^(iter - 1)))

      f = S()
      while !is_unit(coeff(f, 0))
         f = rand(S, 0:0, 1:6, -10:10)
      end
      g = S()
      while !is_unit(coeff(g, 0)) || !is_unit(coeff(f + g, 0))
         g = rand(S, 0:0, 1:6, -10:10)
      end
      f *= x
      g *= x

      @test isequal(exp(f)*exp(g), exp(f + g))
   end
end
