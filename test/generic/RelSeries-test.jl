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

function test_rel_series_constructors()
   print("Generic.RelSeries.constructors...")

   R, x = PowerSeriesRing(JuliaZZ, 30, "x")

   S, t = PolynomialRing(JuliaQQ, "t")
   T, y = PowerSeriesRing(S, 30, "y")

   @test elem_type(R) == Generic.RelSeries{BigInt}
   @test elem_type(Generic.RelSeriesRing{BigInt}) == Generic.RelSeries{BigInt}
   @test parent_type(Generic.RelSeries{BigInt}) == Generic.RelSeriesRing{BigInt}

   @test elem_type(T) == Generic.RelSeries{elem_type(S)}
   @test elem_type(Generic.RelSeriesRing{elem_type(S)}) == Generic.RelSeries{elem_type(S)}
   @test parent_type(Generic.RelSeries{elem_type(S)}) == Generic.RelSeriesRing{elem_type(S)}

   @test isa(R, Generic.RelSeriesRing)

   @test isa(T, Generic.RelSeriesRing)

   a1 = x^3 + 2x + 1
   a2 = (t^2 + 1)*y^2 + (t + 3)*y + O(y^4)

   @test isa(a1, Generic.RelSeries)
   @test isa(a2, Generic.RelSeries)

   b1 = R(a1)
   b2 = T(a2)

   @test isa(b1, Generic.RelSeries)
   @test isa(b2, Generic.RelSeries)

   c1 = R(BigInt[1, 3, 5], 3, 5, 0)
   c2 = T([t + 1, t, S(1)], 3, 5, 0)

   @test isa(c1, Generic.RelSeries)
   @test isa(c2, Generic.RelSeries)

   g1 = R(1)
   h1 = R(JuliaZZ(2))
   k1 = R()

   g2 = T(1)
   h2 = T(JuliaZZ(2))
   h2 = T(JuliaQQ(2, 3))
   k2 = T()

   @test isa(g1, Generic.RelSeries)
   @test isa(g2, Generic.RelSeries)
   @test isa(h1, Generic.RelSeries)
   @test isa(h2, Generic.RelSeries)
   @test isa(k1, Generic.RelSeries)
   @test isa(k2, Generic.RelSeries)

   l = T(t)

   @test isa(l, Generic.RelSeries)

   println("PASS")
end

function test_rel_series_manipulation()
   print("Generic.RelSeries.manipulation...")

   R, t = PolynomialRing(JuliaQQ, "t")
   S, x = PowerSeriesRing(R, 30, "x")

   @test max_precision(S) == 30

   a = 2x + x^3
   b = O(x^4)

   @test pol_length(a) == 3
   @test pol_length(b) == 0

   @test valuation(a) == 1
   @test valuation(b) == 4

   @test precision(a) == 31
   @test precision(b) == 4

   @test isgen(gen(S))

   @test iszero(zero(S))

   @test isone(one(S))

   @test isunit(-1 + x + 2x^2)

   @test isequal(deepcopy(a), a)
   @test isequal(deepcopy(b), b)

   @test normalise(a, 3) == 3

   @test coeff(a, 1) == 2
   @test coeff(b, 7) == 0

   T = ResidueRing(JuliaZZ, 7)
   U, y = PowerSeriesRing(T, 10, "y")

   @test modulus(T) == 7

   println("PASS")
end

function test_rel_series_unary_ops()
   print("Generic.RelSeries.unary_ops...")

   #  Exact ring
   R, x = PowerSeriesRing(JuliaZZ, 10, "x")
   for iter = 1:300
      f = rand(R, 0:12, -10:10)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end

   #  Inexact field
   R, x = PowerSeriesRing(JuliaRealField, 10, "x")
   for iter = 1:300
      f = rand(R, 0:12, -1:1)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end

   # Non-integral domain
   T = ResidueRing(JuliaZZ, 6)
   R, x = PowerSeriesRing(T, 10, "x")
   for iter = 1:300
      f = rand(R, 0:12, 0:5)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end

   println("PASS")
end

function test_rel_series_binary_ops()
   print("Generic.RelSeries.binary_ops...")

   #  Exact ring
   R, x = PowerSeriesRing(JuliaZZ, 10, "x")
   for iter = 1:100
      f = rand(R, 0:12, -10:10)
      g = rand(R, 0:12, -10:10)
      h = rand(R, 0:12, -10:10)
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
   R, x = PowerSeriesRing(JuliaRealField, 10, "x")
   for iter = 1:100
      f = rand(R, 0:12, -1:1)
      g = rand(R, 0:12, -1:1)
      h = rand(R, 0:12, -1:1)
      @test isapprox(f + (g + h), (f + g) + h)
      @test isapprox(f*g, g*f)
      @test isapprox(f*(g*h), (f*g)*h)
      @test isapprox(f - g, -(g - f))
      @test isapprox((f - h) + h, f)
      @test isapprox(f*(g + h), f*g + f*h)
      @test isapprox(f*(g - h), f*g - f*h)
   end

   # Non-integral domain
   T = ResidueRing(JuliaZZ, 6)
   R, x = PowerSeriesRing(T, 10, "x")
   for iter = 1:100
      f = rand(R, 0:12, 0:5)
      g = rand(R, 0:12, 0:5)
      h = rand(R, 0:12, 0:5)
      @test isequal(f + (g + h), (f + g) + h)
      @test isequal(f*g, g*f)
      @test f*(g*h) == (f*g)*h
      @test isequal(f - g, -(g - f))
      @test (f - h) + h == f
      @test f*(g + h) == f*g + f*h
      @test f*(g - h) == f*g - f*h
   end

   println("PASS")
end

function test_rel_series_adhoc_binary_ops()
   print("Generic.RelSeries.adhoc_binary_ops...")

   # Exact ring
   R, x = PowerSeriesRing(JuliaZZ, 10, "x")
   for iter = 1:500
      f = rand(R, 0:12, -10:10)
      c1 = rand(JuliaZZ, -10:10)
      c2 = rand(JuliaZZ, -10:10)
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
   R, x = PowerSeriesRing(JuliaRealField, 10, "x")
   for iter = 1:500
      f = rand(R, 0:12, -1:1)
      c1 = rand(JuliaZZ, -10:10)
      c2 = rand(JuliaZZ, -10:10)
      d1 = rand(JuliaRealField, -1:1)
      d2 = rand(JuliaRealField, -1:1)

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
   R = ResidueRing(JuliaZZ, 6)
   S, x = PowerSeriesRing(R, 10, "x")
   for iter = 1:500
      f = rand(S, 0:12, 0:5)
      c1 = rand(JuliaZZ, -10:10)
      c2 = rand(JuliaZZ, -10:10)
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
   R, x = JuliaZZ["x"]
   S, y = PowerSeriesRing(R, 10, "y")
   for iter = 1:100
      f = rand(S, 0:12, 0:5, -10:10)
      c1 = rand(JuliaZZ, -10:10)
      c2 = rand(JuliaZZ, -10:10)
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

   println("PASS")
end

function test_rel_series_comparison()
   print("Generic.RelSeries.comparison...")

   # Exact ring
   R, x = PowerSeriesRing(JuliaZZ, 10, "x")
   for iter = 1:500
      f = rand(R, 0:12, -10:10)
      g = deepcopy(f)
      h = R()
      while iszero(h)
         h = rand(R, 0:12, -10:10)
      end

      @test f == g
      @test isequal(f, g)
      @test (precision(h) > min(precision(f), precision(g)) || f != g + h)
      @test (precision(h) > min(precision(f), precision(g)) || !isequal(f, g + h))
   end

   # Inexact field
   R, x = PowerSeriesRing(JuliaRealField, 10, "x")
   for iter = 1:500
      f = rand(R, 0:12, -1:1)
      g = deepcopy(f)
      h = R()
      while iszero(h)
         h = rand(R, 0:12, -1:1)
      end

      @test f == g
      @test isequal(f,  g)
      @test (precision(h) > min(precision(f), precision(g)) || f != g + h)
      @test (precision(h) > min(precision(f), precision(g)) || !isequal(f, g + h))
   end

   # Non-integral domain
   R = ResidueRing(JuliaZZ, 6)
   S, x = PowerSeriesRing(R, 10, "x")
   for iter = 1:500
      f = rand(S, 0:12, 0:5)
      g = deepcopy(f)
      h = R()
      while iszero(h)
         h = rand(S, 0:12, 0:5)
      end

      @test f == g
      @test isequal(f,  g)
      @test (precision(h) > min(precision(f), precision(g)) || f != g + h)
      @test (precision(h) > min(precision(f), precision(g)) || !isequal(f, g + h))
   end

   println("PASS")
end

function test_rel_series_adhoc_comparison()
   print("Generic.RelSeries.adhoc_comparison...")

   # Exact ring
   R, x = PowerSeriesRing(JuliaZZ, 10, "x")
   for iter = 1:500
      f = R()
      while f == 0
         f = rand(R, 0:0, -10:10)
      end
      f += rand(R, 1:12, -10:10)
      c1 = rand(JuliaZZ, -10:10)
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
   R, x = PowerSeriesRing(JuliaRealField, 10, "x")
   for iter = 1:500
      f = R()
      while isapprox(f, R())
         f = rand(R, 0:0, -1:1)
      end
      f += rand(R, 1:12, -1:1)
      c1 = rand(JuliaZZ, -10:10)
      d1 = rand(JuliaRealField, -1:1)

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
   R = ResidueRing(JuliaZZ, 6)
   S, x = PowerSeriesRing(R, 10, "x")
   for iter = 1:500
      f = S()
      while f == 0
         f = rand(S, 0:0, 0:5)
      end
      f += rand(S, 1:12, 0:5)
      c1 = rand(JuliaZZ, -10:10)
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
   R, x = JuliaZZ["x"]
   S, y = PowerSeriesRing(R, 10, "y")
   for iter = 1:100
      f = S()
      while f == 0
         f = rand(S, 0:0, 0:5, -10:10)
      end
      f += rand(S, 1:12, 0:5, -10:10)
      c1 = rand(JuliaZZ, -10:10)
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

   println("PASS")
end

function test_rel_series_powering()
   print("Generic.RelSeries.powering...")

   R, t = PolynomialRing(QQ, "t")
   S, x = PowerSeriesRing(R, 30, "x")

   a = 2x + x^3
   b = O(x^4)
   c = 1 + x + 2x^2 + O(x^5)
   d = 2x + x^3 + O(x^4)

   @test isequal(a^12, x^36+24*x^34+264*x^32+1760*x^30+7920*x^28+25344*x^26+59136*x^24+101376*x^22+126720*x^20+112640*x^18+67584*x^16+24576*x^14+4096*x^12+O(x^42))

   @test isequal(b^12, O(x^48))

   @test isequal(c^12, 2079*x^4+484*x^3+90*x^2+12*x+1+O(x^5))

   @test isequal(d^12, 4096*x^12+24576*x^14+O(x^15))

   println("PASS")
end

function test_rel_series_shift()
   print("Generic.RelSeries.shift...")

   R, t = PolynomialRing(QQ, "t")
   S, x = PowerSeriesRing(R, 30, "x")

   a = 2x + x^3
   b = O(x^4)
   c = 1 + x + 2x^2 + O(x^5)
   d = 2x + x^3 + O(x^4)

   @test isequal(shift_left(a, 2), 2*x^3+x^5+O(x^33))

   @test isequal(shift_left(b, 2), O(x^6))

   @test isequal(shift_right(c, 1), 1+2*x+O(x^4))

   @test isequal(shift_right(d, 3), 1+O(x^1))

   println("PASS")
end

function test_rel_series_truncation()
   print("Generic.RelSeries.truncation...")

   R, t = PolynomialRing(QQ, "t")
   S, x = PowerSeriesRing(R, 30, "x")

   a = 2x + x^3
   b = O(x^4)
   c = 1 + x + 2x^2 + O(x^5)
   d = 2x + x^3 + O(x^4)

   @test isequal(truncate(a, 3), 2*x + O(x^3))

   @test isequal(truncate(b, 2), O(x^2))

   @test isequal(truncate(c, 5), 2*x^2+x+1+O(x^5))

   @test isequal(truncate(d, 5), x^3+2*x+O(x^4))

   println("PASS")
end

function test_rel_series_inversion()
   print("Generic.RelSeries.inversion...")

   R, t = PolynomialRing(QQ, "t")
   S, x = PowerSeriesRing(R, 30, "x")

   a = 1 + x + 2x^2 + O(x^5)
   b = S(-1)

   @test isequal(inv(a), -x^4+3*x^3-x^2-x+1+O(x^5))

   @test isequal(inv(b), -1+O(x^30))

   println("PASS")
end

function test_rel_series_exact_division()
   print("Generic.RelSeries.exact_division...")

   R, t = PolynomialRing(QQ, "t")
   S, x = PowerSeriesRing(R, 30, "x")

   a = x + x^3
   b = O(x^4)
   c = 1 + x + 2x^2 + O(x^5)
   d = x + x^3 + O(x^6)

   @test isequal(divexact(a, d), 1+O(x^5))

   @test isequal(divexact(d, a), 1+O(x^5))

   @test isequal(divexact(b, c), O(x^4))

   @test isequal(divexact(d, c), -2*x^5+2*x^4-x^2+x+O(x^6))

   println("PASS")
end

function test_rel_series_adhoc_exact_division()
   print("Generic.RelSeries.adhoc_exact_division...")

   R, t = PolynomialRing(QQ, "t")
   S, x = PowerSeriesRing(R, 30, "x")

   a = x + x^3
   b = O(x^4)
   c = 1 + x + 2x^2 + O(x^5)
   d = x + x^3 + O(x^6)

   @test isequal(divexact(a, 7), fmpz(1)//7*x+fmpz(1)//7*x^3+O(x^31))

   @test isequal(divexact(b, fmpz(11)), 0+O(x^4))

   @test isequal(divexact(c, fmpz(2)), fmpz(1)//2+fmpz(1)//2*x+x^2+O(x^5))

   @test isequal(divexact(d, 9), fmpz(1)//9*x+fmpz(1)//9*x^3+O(x^6))

   @test isequal(divexact(94872394861923874346987123694871329847a, 94872394861923874346987123694871329847), a)

   @test isequal(divexact((t + 1)*a, t + 1), a)

   println("PASS")
end

function test_rel_series_special_functions()
   print("Generic.RelSeries.special_functions...")

   R = ResidueRing(ZZ, 17)
   T, t = PolynomialRing(R, "t")
   S, x = PowerSeriesRing(T, 30, "x")

   @test isequal(exp(x + O(x^10)),  8*x^9+4*x^8+15*x^7+3*x^6+x^5+5*x^4+3*x^3+9*x^2+x+1+O(x^10))

   @test isequal(divexact(x, exp(x + O(x^10)) - 1), x^8+11*x^6+14*x^4+10*x^2+8*x+1+O(x^9))

   println("PASS")
end

function test_gen_rel_series()
   test_rel_series_constructors()
   test_rel_series_manipulation()
   test_rel_series_unary_ops()
   test_rel_series_binary_ops()
   test_rel_series_adhoc_binary_ops()
   test_rel_series_comparison()
   test_rel_series_adhoc_comparison()
   test_rel_series_powering()
   test_rel_series_shift()
   test_rel_series_truncation()
   test_rel_series_exact_division()
   test_rel_series_adhoc_exact_division()
   test_rel_series_inversion()
   test_rel_series_special_functions()

   println("")
end
