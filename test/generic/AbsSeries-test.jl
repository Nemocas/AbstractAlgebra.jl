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

@testset "Generic.AbsSeries.constructors..." begin
   R, x = PowerSeriesRing(ZZ, 30, "x", model=:capped_absolute)

   S, t = PolynomialRing(QQ, "t")
   T, y = PowerSeriesRing(S, 30, "y", model=:capped_absolute)

   @test elem_type(R) == Generic.AbsSeries{BigInt}
   @test elem_type(Generic.AbsSeriesRing{BigInt}) == Generic.AbsSeries{BigInt}
   @test parent_type(Generic.AbsSeries{BigInt}) == Generic.AbsSeriesRing{BigInt}

   @test elem_type(T) == Generic.AbsSeries{elem_type(S)}
   @test elem_type(Generic.AbsSeriesRing{elem_type(S)}) == Generic.AbsSeries{elem_type(S)}
   @test parent_type(Generic.AbsSeries{elem_type(S)}) == Generic.AbsSeriesRing{elem_type(S)}

   @test isa(R, Generic.AbsSeriesRing)

   @test isa(T, Generic.AbsSeriesRing)

   a1 = x^3 + 2x + 1
   a2 = (t^2 + 1)*y^2 + (t + 3)*y + O(y^4)

   @test isa(a1, Generic.AbsSeries)
   @test isa(a2, Generic.AbsSeries)

   b1 = R(a1)
   b2 = T(a2)

   @test isa(b1, Generic.AbsSeries)
   @test isa(b2, Generic.AbsSeries)

   c1 = R(BigInt[1, 3, 5], 3, 5)
   c2 = T([t + 1, t, S(1)], 3, 5)

   @test isa(c1, Generic.AbsSeries)
   @test isa(c2, Generic.AbsSeries)

   g1 = R(1)
   h1 = R(ZZ(2))
   k1 = R()

   g2 = T(1)
   h2 = T(ZZ(2))
   h2 = T(QQ(2, 3))
   k2 = T()

   @test isa(g1, Generic.AbsSeries)
   @test isa(g2, Generic.AbsSeries)
   @test isa(h1, Generic.AbsSeries)
   @test isa(h2, Generic.AbsSeries)
   @test isa(k1, Generic.AbsSeries)
   @test isa(k2, Generic.AbsSeries)

   l = T(t)

   @test isa(l, Generic.AbsSeries)

   @test x in [x, y]
   @test x in [y, x]
   @test !(x in [y])

   @test x in keys(Dict(x => 1))
   @test !(y in keys(Dict(x => 1)))
end

@testset "Generic.AbsSeries.manipulation..." begin
   R, t = PolynomialRing(QQ, "t")
   S, x = PowerSeriesRing(R, 30, "x", model=:capped_absolute)

   @test max_precision(S) == 30

   a = 2x + x^3
   b = O(x^4)

   @test valuation(a) == 1
   @test valuation(b) == 4

   @test precision(a) == 30
   @test precision(b) == 4

   @test isgen(gen(S))

   @test iszero(zero(S))

   @test isone(one(S))

   @test isunit(-1 + x + 2x^2)

   @test isequal(deepcopy(a), a)
   @test isequal(deepcopy(b), b)

   @test normalise(a, 4) == 4

   @test coeff(a, 1) == 2
   @test coeff(b, 7) == 0

   @test_throws DomainError coeff(a, -1)
   @test_throws DomainError coeff(a, -rand(2:100))
   
   @test characteristic(R) == 0

   T = ResidueRing(ZZ, 7)
   U, y = PowerSeriesRing(T, 10, "y", model=:capped_absolute)

   @test modulus(T) == 7
end

@testset "Generic.AbsSeries.unary_ops..." begin
   #  Exact ring
   R, x = PowerSeriesRing(ZZ, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -10:10)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end

   #  Inexact field
   R, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -1:1)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end

   # Non-integral domain
   T = ResidueRing(ZZ, 6)
   R, x = PowerSeriesRing(T, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, 0:5)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end
end

@testset "Generic.AbsSeries.binary_ops..." begin
   #  Exact ring
   R, x = PowerSeriesRing(ZZ, 10, "x", model=:capped_absolute)
   for iter = 1:100
      f = rand(R, 0:12, -10:10)
      g = rand(R, 0:12, -10:10)
      h = rand(R, 0:12, -10:10)
      @test isequal(f + g, g + f)
      @test isequal(f + (g + h), (f + g) + h)
      @test isequal(f*g, g*f)
      @test isequal(f*(g*h), (f*g)*h)
      @test isequal(f - g, -(g - f))
      @test isequal((f - h) + h, f)
      @test isequal(f*(g + h), f*g + f*h)
      @test isequal(f*(g - h), f*g - f*h)
   end

   #  Inexact field
   R, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)
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
   T = ResidueRing(ZZ, 6)
   R, x = PowerSeriesRing(T, 10, "x", model=:capped_absolute)
   for iter = 1:100
      f = rand(R, 0:12, 0:5)
      g = rand(R, 0:12, 0:5)
      h = rand(R, 0:12, 0:5)
      @test isequal(f + (g + h), (f + g) + h)
      @test isequal(f*g, g*f)
      @test isequal(f*(g*h), (f*g)*h)
      @test isequal(f - g, -(g - f))
      @test isequal((f - h) + h, f)
      @test isequal(f*(g + h), f*g + f*h)
      @test isequal(f*(g - h), f*g - f*h)
   end
end

@testset "Generic.AbsSeries.adhoc_binary_ops..." begin
   # Exact ring
   R, x = PowerSeriesRing(ZZ, 10, "x", model=:capped_absolute)
   for iter = 1:500
      f = rand(R, 0:12, -10:10)
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
   R, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)
   for iter = 1:500
      f = rand(R, 0:12, -1:1)
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
   S, x = PowerSeriesRing(R, 10, "x", model=:capped_absolute)
   for iter = 1:500
      f = rand(S, 0:12, 0:5)
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
   S, y = PowerSeriesRing(R, 10, "y", model=:capped_absolute)
   for iter = 1:100
      f = rand(S, 0:12, 0:5, -10:10)
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

@testset "Generic.AbsSeries.comparison..." begin
   # Exact ring
   R, x = PowerSeriesRing(ZZ, 10, "x", model=:capped_absolute)
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
   R, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)
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
   R = ResidueRing(ZZ, 6)
   S, x = PowerSeriesRing(R, 10, "x", model=:capped_absolute)
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
end

@testset "Generic.AbsSeries.adhoc_comparison..." begin
   # Exact ring
   R, x = PowerSeriesRing(ZZ, 10, "x", model=:capped_absolute)
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
   R, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)
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
   R = ResidueRing(ZZ, 6)
   S, x = PowerSeriesRing(R, 10, "x", model=:capped_absolute)
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
   S, y = PowerSeriesRing(R, 10, "y", model=:capped_absolute)
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

@testset "Generic.AbsSeries.powering..." begin
   # Exact ring
   R, x = PowerSeriesRing(ZZ, 10, "x", model=:capped_absolute)

   for iter = 1:100
      f = rand(R, 0:12, -10:10)
      r2 = R(1)

      for expn = 0:10
         r1 = f^expn

         @test (f == 0 && expn == 0 && r1 == 0) || isequal(r1, r2)

         r2 *= f
      end
   end

   # Inexact field
   R, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)

   for iter = 1:100
      f = rand(R, 0:12, -1:1)
      r2 = R(1)

      for expn = 0:4 # cannot set high power here
         r1 = f^expn

         @test (f == 0 && expn == 0 && r1 == 0) || isapprox(r1, r2)

         r2 *= f
      end
   end

   # Non-integral domain
   for iter = 1:100
      n = rand(2:26)

      Zn = ResidueRing(ZZ, n)
      R, x = PowerSeriesRing(Zn, 10, "x", model=:capped_absolute)

      f = rand(R, 0:12, 0:n - 1)
      r2 = R(1)

      for expn = 0:10
         r1 = f^expn

         @test (f == 0 && expn == 0 && r1 == 0) || isequal(r1, r2)

         r2 *= f
      end
   end
end

@testset "Generic.AbsSeries.shift..." begin
   # Exact ring
   R, x = PowerSeriesRing(ZZ, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -10:10)
      s = rand(0:12)
      g = rand(R, 0:0, -10:10) + O(x^s)

      set_prec!(g, max(precision(g), precision(f) + s))

      @test shift_right(shift_left(f, s) + g, s) == f
      @test isequal(shift_left(f, s), x^s*f)
      @test precision(shift_right(f, s)) == max(0, precision(f) - s)
   end
   
   f = rand(R, 0:12, -10:10)
   @test_throws DomainError shift_left(f, -1)
   @test_throws DomainError shift_left(f, -rand(2:100))
   @test_throws DomainError shift_right(f, -1)
   @test_throws DomainError shift_right(f, -rand(2:100))

   # Inexact field
   R, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -1:1)
      s = rand(0:12)
      g = rand(R, 0:0, -1:1) + O(x^s)

      set_prec!(g, max(precision(g), precision(f) + s))

      @test isapprox(shift_right(shift_left(f, s) + g, s), f)
      @test isapprox(shift_left(f, s), x^s*f)
      @test precision(shift_right(f, s)) == max(0, precision(f) - s)
   end
   
   f = rand(R, 0:12, -1:1)
   @test_throws DomainError shift_left(f, -1)
   @test_throws DomainError shift_left(f, -rand(2:100))
   @test_throws DomainError shift_right(f, -1)
   @test_throws DomainError shift_right(f, -rand(2:100))

   # Non-integral domain
   T = ResidueRing(ZZ, 6)
   R, x = PowerSeriesRing(T, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, 0:5)
      s = rand(0:12)
      g = rand(R, 0:0, 0:5) + O(x^s)

      set_prec!(g, max(precision(g), precision(f) + s))

      @test shift_right(shift_left(f, s) + g, s) == f
      @test isequal(shift_left(f, s), x^s*f)
      @test precision(shift_right(f, s)) == max(0, precision(f) - s)
   end
   
   f = rand(R, 0:12, 0:5)
   @test_throws DomainError shift_left(f, -1)
   @test_throws DomainError shift_left(f, -rand(2:100))
   @test_throws DomainError shift_right(f, -1)
   @test_throws DomainError shift_right(f, -rand(2:100))
end

@testset "Generic.AbsSeries.truncation..." begin
   # Exact ring
   R, x = PowerSeriesRing(ZZ, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -10:10)
      s = rand(0:12)

      @test truncate(f, s) == f
      @test isequal(truncate(f, s), f + O(x^s))
      @test precision(truncate(f, s)) == min(precision(f), s)
   end
   
   f = rand(R, 0:12, -10:10)
   @test_throws DomainError truncate(f, -1)
   @test_throws DomainError truncate(f, -rand(2:100))
   @test_throws DomainError truncate(f, -1)
   @test_throws DomainError truncate(f, -rand(2:100))

   # Inexact field
   R, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -1:1)
      s = rand(0:12)

      @test truncate(f, s) == f
      @test isapprox(truncate(f, s), f + O(x^s))
      @test precision(truncate(f, s)) == min(precision(f), s)
   end
   
   f = rand(R, 0:12, -1:1)
   @test_throws DomainError truncate(f, -1)
   @test_throws DomainError truncate(f, -rand(2:100))
   @test_throws DomainError truncate(f, -1)
   @test_throws DomainError truncate(f, -rand(2:100))

   # Non-integral domain
   T = ResidueRing(ZZ, 6)
   R, x = PowerSeriesRing(T, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, 0:5)
      s = rand(0:12)

      @test truncate(f, s) == f
      @test isequal(truncate(f, s), f + O(x^s))
      @test precision(truncate(f, s)) == min(precision(f), s)
   end
   
   f = rand(R, 0:12, 0:5)
   @test_throws DomainError truncate(f, -1)
   @test_throws DomainError truncate(f, -rand(2:100))
   @test_throws DomainError truncate(f, -1)
   @test_throws DomainError truncate(f, -rand(2:100))   
end

@testset "Generic.AbsSeries.inversion..." begin

    # Exact ring
    R, x = PowerSeriesRing(ZZ, 10, "x", model=:capped_absolute)
    for iter = 1:300
       f = R()
       while !isunit(f)
          f = rand(R, 0:0, -10:10)
       end

       @test f*inv(f) == 1
    end

    # Inexact field
    R, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)
    for iter = 1:300
       f = R()
       while coeff(f, 0) == 0
          f = rand(R, 0:0, -1:1)
       end

       @test isapprox(f*inv(f), R(1))
    end

    # Non-integral domain
    T = ResidueRing(ZZ, 6)
    R, x = PowerSeriesRing(T, 10, "x", model=:capped_absolute)
    for iter = 1:300
       f = R()
       while !isunit(f)
          f = rand(R, 0:0, 0:5)
       end

       @test f*inv(f) == 1
    end
 end

@testset "Generic.AbsSeries.square_root..." begin
    # Exact ring
    R, x = PowerSeriesRing(ZZ, 10, "x", model=:capped_absolute)
    for iter = 1:300
       f = rand(R, 0:10, -10:10)
       g = f^2

       @test isequal(sqrt(g)^2, g)
    end

    # Inexact field
    R, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)
    for iter = 1:300
       f = rand(R, 0:10, -1:1)
       g = f^2

       @test isapprox(sqrt(g)^2, g)
    end
end

@testset "Generic.AbsSeries.exact_division..." begin
   # Exact ring
   R, x = PowerSeriesRing(ZZ, 10, "x", model=:capped_absolute)
   for iter = 1:300
      s = rand(0:9)
      f = rand(R, s:s, -10:10)
      while valuation(f) != s || !isunit(coeff(f, s))
         f = rand(R, s:s, -10:10)
      end
      g = rand(R, s:s, -10:10)
      while valuation(g) != s || !isunit(coeff(g, s))
         g = rand(R, s:s, -10:10)
      end

      @test divexact(f, g)*g == f
   end

   # Inexact field
   R, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)
   for iter = 1:300
      s = rand(0:9)
      f = rand(R, s:s, -1:1)
      while valuation(f) != s
         f = rand(R, s:s, -1:1)
      end
      g = rand(R, s:s, -1:1)
      while valuation(g) != s || coeff(g, s) == 0
         g = rand(R, s:s, -1:1)
      end

      @test isapprox(divexact(f, g)*g, f)
   end

   # Non-integral domain
   T = ResidueRing(ZZ, 6)
   R, x = PowerSeriesRing(T, 10, "x", model=:capped_absolute)
   for iter = 1:300
      s = rand(0:9)
      f = rand(R, s:s, 0:5)
      while valuation(f) != s || !isunit(coeff(f, s))
         f = rand(R, s:s, 0:5)
      end
      g = rand(R, s:s, 0:5)
      while valuation(g) != s || !isunit(coeff(g, s))
         g = rand(R, s:s, 0:5)
      end

      @test divexact(f, g)*g == f
   end
end

@testset "Generic.AbsSeries.adhoc_exact_division..." begin
   # Exact field
   R, x = PowerSeriesRing(ZZ, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -10:10)
      c = ZZ()
      while c == 0
         c = rand(ZZ, -10:10)
      end

      @test isequal(divexact(f*c, c), f)
   end

   # Inexact field
   R, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -1:1)
      c = RealField()
      while isapprox(c, 0)
         c = rand(RealField, -1:1)
      end

      @test isapprox(divexact(f*c, c), f)
   end

   # Non-integral domain
   T = ResidueRing(ZZ, 6)
   R, x = PowerSeriesRing(T, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:0, 0:5)
      c = T()
      while !isunit(c)
         c = rand(T, 0:5)
      end

      @test isequal(divexact(f*c, c), f)
   end
end

@testset "Generic.AbsSeries.special_functions..." begin
   # Exact field
   S, x = PowerSeriesRing(QQ, 10, "x", model=:capped_absolute)

   for iter = 1:100
      @test exp(x + O(x^iter)) == exp(x + O(x^(iter - 1)))

      f = S()
      while !isunit(f)
         f = rand(S, 0:0, -10:10)
      end
      g = S()
      while !isunit(g) || !isunit(f + g)
         g = rand(S, 0:0, -10:10)
      end
      f *= x
      g *= x

      @test isequal(exp(f)*exp(g), exp(f + g))
   end

   # Inexact field
   S, x = PowerSeriesRing(RealField, 10, "x", model=:capped_absolute)

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

      @test isapprox(exp(f)*exp(g), exp(f + g))
   end

   # Non-integral domain
   R = ResidueRing(ZZ, 143)
   S, x = PowerSeriesRing(R, 10, "x", model=:capped_absolute)

   for iter = 1:10
      @test exp(x + O(x^iter)) == exp(x + O(x^(iter - 1)))

      f = S()
      while !isunit(coeff(f, 0))
         f = rand(S, 0:0, -10:10)
      end
      g = S()
      while !isunit(coeff(g, 0)) || !isunit(coeff(f + g, 0))
         g = rand(S, 0:0, -10:10)
      end
      f *= x
      g *= x

      @test isequal(exp(f)*exp(g), exp(f + g))
   end
end
