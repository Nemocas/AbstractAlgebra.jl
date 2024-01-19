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

@testset "Generic.AbsSeries.types" begin
   @test abs_series_type(BigInt) == Generic.AbsSeries{BigInt}
   @test abs_series_type(Rational{BigInt}) == Generic.AbsSeries{Rational{BigInt}}
end

@testset "Generic.AbsSeries.constructors" begin
   S1 = AbsPowerSeriesRing(ZZ, 10)
   S2 = AbsPowerSeriesRing(ZZ, 10)

   @test S1 !== S2
   @test isa(S1, Generic.AbsPowerSeriesRing)

   R, x = power_series_ring(ZZ, 30, "x", model=:capped_absolute)

   S, t = polynomial_ring(QQ, "t")
   T, y = power_series_ring(S, 30, "y", model=:capped_absolute)

   @test elem_type(R) == Generic.AbsSeries{BigInt}
   @test elem_type(Generic.AbsPowerSeriesRing{BigInt}) == Generic.AbsSeries{BigInt}
   @test parent_type(Generic.AbsSeries{BigInt}) == Generic.AbsPowerSeriesRing{BigInt}

   @test elem_type(T) == Generic.AbsSeries{elem_type(S)}
   @test elem_type(Generic.AbsPowerSeriesRing{elem_type(S)}) == Generic.AbsSeries{elem_type(S)}
   @test parent_type(Generic.AbsSeries{elem_type(S)}) == Generic.AbsPowerSeriesRing{elem_type(S)}

   @test isa(R, Generic.AbsPowerSeriesRing)

   @test isa(T, Generic.AbsPowerSeriesRing)

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
   c3 = R([1, 3, 5], 3, 5)

   @test isa(c1, Generic.AbsSeries)
   @test isa(c2, Generic.AbsSeries)
   @test isa(c3, Generic.AbsSeries)

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

   R, x = power_series_ring(ZZ, 30, "x", model=:capped_absolute)
   S, x = power_series_ring(ZZ, 30, "x", model=:capped_absolute)
   @test R === S
   S, x = power_series_ring(ZZ, 30, "x", model=:capped_absolute, cached = false)
   @test R !== S
end

@testset "Generic.AbsSeries.manipulation" begin
   R, t = polynomial_ring(QQ, "t")
   S, x = power_series_ring(R, 30, "x", model=:capped_absolute)

   @test max_precision(S) == 30

   a = 2x + x^3
   b = O(x^4)

   @test valuation(a) == 1
   @test valuation(b) == 4

   @test precision(a) == 30
   @test precision(b) == 4

   @test is_gen(gen(S))

   @test iszero(zero(S))

   @test isone(one(S))

   @test is_unit(-1 + x + 2x^2)

   @test isequal(deepcopy(a), a)
   @test isequal(deepcopy(b), b)

   @test normalise(a, 4) == 4

   @test coeff(a, 1) == 2
   @test polcoeff(a, 1) == 2
   @test coeff(b, 7) == 0
   @test polcoeff(b, 7) == 0

   @test_throws DomainError coeff(a, -1)
   @test_throws DomainError polcoeff(a, -1)
   @test_throws DomainError coeff(a, -rand(2:100))
   @test_throws DomainError polcoeff(a, -rand(2:100))
   
   @test characteristic(R) == 0

   T, = residue_ring(ZZ, 7)
   U, y = power_series_ring(T, 10, "y", model=:capped_absolute)

   @test modulus(T) == 7
end

@testset "Generic.AbsSeries.similar" begin
   R, x = power_series_ring(ZZ, 10, "x"; model=:capped_absolute)

   for iters = 1:10
      f = rand(R, 0:10, -10:10)

      g = similar(f, QQ, "y")
      h = similar(f, "y")
      k = similar(f)
      m = similar(f, QQ, 5)
      n = similar(f, 5)

      @test isa(g, AbsPowerSeriesRingElem)
      @test isa(h, AbsPowerSeriesRingElem)
      @test isa(k, AbsPowerSeriesRingElem)
      @test isa(m, AbsPowerSeriesRingElem)
      @test isa(n, AbsPowerSeriesRingElem)

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
   end
end

@testset "Generic.AbsSeries.abs_series" begin
   f = abs_series(ZZ, [1, 2, 3], 3, 5, "y")

   @test isa(f, AbsPowerSeriesRingElem)
   @test base_ring(f) === ZZ
   @test coeff(f, 0) == 1
   @test coeff(f, 2) == 3
   @test parent(f).S == :y

   g = abs_series(ZZ, [1, 2, 3], 3, 5)

   @test isa(g, AbsPowerSeriesRingElem)
   @test base_ring(g) === ZZ
   @test coeff(g, 0) == 1
   @test coeff(g, 2) == 3
   @test parent(g).S == :x

   h = abs_series(ZZ, [1, 2, 3], 2, 5)
   k = abs_series(ZZ, [1, 2, 3], 1, 6, cached=false)
   m = abs_series(ZZ, [1, 2, 3], 3, 9, cached=false)

   @test parent(h) === parent(g)
   @test parent(k) !== parent(m)

   p = abs_series(ZZ, BigInt[], 0, 4)
   q = abs_series(ZZ, [], 0, 6)

   @test isa(p, AbsPowerSeriesRingElem)
   @test isa(q, AbsPowerSeriesRingElem)

   @test length(p) == 0
   @test pol_length(p) == 0
   @test length(q) == 0
   @test pol_length(q) == 0

   r = abs_series(QQ, BigInt[1, 2, 3], 3, 5)

   @test isa(r, AbsPowerSeriesRingElem)

   s = abs_series(ZZ, [1, 2, 3], 3, 5; max_precision=10)
   
   @test max_precision(parent(s)) == 10
end

@testset "Generic.AbsSeries.unary_ops" begin
   #  Exact ring
   R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -10:10)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end

   #  Inexact field
   R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -1:1)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = power_series_ring(T, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, 0:5)

      @test isequal(-(-f), f)
      @test iszero(f + (-f))
   end
end

@testset "Generic.AbsSeries.binary_ops" begin
   #  Exact ring
   R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
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
   R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)
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
   T, = residue_ring(ZZ, 6)
   R, x = power_series_ring(T, 10, "x", model=:capped_absolute)
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

@testset "Generic.AbsSeries.adhoc_binary_ops" begin
   # Exact ring
   R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
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
   R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)
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
   R, = residue_ring(ZZ, 6)
   S, x = power_series_ring(R, 10, "x", model=:capped_absolute)
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
   S, y = power_series_ring(R, 10, "y", model=:capped_absolute)
   for iter = 1:100
      f = rand(S, 0:12, -1:5, -10:10)
      c1 = rand(ZZ, -10:10)
      c2 = rand(ZZ, -10:10)
      d1 = rand(R, -1:5, -10:10)
      d2 = rand(R, -1:5, -10:10)

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

@testset "Generic.AbsSeries.comparison" begin
   # Exact ring
   R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
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
   R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)
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
   R, = residue_ring(ZZ, 6)
   S, x = power_series_ring(R, 10, "x", model=:capped_absolute)
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

@testset "Generic.AbsSeries.adhoc_comparison" begin
   # Exact ring
   R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
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
   R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)
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
   S, x = power_series_ring(R, 10, "x", model=:capped_absolute)
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
   S, y = power_series_ring(R, 10, "y", model=:capped_absolute)
   for iter = 1:100
      f = S()
      while f == 0
         f = rand(S, 0:0, 0:5, -10:10)
      end
      f += rand(S, 1:12, -1:5, -10:10)
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

@testset "Generic.AbsSeries.powering" begin
   # Exact ring
   R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)

   for iter = 1:100
      f = rand(R, 0:12, -10:10)
      r2 = R(1)

      for expn = 0:10
         r1 = f^expn

         @test isequal(r1, r2)

         r2 *= f
      end
   end

   # Inexact field
   R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)

   for iter = 1:100
      f = rand(R, 0:12, -1:1)
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
      R, x = power_series_ring(Zn, 10, "x", model=:capped_absolute)

      f = rand(R, 0:12, 0:n - 1)
      r2 = R(1)

      for expn = 0:10
         r1 = f^expn

         @test isequal(r1, r2)

         r2 *= f
      end
   end
end

@testset "Generic.AbsSeries.shift" begin
   # Exact ring
   R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -10:10)
      s = rand(0:12)
      g = rand(R, 0:0, -10:10) + O(x^s)

      g = set_precision!(g, max(precision(g), precision(f) + s))

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
   R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -1:1)
      s = rand(0:12)
      g = rand(R, 0:0, -1:1) + O(x^s)

      g = set_precision!(g, max(precision(g), precision(f) + s))

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
   T, = residue_ring(ZZ, 6)
   R, x = power_series_ring(T, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, 0:5)
      s = rand(0:12)
      g = rand(R, 0:0, 0:5) + O(x^s)

      g = set_precision!(g, max(precision(g), precision(f) + s))

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

@testset "Generic.AbsSeries.truncation" begin
   # Exact ring
   R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
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
   R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)
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
   T, = residue_ring(ZZ, 6)
   R, x = power_series_ring(T, 10, "x", model=:capped_absolute)
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

@testset "Generic.AbsSeries.inversion" begin

    # Exact ring
    R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
    for iter = 1:300
       f = R()
       while !is_unit(f)
          f = rand(R, 0:0, -10:10)
       end

       @test f*inv(f) == 1
    end

    # Exact field
    for prec = 1:10
       R, x = power_series_ring(QQ, prec, "x"; model=:capped_absolute)
       for iter = 1:30
          f = R()
          while valuation(f) != 0
             f = rand(R, 0:0, -10:10)
          end

          @test f*inv(f) == 1
       end
    end

    # Inexact field
    R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)
    for iter = 1:300
       f = R()
       while coeff(f, 0) == 0
          f = rand(R, 0:0, -1:1)
       end

       @test isapprox(f*inv(f), R(1))
    end

    # Non-integral domain
    T, = residue_ring(ZZ, 6)
    R, x = power_series_ring(T, 10, "x", model=:capped_absolute)
    for iter = 1:300
       f = R()
       while !is_unit(f)
          f = rand(R, 0:0, 0:5)
       end

       @test f*inv(f) == 1
    end
end

@testset "Generic.AbsSeries.compose" begin
    # Exact ring
    R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
    for iter = 1:300
        f1 = rand(R, 0:10, -10:10)
        f2 = rand(R, 0:10, -10:10)
     
        g = rand(R, 1:10, -10:10)

        @test compose(f1 + f2, g) == compose(f1, g) + compose(f2, g)
        @test compose(x, g) == g
        @test compose(x^2, g) == g^2
        @test compose(R(), g) == R()
    end

    S, y = power_series_ring(ZZ, 10, "y", model=:capped_absolute)
    for iter = 1:300
        f1 = rand(R, 0:10, -10:10)
        f2 = rand(R, 0:10, -10:10)

        g = rand(S, 1:10, -10:10)

        @test compose(f1 + f2, g) == compose(f1, g) + compose(f2, g)
        @test compose(R(), g) == S()
    end

    # Inexact field
    R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)
    for iter = 1:300
        f1 = rand(R, 0:10, -10:10)
        f2 = rand(R, 0:10, -10:10)

        g = rand(R, 1:10, -10:10)

        @test isapprox(compose(f1 + f2, g), compose(f1, g) + compose(f2, g))
        @test isapprox(compose(x, g), g)
        @test isapprox(compose(x^2, g), g^2)
        @test isapprox(compose(R(), g), R())
    end

    S, y = power_series_ring(RealField, 10, "y", model=:capped_absolute)
    for iter = 1:300
        f1 = rand(R, 0:10, -10:10)
        f2 = rand(R, 0:10, -10:10)

        g = rand(S, 1:10, -10:10)
        @test isapprox(compose(f1 + f2, g), compose(f1, g) + compose(f2, g))
        @test isapprox(compose(R(), g), S())
    end

    # Non-integral domain
    T, = residue_ring(ZZ, 6)
    R, x = power_series_ring(T, 10, "x", model=:capped_absolute)
    for iter = 1:300
        f1 = rand(R, 0:10, -10:10)
        f2 = rand(R, 0:10, -10:10)

        g = rand(R, 1:10, -10:10)

        @test compose(f1 + f2, g) == compose(f1, g) + compose(f2, g)
        @test compose(x, g) == g
        @test compose(x^2, g) == g^2
        @test compose(R(), g) == R()
    end

    S, y = power_series_ring(T, 10, "y", model=:capped_absolute)
    for iter = 1:300
        f1 = rand(R, 0:10, -10:10)
        f2 = rand(R, 0:10, -10:10)

        g = rand(S, 1:10, -10:10)

        @test compose(f1 + f2, g) == compose(f1, g) + compose(f2, g)
        @test compose(R(), g) == S()
    end
end

@testset "Generic.AbsSeries.square_root" begin
    # Exact ring
    R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
    for iter = 1:300
       f = rand(R, 0:10, -10:10)
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
    R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)
    for iter = 1:300
       f = rand(R, 0:10, -1:1)
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
    for p in [2, 7, 19, 65537, ZZ(2), ZZ(7), ZZ(19), ZZ(65537)]
        R, = residue_field(ZZ, p)

        S, x = power_series_ring(R, 10, "x", model=:capped_absolute)

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

    S, x = power_series_ring(T, 10, "x", model=:capped_absolute)

    f = 1 + y^2*x^2 + (y^2 + y + 1)*x^4 + O(x^10)

    @test_throws ErrorException sqrt(f)

    f = x + y^4*x^4 + (y^4 + y^2 + 1)*x^8 + O(x^11)

    @test_throws ErrorException sqrt(f)

    f = zero(S)

    @test iszero(sqrt(f))
end

@testset "Generic.AbsSeries.exact_division" begin
   # Exact ring
   R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
   for iter = 1:300
      s = rand(0:9)
      f = rand(R, s:s, -10:10)
      while valuation(f) != s || !is_unit(coeff(f, s))
         f = rand(R, s:s, -10:10)
      end
      g = rand(R, s:s, -10:10)
      while valuation(g) != s || !is_unit(coeff(g, s))
         g = rand(R, s:s, -10:10)
      end

      @test divexact(f, g)*g == f
   end

   # Exact field
   R, x = power_series_ring(QQ, 10, "x", model=:capped_absolute)
   for iter = 1:300
      s = rand(0:9)
      f = rand(R, s:s, -10:10)
      while valuation(f) != s || !is_unit(coeff(f, s))
         f = rand(R, s:s, -10:10)
      end
      g = rand(R, s:s, -10:10)
      while valuation(g) != s || !is_unit(coeff(g, s))
         g = rand(R, s:s, -10:10)
      end

      @test divexact(f, g)*g == f
   end

   # Inexact field
   R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)
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
   T, = residue_ring(ZZ, 6)
   R, x = power_series_ring(T, 10, "x", model=:capped_absolute)
   for iter = 1:300
      s = rand(0:9)
      f = rand(R, s:s, 0:5)
      while valuation(f) != s || !is_unit(coeff(f, s))
         f = rand(R, s:s, 0:5)
      end
      g = rand(R, s:s, 0:5)
      while valuation(g) != s || !is_unit(coeff(g, s))
         g = rand(R, s:s, 0:5)
      end

      @test divexact(f, g)*g == f
   end
end

@testset "Generic.AbsSeries.adhoc_exact_division" begin
   # Exact field
   R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -10:10)
      c = ZZ()
      while c == 0
         c = rand(ZZ, -10:10)
      end

      @test isequal(divexact(f*c, c), f)
   end

   # Inexact field
   R, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:12, -1:1)
      c = RealField()
      while isapprox(c, 0)
         c = rand(RealField, -1:1)
      end

      @test isapprox(divexact(f*c, c), f)
   end

   # Non-integral domain
   T, = residue_ring(ZZ, 6)
   R, x = power_series_ring(T, 10, "x", model=:capped_absolute)
   for iter = 1:300
      f = rand(R, 0:0, 0:5)
      c = T()
      while !is_unit(c)
         c = rand(T, 0:5)
      end

      @test isequal(divexact(f*c, c), f)
   end
end

@testset "Generic.AbsSeries.derivative_integral" begin
   # Exact field
   S, x = power_series_ring(QQ, 10, "x"; model=:capped_absolute)

   for iter = 1:100
      f = rand(S, 0:10, -10:10)

      const_coeff = S(coeff(f, 0))
      set_precision!(const_coeff, precision(f))

      @test isequal(integral(derivative(f)) + const_coeff, f)
   end

   # Inexact field
   S, x = power_series_ring(RealField, 10, "x"; model=:capped_absolute)

   for iter = 1:100
      f = rand(S, 0:10, -10:10)

      const_coeff = S(coeff(f, 0))
      set_precision!(const_coeff, precision(f))

      @test isapprox(integral(derivative(f)) + const_coeff, f)
   end
end

@testset "Generic.AbsSeries.special_functions" begin
   # Exact field
   S, x = power_series_ring(QQ, 10, "x", model=:capped_absolute)

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
   S, x = power_series_ring(RealField, 10, "x", model=:capped_absolute)

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

      @test isapprox(log(exp(f)), f)

      while coeff(f, 0) <= 0
         f = rand(S, 0:0, -10:10)
      end

      @test isapprox(exp(log(f)), f)

   end

   # Non-integral domain
   R, = residue_ring(ZZ, 143)
   S, x = power_series_ring(R, 10, "x", model=:capped_absolute)

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

@testset "Generic.AbsSeries.change_base_ring" begin
   Zx, x = power_series_ring(ZZ, 10, "x"; model=:capped_absolute)
   @test 1 == map_coefficients(sqrt, x^0)
   p = Zx([i for i in 1:10], 10, 11)
   q = Zx([i for i in 10:-1:1], 10, 11)
   pq = p * q
   for R in [QQ, GF(2), GF(13), ZZ]
      pR = change_base_ring(R, p)
      qR = change_base_ring(R, q, parent = parent(pR))
      @test parent(qR) === parent(pR)
      pqR = change_base_ring(R, pq, parent = parent(pR))
      @test pR * qR == pqR
   end

   ps = map_coefficients(z -> z^2, p)
   @test ps == Zx([i^2 for i in 1:10], 10, 11)

   f = x^2 + 3x^3 + 2x^6
   @test map_coefficients(one, f) == x^2 + x^3 + x^6
   f2 = map_coefficients(t -> t + 2, f)
   @test f2 == 3x^2 + 5x^3 + 4x^6
   for i in [0, 1, 4, 5]
      @test coeff(f2, i) !== coeff(f, i)
   end

   F = GF(11)
   P, y = power_series_ring(F, 10, "x"; model=:capped_absolute)
   @test map_coefficients(t -> F(t) + 2, f) == 3y^2 + 5y^3 + 4y^6
end

@testset "Generic.AbsSeries.unsafe_operators" begin
   # Exact ring
   R, x = power_series_ring(ZZ, 10, "x", model=:capped_absolute)
   
   for iter = 1:300
      f = rand(R, 0:9, -10:10)
      g = rand(R, 0:9, -10:10)
      f0 = deepcopy(f)
      g0 = deepcopy(g)

      h = rand(R, 0:9, -10:10)

      k = f + g
      h = add!(h, f, g)
      @test isequal(h, k)
      @test isequal(f, f0)
      @test isequal(g, g0)

      f1 = deepcopy(f)
      f1 = add!(f1, f1, g)
      @test isequal(f1, k)
      @test isequal(g, g0)

      g1 = deepcopy(g)
      g1 = add!(g1, f, g1)
      @test isequal(g1, k)
      @test isequal(f, f0)

      f1 = deepcopy(f)
      f1 = addeq!(f1, g)
      @test isequal(h, k)
      @test isequal(g, g0)

      k = f*g
      h = mul!(h, f, g)
      @test isequal(h, k)
      @test isequal(f, f0)
      @test isequal(g, g0)      
     
      f1 = deepcopy(f)
      f1 = mul!(f1, f1, g)
      @test isequal(f1, k)
      @test isequal(g, g0)

      g1 = deepcopy(g)
      g1 = mul!(g1, f, g1)
      @test isequal(g1, k)
      @test isequal(f, f0)

      h = zero!(h)
      @test isequal(h, R())
   end
end

@testset "Generic.AbsSeries.euclidean" begin
   R, x = power_series_ring(QQ, 20, "x", model = :capped_absolute)
   S, y = power_series_ring(GF(5), 20, "y", model = :capped_absolute)

   T, z = power_series_ring(GF(7), 20, "z", model = :capped_absolute)
   @test_throws ErrorException divrem(y, z)

   for (T, t) in [ (R, x), (S, y) ]
      @test divrem(t^2, t) == (t, zero(T))
      @test divrem(t + 1, t) == (zero(T), t + 1)
      @test_throws DivideError divrem(t, zero(T))

      @test gcd(t^2*(t - 1), t^2) == t^2
      @test gcd((t + 1)*(t - 1), t + 1) == one(T)
      @test gcdx(t^2*(t - 1), t^2) == (t^2, zero(T), one(T))
      g, u, v = gcdx((t + 1)*(t - 1), t + 1)
      @test u*(t + 1)*(t - 1) + v*(t + 1) == g

      M = matrix(T, 4, 3, [ 0 0 0; t^3 + 1 t^2 0; 0 t^2 t^5; t^4 + 1 t^2 t^5 + t^3 ])
      @test is_hnf(hnf(M))
      H, U = hnf_with_transform(M)
      @test is_hnf(H)
      @test U*M == H
      @test is_unit(det(U))

      @test is_snf(snf(M))
      S, U, V = snf_with_transform(M)
      @test is_snf(S)
      @test U*M*V == S
      @test is_unit(det(U))
      @test is_unit(det(V))
   end
 end
