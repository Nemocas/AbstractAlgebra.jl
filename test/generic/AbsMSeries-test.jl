@testset "Generic.AbsMSeries.constructors" begin
   S, x = polynomial_ring(ZZ, "x")

   for nvars in 1:5
      prec = [rand(0:10) for i in 1:nvars]
      R, gens = power_series_ring(QQ, prec, ["x$(i)" for i in 1:nvars])
   
      f = rand(R, 0:10, -10:10)

      @test parent_type(f) == Generic.AbsMSeriesRing{Rational{BigInt},
                                               Generic.MPoly{Rational{BigInt}}}

      @test elem_type(R) == Generic.AbsMSeries{Rational{BigInt},
                                               Generic.MPoly{Rational{BigInt}}}

      @test R isa Generic.AbsMSeriesRing

      @test f isa Generic.AbsMSeries

      @test O(gens[rand(1:nvars)]^rand(1:10)) isa Generic.AbsMSeries

      @test O(R, 3) isa Generic.AbsMSeries
      @test O(R, 0) isa Generic.AbsMSeries

      @test R() isa Generic.AbsMSeries

      @test R(rand(-10:10)) isa Generic.AbsMSeries

      @test R(rand(QQ, -10:10)) isa Generic.AbsMSeries

      p = rand(R.poly_ring, 0:10, 0:10, -10:10)
      p = Generic.truncate_poly(p, prec)

      @test R(p, prec) isa Generic.AbsMSeries

      R, gens = power_series_ring(S, prec, ["x$(i)" for i in 1:nvars])
     
      @test R(ZZ(2)) isa Generic.AbsMSeries

      @test R(x) isa Generic.AbsMSeries
   end

   R1, (x1, y1) = power_series_ring(ZZ, [3, 3], ["x", "y"])
   R2, (x2, y2) = power_series_ring(ZZ, 3, ["x", "y"])

   @test R1 === R2

   # weighted

   for nvars in 1:5
      weights = [rand(1:4) for i in 1:nvars]
      prec = rand(0:40)
      R, gens = power_series_ring(QQ, weights, prec, ["x$(i)" for i in 1:nvars])
   
      f = rand(R, 0:10, -10:10)

      @test parent_type(f) == Generic.AbsMSeriesRing{Rational{BigInt},
                                               Generic.MPoly{Rational{BigInt}}}

      @test elem_type(R) == Generic.AbsMSeries{Rational{BigInt},
                                               Generic.MPoly{Rational{BigInt}}}

      @test R isa Generic.AbsMSeriesRing

      @test f isa Generic.AbsMSeries

      @test R() isa Generic.AbsMSeries

      @test R(rand(-10:10)) isa Generic.AbsMSeries

      @test R(rand(QQ, -10:10)) isa Generic.AbsMSeries

      p = rand(R.poly_ring, 0:10, 0:10, -10:10)
      p = Generic.truncate_poly(p, weights, prec)

      @test R(p, [0 for i in 1:nvars]) isa Generic.AbsMSeries

      R, gens = power_series_ring(S, weights, prec, ["x$(i)" for i in 1:nvars])
     
      @test R(ZZ(2)) isa Generic.AbsMSeries

      @test R(x) isa Generic.AbsMSeries
   end
end

@testset "Generic.AbsMSeries.printing" begin
   R, (x, y) = power_series_ring(ZZ, [5, 3], ["x", "y"])

   @test string(zero(R)) == "O(y^3) + O(x^5)"
   @test string(one(R)) == "1 + O(y^3) + O(x^5)"
   @test string(x) == "x + O(y^3) + O(x^5)"
   @test string(y) == "y + O(y^3) + O(x^5)"
   @test string(x^2 - y) == "-y + x^2 + O(y^3) + O(x^5)"

   # weighted 

   R, (x, y) = power_series_ring(ZZ, [2, 3], 10, ["x", "y"])

   @test string(zero(R)) == "O(10)"
   @test string(one(R)) == "1 + O(10)"
   @test string(x) == "x + O(10)"
   @test string(y) == "y + O(10)"
   @test string(x^2 - y) == "-y + x^2 + O(10)"

   @test !occursin("\n", sprint(show, R))
end

@testset "Generic.AbsMSeries.manipulation" begin
   R, (x, y) = power_series_ring(ZZ, [5, 3], ["x", "y"])

   f = x^2 + 2x*y + y + 1 + O(x^3) + O(y^2)

   @test Generic.poly(f) isa Generic.MPoly{BigInt}

   @test length(f) == 4

   @test nvars(R) == 2

   @test precision(f) == [3, 2]

   f = set_precision!(f, [2, 2])

   @test precision(f) == [2, 2]

   @test max_precision(R) == [5, 3]

   @test valuation(x^2 + 2*x^2*y + y) == [2, 1]

   @test valuation(x^2) == [2, 3]

   @test valuation(R()) == [5, 3]

   @test coeff(f, 3) == 2
   @test coeff(f, 1) == 1
   @test_throws BoundsError coeff(f, 0)
   @test_throws BoundsError coeff(f, 5)

   @test iszero(zero(R))
   @test isone(one(R))
   @test !iszero(one(R))
   @test !isone(zero(R))
   @test zero(R) == 0
   @test one(R) == 1

   @test is_unit(1 + y + x)
   @test !is_unit(x)
   @test !is_unit(2 + x)

   @test is_gen(gen(R, 1))
   @test is_gen(gen(R, 2))
   @test_throws ArgumentError gen(R, 0)
   @test_throws ArgumentError gen(R, 5)
   @test !is_gen(x^2)
   @test !is_gen(R(1))

   @test gen(R,1) == R[1]

   @test gens(R) == [x, y]

   @test parent(x) == R

   @test symbols(R) == [:x, :y]

   @test base_ring(x) == ZZ
   @test base_ring(R) == ZZ

   @test deepcopy(f) == f

   @test hash(f) isa UInt

   @test characteristic(R) == 0

   S, = residue_ring(ZZ, 7)
   R, (x, y) = power_series_ring(S, [2, 3], ["x", "y"])

   @test characteristic(R) == 7

   # weighted

   R, (x, y) = power_series_ring(ZZ, [2, 3], 10, ["x", "y"])

   f = x^2 + 2x*y + y + 1

   @test Generic.poly(f) isa Generic.MPoly{BigInt}

   @test length(f) == 4

   @test nvars(R) == 2

   @test precision(f) == 10

   @test coeff(f, 3) == 2
   @test coeff(f, 1) == 1
   @test_throws BoundsError coeff(f, 0)
   @test_throws BoundsError coeff(f, 5)

   @test iszero(zero(R))
   @test isone(one(R))
   @test !iszero(one(R))
   @test !isone(zero(R))
   @test zero(R) == 0
   @test one(R) == 1

   @test is_unit(1 + y + x)
   @test !is_unit(x)
   @test !is_unit(2 + x)

   @test is_gen(gen(R, 1))
   @test is_gen(gen(R, 2))
   @test_throws ArgumentError gen(R, 0)
   @test_throws ArgumentError gen(R, 5)
   @test !is_gen(x^2)
   @test !is_gen(R(1))

   @test gen(R,1) == R[1]

   @test gens(R) == [x, y]

   @test parent(x) == R

   @test symbols(R) == [:x, :y]

   @test base_ring(x) == ZZ
   @test base_ring(R) == ZZ

   @test deepcopy(f) == f

   @test hash(f) isa UInt

   @test characteristic(R) == 0

   S, = residue_ring(ZZ, 7)
   R, (x, y) = power_series_ring(S, [2, 3], 10, ["x", "y"])

   @test characteristic(R) == 7
end

@testset "Generic.AbsMSeries.truncation" begin
   for nvars in 1:5
      prec = [rand(0:10) for i in 1:nvars]
      R, gens = power_series_ring(ZZ, prec, ["x$(i)" for i in 1:nvars])

      for iters = 1:100
         f = rand(R, 0:10, -10:10)

         prec2 = [rand(0:10) for i in 1:nvars]
         g = truncate(f, prec2)
         @test g == f
         @test precision(g) == min.(prec, prec2)

         prec3 = rand(0:10)
         g = truncate(f, prec3)
         @test g == f
         @test precision(g) == [min(x, prec3) for x in prec]
      end
   end
end

@testset "Generic.AbsMSeries.unary_ops" begin
   for nvars in 1:5
      prec = [rand(0:10) for i in 1:nvars]
      R, gens = power_series_ring(ZZ, prec, ["x$(i)" for i in 1:nvars])

      for iters = 1:100
         f = rand(R, 0:10, -10:10)

         @test isequal(f, -(-f))
         @test iszero(f + (-f))
      end
   end

   # weighted

   for nvars in 1:5
      weights = [rand(1:4) for i in 1:nvars]
      prec = rand(0:20)
      R, gens = power_series_ring(ZZ, weights, prec, ["x$(i)" for i in 1:nvars])

      for iters = 1:100
         f = rand(R, 0:10, -10:10)

         @test isequal(f, -(-f))
         @test iszero(f + (-f))
      end
   end
end

@testset "Generic.AbsMSeries.binary_ops" begin
   for nvars in 1:5
      prec = [rand(0:10) for i in 1:nvars]
      R, gens = power_series_ring(ZZ, prec, ["x$(i)" for i in 1:nvars])

      for iters = 1:100
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
   end

   # weighted

   for nvars in 1:5
      weights = [rand(1:4) for i in 1:nvars]
      prec = rand(0:20)
      R, gens = power_series_ring(ZZ, weights, prec, ["x$(i)" for i in 1:nvars])

      for iters = 1:100
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
   end
end

@testset "Generic.AbsMSeries.adhoc_binary_ops" begin
   S, x = polynomial_ring(ZZ, "x")

   for nvars in 1:5
      prec = [rand(0:10) for i in 1:nvars]
      R, gens = power_series_ring(ZZ, prec, ["x$(i)" for i in 1:nvars])

      for iters = 1:100
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
   
      prec = [rand(0:10) for i in 1:nvars]
      R, gens = power_series_ring(S, prec, ["x$(i)" for i in 1:nvars])

      for iters = 1:100
         f = rand(R, 0:12, -1:2, -10:10)
         c1 = rand(S, -1:2, -10:10)
         c2 = rand(S, -1:2, -10:10)

         @test isequal(c1*f - c2*f, (c1 - c2)*f)
         @test isequal(c1*f + c2*f, (c1 + c2)*f)

         @test isequal(f*c1 - f*c2, f*(c1 - c2))
         @test isequal(f*c1 + f*c2, f*(c1 + c2))

         @test isequal(f + c1, c1 + f)
         @test isequal(c1 * f, f * c1)
      end
   end

   # weighted

   for nvars in 1:5
      weights = [rand(1:4) for i in 1:nvars]
      prec = rand(0:20)
      R, gens = power_series_ring(ZZ, weights, prec, ["x$(i)" for i in 1:nvars])

      for iters = 1:100
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
   
      weights = [rand(1:4) for i in 1:nvars]
      prec = rand(0:20)
      R, gens = power_series_ring(S, weights, prec, ["x$(i)" for i in 1:nvars])

      for iters = 1:100
         f = rand(R, 0:12, -1:2, -10:10)
         c1 = rand(S, -1:2, -10:10)
         c2 = rand(S, -1:2, -10:10)

         @test isequal(c1*f - c2*f, (c1 - c2)*f)
         @test isequal(c1*f + c2*f, (c1 + c2)*f)

         @test isequal(f*c1 - f*c2, f*(c1 - c2))
         @test isequal(f*c1 + f*c2, f*(c1 + c2))

         @test isequal(f + c1, c1 + f)
         @test isequal(c1 * f, f * c1)
      end
   end
end

@testset "Generic.AbsMSeries.powering" begin
   for nvars in 1:5
      prec = [rand(0:10) for i in 1:nvars]
      R, gens = power_series_ring(ZZ, prec, ["x$(i)" for i in 1:nvars])
      
      for iter = 1:100
         f = rand(R, 0:12, -10:10)
         r2 = R(1)

         for expn = 0:10
            r1 = f^expn

            @test isequal(r1, r2)

            r2 *= f
         end
      end
   end

   # weighted

   for nvars in 1:5
      weights = [rand(1:4) for i in 1:nvars]
      prec = rand(0:20)
      R, gens = power_series_ring(ZZ, weights, prec, ["x$(i)" for i in 1:nvars])
      
      for iter = 1:100
         f = rand(R, 0:12, -10:10)
         r2 = R(1)

         for expn = 0:10
            r1 = f^expn

            @test isequal(r1, r2)

            r2 *= f
         end
      end
   end
end

@testset "Generic.AbsMSeries.comparison" begin
   for nvars in 1:5
      prec = [rand(0:10) for i in 1:nvars]
      R, gens = power_series_ring(ZZ, prec, ["x$(i)" for i in 1:nvars])

      for iter = 1:100
         f = rand(R, 0:12, -10:10)

         @test f == deepcopy(f)
         @test isequal(f, deepcopy(f))

         prec2 = [rand(0:10) for i in 1:nvars]

         g = truncate(f, prec2)

         @test f == g
      end
   end

   # weighted

   for nvars in 1:5
      weights = [rand(1:4) for i in 1:nvars]
      prec = rand(0:20)
      R, gens = power_series_ring(ZZ, weights, prec, ["x$(i)" for i in 1:nvars])

      for iter = 1:100
         f = rand(R, 0:12, -10:10)

         @test f == deepcopy(f)
         @test isequal(f, deepcopy(f))
      end
   end
end

@testset "Generic.AbsMSeries.adhoc_comparison" begin
   for nvars in 1:5
      prec = [rand(1:10) for i in 1:nvars]
      R, gens = power_series_ring(ZZ, prec, ["x$(i)" for i in 1:nvars])

      for iter = 1:100
         c = rand(ZZ, -10:10)

         @test c == R(c)
         @test R(c) == c

         c = rand(zz, -10:10)

         @test c == R(c)
         @test R(c) == c
      end
   end

   # weighted

   for nvars in 1:5
      weights = [rand(1:4) for i in 1:nvars]
      prec = rand(0:20)
      R, gens = power_series_ring(ZZ, weights, prec, ["x$(i)" for i in 1:nvars])

      for iter = 1:100
         c = rand(ZZ, -10:10)

         @test c == R(c)
         @test R(c) == c

         c = rand(zz, -10:10)

         @test c == R(c)
         @test R(c) == c
      end
   end
end

@testset "Generic.AbsMSeries.inversion" begin
   for nvars in 1:4
      prec = [rand(1:8) for i in 1:nvars]
      R, gens = power_series_ring(QQ, prec, ["x$(i)" for i in 1:nvars])

      for iter = 1:20*nvars
         f = rand(R, 0:8, -10:10)
         if is_unit(f)
            @test isequal(f*inv(f), one(R))
         end
      end
   end

   # weighted

   for nvars in 1:4
      weights = [rand(1:4) for i in 1:nvars]
      prec = rand(0:20)
      R, gens = power_series_ring(QQ, weights, prec, ["x$(i)" for i in 1:nvars])
      for iter = 1:20*nvars
         f = rand(R, 0:8, -10:10)
         if is_unit(f)
            @test isequal(f*inv(f), one(R))
         end
      end
   end
end

@testset "Generic.AbsMSeries.exact_division" begin
   for nvars in 1:4
      prec = [rand(1:8) for i in 1:nvars]
      R, gens = power_series_ring(QQ, prec, ["x$(i)" for i in 1:nvars])

      for iter = 1:20*nvars
         f = rand(R, 0:8, -10:10)
         g = rand(R, 0:8, -10:10)
         if is_unit(g)
            @test isequal(divexact(f, g)*g, f)
         end
      end
   end

   # weighted

   for nvars in 1:4
      weights = [rand(1:4) for i in 1:nvars]
      prec = rand(0:20)
      R, gens = power_series_ring(QQ, weights, prec, ["x$(i)" for i in 1:nvars])

      for iter = 1:20*nvars
         f = rand(R, 0:8, -10:10)
         g = rand(R, 0:8, -10:10)
         if is_unit(g)
            @test isequal(divexact(f, g)*g, f)
         end
      end
   end
end

@testset "Generic.AbsMSeries.evaluation" begin
   R, (x, y) = power_series_ring(ZZ, [10, 10], ["x", "y"])

   f = x^2 + y

   @test evaluate(f, [x^3 + 1, y^2]) == 1 + y^2 + 2*x^3 + x^6 + O(y^10) + O(x^10)

   for iters = 1:30
      f = rand(R, 10:10, -10:10)
      g = rand(R, 10:10, -10:10)
      h = rand(R, 10:10, -10:10)
      k = rand(R, 10:10, -10:10)

      @test evaluate(f + k, [g, h]) == evaluate(f, [g, h]) + evaluate(k, [g, h])

      @test evaluate(f, [1], [g]) == evaluate(f, [g, y])
      @test evaluate(f, [2], [h]) == evaluate(f, [x, h])
      @test evaluate(f, [x], [g]) == evaluate(f, [g, y])
      @test evaluate(f, [y], [h]) == evaluate(f, [x, h])
   end

   # weighted

   R, (x, y) = power_series_ring(ZZ, [2, 3], 30, ["x", "y"])

   f = x^2 + y

   @test evaluate(f, [x^3 + 1, y^2]) == 1 + y^2 + 2*x^3 + x^6

   for iters = 1:30
      f = rand(R, 10:10, -10:10)
      g = rand(R, 10:10, -10:10)
      h = rand(R, 10:10, -10:10)
      k = rand(R, 10:10, -10:10)

      @test evaluate(f + k, [g, h]) == evaluate(f, [g, h]) + evaluate(k, [g, h])

      @test evaluate(f, [1], [g]) == evaluate(f, [g, y])
      @test evaluate(f, [2], [h]) == evaluate(f, [x, h])
      @test evaluate(f, [x], [g]) == evaluate(f, [g, y])
      @test evaluate(f, [y], [h]) == evaluate(f, [x, h])
   end
end
