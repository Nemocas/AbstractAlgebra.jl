@testset "Generic.AbsMSeries.constructors" begin
   S, x = PolynomialRing(ZZ, "x")

   for nvars in 1:5
      prec = [rand(0:10) for i in 1:nvars]
      R, gens = PowerSeriesRing(QQ, prec, ["x$(i)" for i in 1:nvars])
   
      f = rand(R, 0:10, -10:10)

      @test parent_type(f) == Generic.AbsMSeriesRing{Rational{BigInt}}

      @test elem_type(R) == Generic.AbsMSeries{Rational{BigInt}}

      @test R isa Generic.AbsMSeriesRing

      @test f isa Generic.AbsMSeries

      @test O(gens[rand(1:nvars)]^rand(1:10)) isa Generic.AbsMSeries

      @test R() isa Generic.AbsMSeries

      @test R(rand(-10:10)) isa Generic.AbsMSeries

      @test R(rand(QQ, -10:10)) isa Generic.AbsMSeries

      p = rand(R.poly_ring, 0:10, 0:10, -10:10)
      p = Generic.truncate_poly(p, prec)

      @test R(p, prec) isa Generic.AbsMSeries

      R, gens = PowerSeriesRing(S, prec, ["x$(i)" for i in 1:nvars])
     
      @test R(ZZ(2)) isa Generic.AbsMSeries

      @test R(x) isa Generic.AbsMSeries
   end
end

@testset "Generic.AbsMSeries.printing" begin
   R, (x, y) = PowerSeriesRing(ZZ, [5, 3], ["x", "y"])

   @test string(zero(R)) == "O(y^3) + O(x^5)"
   @test string(one(R)) == "1 + O(y^3) + O(x^5)"
   @test string(x) == "x + O(y^3) + O(x^5)"
   @test string(y) == "y + O(y^3) + O(x^5)"
   @test string(x^2 - y) == "-y + x^2 + O(y^3) + O(x^5)"
end

@testset "Generic.AbsMSeries.manipulation" begin
   R, (x, y) = PowerSeriesRing(ZZ, [5, 3], ["x", "y"])

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

   @test isunit(1 + y + x)
   @test !isunit(x)
   @test !isunit(2 + x)

   @test isgen(gen(R, 1))
   @test isgen(gen(R, 2))
   @test !isgen(x^2)
   @test !isgen(R(1))

   @test gens(R) == [x, y]

   @test parent(x) == R

   @test vars(R) == [:x, :y]

   @test base_ring(x) == ZZ
   @test base_ring(R) == ZZ

   @test deepcopy(f) == f

   @test hash(f) isa UInt

   @test characteristic(R) == 0

   S = ResidueRing(ZZ, 7)
   R, (x, y) = PowerSeriesRing(S, [2, 3], ["x", "y"])

   @test characteristic(R) == 7
end

@testset "Generic.AbsMSeries.truncation" begin
   for nvars in 1:5
      prec = [rand(0:10) for i in 1:nvars]
      R, gens = PowerSeriesRing(ZZ, prec, ["x$(i)" for i in 1:nvars])

      f = rand(R, 0:10, -10:10)

      prec2 = [rand(0:10) for i in 1:nvars]

      g = truncate(f, prec2)

      @test precision(g) == min.(prec, prec2)
   end
end

@testset "Generic.AbsMSeries.unary_ops" begin
   for nvars in 1:5
      prec = [rand(0:10) for i in 1:nvars]
      R, gens = PowerSeriesRing(ZZ, prec, ["x$(i)" for i in 1:nvars])

      f = rand(R, 0:10, -10:10)

      @test isequal(f, -(-f))
      @test iszero(f + (-f))
   end
end
