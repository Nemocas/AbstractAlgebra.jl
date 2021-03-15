@testset "Generic.AbsMSeries.constructors" begin
   S, x = PolynomialRing(ZZ, "x")

   for nvars in 1:5
      prec = [rand(0:10) for i in 1:nvars]
      R, gens = PowerSeriesRing(QQ, prec, ["x$(i)" for i in 1:nvars])
      
      f = rand(R, 0:10, -10:10)

      @test R isa Generic.AbsMSeriesRing

      @test f isa Generic.AbsMSeries

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
