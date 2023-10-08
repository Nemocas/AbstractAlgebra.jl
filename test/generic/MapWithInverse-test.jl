import AbstractAlgebra.PrettyPrinting

@testset "Generic.MapWithInverse.constructors" begin
   f = map_from_func(x -> x + 1, ZZ, ZZ)
   g = map_from_func(x -> x - 1, ZZ, ZZ)

   s = map_with_section(f, g)

   @test domain(s) == AbstractAlgebra.JuliaZZ
   @test codomain(s) == AbstractAlgebra.JuliaZZ

   @test image_fn(s)(ZZ(1)) == 2
   @test inverse_fn(s)(ZZ(1)) == 0

   @test image_map(s)(ZZ(1)) == 2
   @test section_map(s)(ZZ(1)) == 0

   @test isa(s, Map(Generic.MapWithSection))

   t = map_with_retraction(f, g)

   @test domain(t) == AbstractAlgebra.JuliaZZ
   @test codomain(t) == AbstractAlgebra.JuliaZZ

   @test image_fn(t)(ZZ(1)) == 2
   @test inverse_fn(t)(ZZ(1)) == 0

   @test image_map(t)(ZZ(1)) == 2
   @test retraction_map(t)(ZZ(1)) == 0

   @test isa(t, Map(Generic.MapWithRetraction))

   u = map_with_section_from_func(x -> x + 1, x -> x - 1, ZZ, ZZ)

   @test u(ZZ(1)) == 2

   v = map_with_section_from_func(x -> x + 1, x -> x - 1, ZZ, ZZ)

   @test u(ZZ(1)) == 2
end

@testset "Generic.MapWithInverse.composition" begin
   f = map_from_func(x -> x + 1, ZZ, ZZ)
   g = map_from_func(x -> x - 1, ZZ, ZZ)
   h = map_from_func(x -> x + 2, ZZ, ZZ)
   k = map_from_func(x -> x - 2, ZZ, ZZ)

   s = map_with_section(f, g)
   t = map_with_section(h, k)

   u = compose(s, t)

   @test u(ZZ(1)) == 4

   s = map_with_retraction(f, g)
   t = map_with_retraction(h, k)

   u = compose(s, t)

   @test u(ZZ(1)) == 4
end

@testset "Generic.MapWithInverse.inv" begin
   f = map_from_func(x -> x + 1, ZZ, ZZ)
   g = map_from_func(x -> x - 1, ZZ, ZZ)
   h = map_from_func(x -> x + 2, ZZ, ZZ)
   k = map_from_func(x -> x - 2, ZZ, ZZ)

   s = map_with_section(f, g)
   t = map_with_section(h, k)

   u = compose(s, t)

   v = inv(u)
   w = inv(s)

   @test v(ZZ(1)) == -2
   @test w(ZZ(1)) == 0

   s = map_with_retraction(f, g)
   t = map_with_retraction(h, k)

   u = compose(s, t)

   v = inv(u)
   w = inv(s)

   @test v(ZZ(1)) == -2
   @test w(ZZ(1)) == 0
end

@testset "Generic.MapWithInverse.printing" begin
  F = GF(5)
  u = map_with_section_from_func(x -> F(x + 1), x -> lift(x - 1), ZZ, F)
  str = """
        Map with section
          from integers
          to finite field F_5"""
  @test PrettyPrinting.detailed(u) == str
  @test PrettyPrinting.oneline(u) == "Map: integers -> finite field F_5"
  @test PrettyPrinting.supercompact(u) == "Map with section"

  v = map_with_retraction_from_func(x -> QQ(x + 1), x -> ZZ(x - 1), ZZ, QQ)
  str = """
        Map with retraction
          from integers
          to rationals"""
  @test PrettyPrinting.detailed(v) == str
  @test PrettyPrinting.oneline(v) == "Map: integers -> rationals"
  @test PrettyPrinting.supercompact(v) == "Map with retraction"
end
