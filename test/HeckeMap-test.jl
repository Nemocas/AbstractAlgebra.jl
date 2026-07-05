import AbstractAlgebra.PrettyPrinting

@testset "MapFromFunc" begin
   f = map_from_func(ZZ, ZZ, x -> x + 1)
   g = map_from_func(ZZ, QQ, x -> QQ(x))

   @test isa(f, Map(MapFromFunc))
   @test isa(g, Map(MapFromFunc))

   @test domain(f) == AbstractAlgebra.JuliaZZ
   @test codomain(f) == AbstractAlgebra.JuliaZZ
   @test domain(g) == AbstractAlgebra.JuliaZZ
   @test codomain(g) == AbstractAlgebra.JuliaQQ

   @test image(f, ZZ(1)) == 2
   @test image(g, ZZ(2)) == QQ(2)

   @test image_fn(f)(ZZ(1)) == 2
   @test image_fn(g)(ZZ(2)) == QQ(2)

   h = compose(f, g)

   @test isa(h, Map(Generic.CompositeMap))

   k = f*g

   for i in 1:10
      @test h(ZZ(i)) == k(ZZ(i))
      @test image(h, ZZ(i)) == image(k, ZZ(i))
   end

   @test domain(h) == AbstractAlgebra.JuliaZZ
   @test codomain(h) == AbstractAlgebra.JuliaQQ

   @test image(h, ZZ(1)) == QQ(2)
   # no `image_fn` defined for Generic.CompositeMap - is this a problem?
   # @test image_fn(h)(ZZ(1)) == QQ(2)

   @test map1(h) === f
   @test map2(h) === g
end


@testset "MapFromFunc IdentityMap" begin
   f = map_from_func(ZZ, QQ, x -> QQ(x + 1))
   g = identity_map(ZZ)
   h = identity_map(QQ)

   @test isa(g, Map(IdentityMap))
   @test isa(h, Map(IdentityMap))

   @test compose(g, f) === f
   @test compose(f, h) === f

   @test domain(g) == AbstractAlgebra.JuliaZZ
   @test codomain(g) == AbstractAlgebra.JuliaZZ
   @test domain(h) == AbstractAlgebra.JuliaQQ
   @test codomain(h) == AbstractAlgebra.JuliaQQ

   for i = 1:10
      @test g(ZZ(i)) == ZZ(i)
      @test h(ZZ(i)//(i + 1)) == ZZ(i)//(i + 1)
      @test image(g, ZZ(i)) == ZZ(i)
      @test image(h, ZZ(i)//(i + 1)) == ZZ(i)//(i + 1)
   end

   @test compose(g, g) === g
   @test compose(h, h) === h

   @test inv(g) === g
   @test inv(h) === h
end


@testset "MapFromFunc printing" begin
  u = map_from_func(ZZ, QQ, x -> QQ(x + 1))
  str = """
        Map defined by a Julia function
          from integers
          to rationals"""
  @test PrettyPrinting.repr_detailed(u) == str
  @test PrettyPrinting.repr_oneline(u) == "Map: integers -> rationals"
  @test PrettyPrinting.repr_terse(u) == "Map"

  f = map_from_func(ZZ, ZZ, x -> x + 1)
  g = map_from_func(ZZ, QQ, x -> QQ(x))
  v = compose(f, g)
  str = """
        Composite map
          from integers
          to rationals
        which is the composite of
          Map: integers -> integers
          Map: integers -> rationals"""
  @test PrettyPrinting.repr_detailed(v) == str
  @test PrettyPrinting.repr_oneline(v) == "Map: integers -> integers -> rationals"
  @test PrettyPrinting.repr_terse(v) == "Composite map"
end


@testset "MapFromFunc broadcasting" begin
   f = map_from_func(ZZ, ZZ, x -> x + 1)
   @test f.([ZZ(1), ZZ(2), ZZ(3)]) == [ZZ(2), ZZ(3), ZZ(4)]
   @test [ZZ(1), ZZ(2), ZZ(3)] .|> f == [ZZ(2), ZZ(3), ZZ(4)]
end


@testset "MapFromFunc with inverse" begin
  f = map_from_func(ZZ, ZZ, x -> x+1, x -> x-1)

  @test image(f, ZZ(1)) == 2
  @test preimage(f, ZZ(1)) == 0

  @test image_fn(f)(ZZ(1)) == 2
  @test inverse_fn(f)(ZZ(1)) == 0

  @test preimage(f, image(f, ZZ(1))) == ZZ(1)

  finv = inv(f)
  y = finv(ZZ(1))
  @test preimage(finv, y) == ZZ(1)
end


@testset "MapFromFunc attributes" begin
  f = map_from_func(ZZ, ZZ, x -> x+1, x -> x-1)

  @test is_attribute_storing(f)

  set_attribute!(f, :has_inverse, true)
  @test get_attribute(f, :has_inverse)
  @test get_attribute(f.header, :has_inverse)
end


@testset "MapHeader constructions" begin
   M1 = AbstractAlgebra.MapHeader{typeof(ZZ), typeof(QQ)}()
   M2 = AbstractAlgebra.MapHeader{typeof(ZZ), typeof(QQ)}(ZZ, QQ)
   M3 = AbstractAlgebra.MapHeader(ZZ, QQ)
   @test typeof(M1) == typeof(M2)
   @test typeof(M2) == typeof(M3)
   @test M2.domain == M3.domain
   @test M2.codomain == M3.codomain
end
