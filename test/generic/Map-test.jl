import AbstractAlgebra.PrettyPrinting

module MyMapMod
   using AbstractAlgebra
   mutable struct MyMap <: Map{AbstractAlgebra.Integers{BigInt}, AbstractAlgebra.Integers{BigInt}, SetMap, MyMap}
      a::Int
   end

   Generic.domain(f::MyMap) = AbstractAlgebra.JuliaZZ
   Generic.codomain(f::MyMap) = AbstractAlgebra.JuliaZZ

   a(f::MyMap) = f.a

   (f::MyMap)(x) =  a(f)*(x + 1)
end

@testset "Generic.Map.FunctionalMap" begin
   f = map_from_func(x -> x + 1, ZZ, ZZ)
   g = map_from_func(x -> QQ(x), ZZ, QQ)

   @test isa(f, Map(FunctionalMap))
   @test isa(g, Map(FunctionalMap))

   @test domain(f) == AbstractAlgebra.JuliaZZ
   @test codomain(f) == AbstractAlgebra.JuliaZZ
   @test domain(g) == AbstractAlgebra.JuliaZZ
   @test codomain(g) == AbstractAlgebra.JuliaQQ

   @test image(f, ZZ(1)) == 2
   @test image(g, ZZ(2)) == QQ(2)

   @test image_fn(f)(ZZ(1)) == 2
   @test image_fn(g)(ZZ(2)) == QQ(2)
end

@testset "Generic.Map.FunctionalCompositeMap" begin
   f = map_from_func(x -> x + 1, ZZ, ZZ)
   g = map_from_func(x -> QQ(x), ZZ, QQ)

   h = compose(f, g)

   @test isa(h, Map(FunctionalMap))
   @test isa(h, Map(Generic.FunctionalCompositeMap))

   k = f*g

   for i in 1:10
      @test h(ZZ(i)) == k(ZZ(i))
      @test image(h, ZZ(i)) == image(k, ZZ(i))
   end

   @test domain(h) == AbstractAlgebra.JuliaZZ
   @test codomain(h) == AbstractAlgebra.JuliaQQ

   @test image(h, ZZ(1)) == QQ(2)
   @test image_fn(h)(ZZ(1)) == QQ(2)

   @test map1(h) === f
   @test map2(h) === g
end

@testset "Generic.Map.CompositeMap" begin
   f = map_from_func(x -> x + 1, ZZ, ZZ)

   s = MyMapMod.MyMap(2)

   t = compose(f, s)

   @test isa(t, Map(Generic.CompositeMap))

   for i in 1:10
      @test t(ZZ(i)) == 2*(i + 2)
      @test image(t, ZZ(i)) == 2*(i + 2)
   end

   @test domain(t) == AbstractAlgebra.JuliaZZ
   @test codomain(t) == AbstractAlgebra.JuliaZZ

   @test map1(t) === f
   @test map2(t) === s

   V = free_module(QQ, 2)
   f = AbstractAlgebra.ModuleHomomorphism(V, V, QQ[1 1; 1 2])
   g = identity_map(V)
   h = f * g
   @test domain(h) == V
   @test codomain(h) == V
end

@testset "Generic.Map.IdentityMap" begin
   f = map_from_func(x -> QQ(x + 1), ZZ, QQ)
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

@testset "Generic.Map.printing" begin
  id = identity_map(ZZ)
  str = """
        Identity map
          of integers"""
  @test PrettyPrinting.repr_detailed(id) == str
  @test PrettyPrinting.repr_oneline(id) == "Identity map of integers"
  @test PrettyPrinting.repr_terse(id) == "Identity map"

  u = map_from_func(x -> QQ(x + 1), ZZ, QQ)
  str = """
        Map defined by a Julia function
          from integers
          to rationals"""
  @test PrettyPrinting.repr_detailed(u) == str
  @test PrettyPrinting.repr_oneline(u) == "Map: integers -> rationals"
  @test PrettyPrinting.repr_terse(u) == "Map defined by a Julia function"

  f = map_from_func(x -> x + 1, ZZ, ZZ)
  g = map_from_func(x -> QQ(x), ZZ, QQ)
  v = compose(f, g)
  str = """
        Functional composite map
          from integers
          to rationals
        which is the composite of
          Map: integers -> integers
          Map: integers -> rationals"""
  @test PrettyPrinting.repr_detailed(v) == str
  @test PrettyPrinting.repr_oneline(v) == "Map: integers -> integers -> rationals"
  @test PrettyPrinting.repr_terse(v) == "Functional composite map"

  f = map_from_func(x -> x + 1, ZZ, ZZ)
  s = MyMapMod.MyMap(2)
  t = compose(f, s)
  str = """
        Composite map
          from integers
          to integers
        which is the composite of
          Map: integers -> integers
          Map: integers -> integers"""
  @test PrettyPrinting.repr_detailed(t) == str
  @test PrettyPrinting.repr_oneline(t) == "Map: integers -> integers -> integers"
  @test PrettyPrinting.repr_terse(t) == "Composite map"

end

@testset "Generic.Map.broadcasting" begin
   id = identity_map(ZZ)
   @test id.([ZZ(1), ZZ(2), ZZ(3)]) == [ZZ(1), ZZ(2), ZZ(3)]
   @test [ZZ(1), ZZ(2), ZZ(3)] .|> id == [ZZ(1), ZZ(2), ZZ(3)]

   f = map_from_func(x -> x + 1, ZZ, ZZ)
   @test f.([ZZ(1), ZZ(2), ZZ(3)]) == [ZZ(2), ZZ(3), ZZ(4)]
   @test [ZZ(1), ZZ(2), ZZ(3)] .|> f == [ZZ(2), ZZ(3), ZZ(4)]

   s = MyMapMod.MyMap(2)
   @test s.([ZZ(1), ZZ(2), ZZ(3)]) == [ZZ(4), ZZ(6), ZZ(8)]
   @test [ZZ(1), ZZ(2), ZZ(3)] .|> s == [ZZ(4), ZZ(6), ZZ(8)]

   t = compose(f, s)
   @test t.([ZZ(1), ZZ(2), ZZ(3)]) == [ZZ(6), ZZ(8), ZZ(10)]
   @test [ZZ(1), ZZ(2), ZZ(3)] .|> t == [ZZ(6), ZZ(8), ZZ(10)]
end
