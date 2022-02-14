mutable struct MyMap <: Map{AbstractAlgebra.Integers{BigInt}, AbstractAlgebra.Integers{BigInt}, SetMap, MyMap}
   a::Int
end

Generic.domain(f::MyMap) = AbstractAlgebra.JuliaZZ
Generic.codomain(f::MyMap) = AbstractAlgebra.JuliaZZ

a(f::MyMap) = f.a

(f::MyMap)(x) =  a(f)*(x + 1)

@testset "Generic.Map.FunctionalMap" begin
   f = map_from_func(x -> x + 1, ZZ, ZZ)
   g = map_from_func(x -> QQ(x), ZZ, QQ)

   @test isa(f, Map(FunctionalMap))
   @test isa(g, Map(FunctionalMap))

   @test domain(f) == AbstractAlgebra.JuliaZZ
   @test codomain(f) == AbstractAlgebra.JuliaZZ
   @test domain(g) == AbstractAlgebra.JuliaZZ
   @test codomain(g) == AbstractAlgebra.JuliaQQ

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
   end

   @test domain(h) == AbstractAlgebra.JuliaZZ
   @test codomain(h) == AbstractAlgebra.JuliaQQ

   @test image_fn(h)(ZZ(1)) == QQ(2)

   @test map1(h) === f
   @test map2(h) === g
end

@testset "Generic.Map.CompositeMap" begin
   f = map_from_func(x -> x + 1, ZZ, ZZ)

   s = MyMap(2)

   t = compose(f, s)

   @test isa(t, Map(Generic.CompositeMap))

   for i in 1:10
      @test t(ZZ(i)) == 2*(i + 2)
   end

   @test domain(t) == AbstractAlgebra.JuliaZZ
   @test codomain(t) == AbstractAlgebra.JuliaZZ

   @test map1(t) === f
   @test map2(t) === s

   V = FreeModule(QQ, 2)
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
   end

   @test compose(g, g) === g
   @test compose(h, h) === h

   @test inv(g) === g
   @test inv(h) === h
end
