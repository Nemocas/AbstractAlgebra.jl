@testset "Generic.MapCache.constructors..." begin
   f = cached(map_from_func(x -> x + 1, ZZ, ZZ))
   g = cached(map_from_func(x -> QQ(x), ZZ, QQ), limit=2)
   h = cached(map_from_func(x -> x + 2, ZZ, ZZ), enabled=false)

   @test isa(f, Map(Generic.MapCache))
   @test isa(g, Map(Generic.MapCache))
   @test isa(h, Map(Generic.MapCache))

   @test f(ZZ(1)) == 2

   for i = 1:10
      @test g(ZZ(i)) == QQ(i)
   end

   @test h(ZZ(1)) == 3

   @test domain(f) == AbstractAlgebra.JuliaZZ
   @test codomain(g) == AbstractAlgebra.JuliaQQ
   @test domain(h) == AbstractAlgebra.JuliaZZ

   @test image_fn(f)(ZZ(1)) == 2
end

@testset "Generic.MapCache.enable_disable..." begin
   f = cached(map_from_func(x -> x + 1, ZZ, ZZ))

   @test f(ZZ(1)) == 2

   disable_cache!(f)

   @test f(ZZ(1)) == 2
   @test f(ZZ(2)) == 3

   enable_cache!(f)

   @test f(ZZ(1)) == 2
   @test f(ZZ(3)) == 4

   g = cached(map_from_func(x -> x + 1, ZZ, ZZ), enabled=false)

   @test g(ZZ(1)) == 2

   enable_cache!(g)

   @test g(ZZ(1)) == 2
   @test g(ZZ(2)) == 3

   disable_cache!(g)

   @test g(ZZ(1)) == 2
   @test g(ZZ(3)) == 4
end

@testset "Generic.MapCache.limit..." begin
   g = cached(map_from_func(x -> QQ(x), ZZ, QQ), limit=2)

   set_limit!(g, 3)

   for i = 1:10
      @test g(ZZ(i)) == QQ(i)
   end
end
