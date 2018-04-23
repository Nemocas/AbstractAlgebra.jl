function test_gen_map_cache_constructors()
   print("Generic.MapCache.constructors...")

   f = cached(map_from_func(ZZ, ZZ, x -> x + 1))
   g = cached(map_from_func(ZZ, QQ, x -> QQ(x)), limit=2)
   h = cached(map_from_func(ZZ, ZZ, x -> x + 2), enabled=false)

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

   println("PASS")
end

function test_gen_map_cache_enable_disable()
   print("Generic.MapCache.enable_disable...")

   f = cached(map_from_func(ZZ, ZZ, x -> x + 1))

   @test f(ZZ(1)) == 2
   
   disable_cache!(f)
   
   @test f(ZZ(1)) == 2
   @test f(ZZ(2)) == 3

   enable_cache!(f)

   @test f(ZZ(1)) == 2
   @test f(ZZ(3)) == 4

   g = cached(map_from_func(ZZ, ZZ, x -> x + 1), enabled=false)

   @test g(ZZ(1)) == 2

   enable_cache!(g)

   @test g(ZZ(1)) == 2
   @test g(ZZ(2)) == 3

   disable_cache!(g)

   @test g(ZZ(1)) == 2
   @test g(ZZ(3)) == 4
   
   println("PASS")
end

function test_gen_map_cache_limit()
   print("Generic.MapCache.limit...")

   g = cached(map_from_func(ZZ, QQ, x -> QQ(x)), limit=2)

   set_limit!(g, 3)

   for i = 1:10
      @test g(ZZ(i)) == QQ(i)
   end

   println("PASS")
end

function test_gen_map_cache()
   test_gen_map_cache_constructors()
   test_gen_map_cache_enable_disable()
   test_gen_map_cache_limit()

   println("")
end

