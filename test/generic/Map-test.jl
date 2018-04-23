type MyMap <: Map{AbstractAlgebra.Integers{BigInt}, AbstractAlgebra.Integers{BigInt}, SetMap, MyMap}
   a::Int
end

Generic.domain(f::Map(MyMap)) = AbstractAlgebra.JuliaZZ
Generic.codomain(f::Map(MyMap)) = AbstractAlgebra.JuliaZZ

a(f::Map(MyMap)) = Generic.get_field(f, :a)

(f::MyMap)(x) =  a(f)*(x + 1)

function test_gen_map_functional()
   print("Generic.Map.FunctionalMap...")

   f = map_from_func(ZZ, ZZ, x -> x + 1)
   g = map_from_func(ZZ, QQ, x -> QQ(x))

   @test isa(f, Map(FunctionalMap))
   @test isa(g, Map(FunctionalMap))

   @test domain(f) == AbstractAlgebra.JuliaZZ
   @test codomain(f) == AbstractAlgebra.JuliaZZ
   @test domain(g) == AbstractAlgebra.JuliaZZ
   @test codomain(g) == AbstractAlgebra.JuliaQQ

   @test image_fn(f)(ZZ(1)) == 2
   @test image_fn(g)(ZZ(2)) == QQ(2)

   println("PASS")
end

function test_gen_map_functional_composite()
   print("Generic.Map.FunctionalCompositeMap...")

   f = map_from_func(ZZ, ZZ, x -> x + 1)
   g = map_from_func(ZZ, QQ, x -> QQ(x))

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

   println("PASS")
end

function test_gen_map_composite()
   print("Generic.Map.CompositeMap...")

   f = map_from_func(ZZ, ZZ, x -> x + 1)

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

   println("PASS")
end

function test_gen_map()
   test_gen_map_functional()
   test_gen_map_functional_composite()
   test_gen_map_composite()

   println("")
end

