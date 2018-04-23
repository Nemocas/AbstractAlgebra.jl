type MyMap <: Map{AbstractAlgebra.Integers{BigInt}, AbstractAlgebra.Integers{BigInt}, SetMap, MyMap}
   a::Int
end

Generic.domain(f::Map(MyMap)) = AbstractAlgebra.JuliaZZ
Generic.codomain(f::Map(MyMap)) = AbstractAlgebra.JuliaZZ

a(f::Map(MyMap)) = Generic.get_field(f, :a)

(f::MyMap)(x) =  a(f)*(x + 1)

function test_gen_map_functional()
   print("Generic.Map.functional...")

   f = map_from_func(ZZ, ZZ, x -> x + 1)
   g = map_from_func(ZZ, QQ, x -> QQ(x))

   @test isa(f, Map(FunctionalMap))
   @test isa(g, Map(FunctionalMap))

   h = compose(f, g)

   @test isa(h, Map(FunctionalMap))

   k = f*g

   for i in 1:10
      @test h(ZZ(i)) == k(ZZ(i))
   end

   s = MyMap(2)

   t = compose(f, s)

   @test isa(t, Map(Generic.CompositeMap))     

   for i in 1:10
      @test t(ZZ(i)) == 2*(i + 2)
   end

   println("PASS")
end

function test_gen_map()
   test_gen_map_functional()

   println("")
end

