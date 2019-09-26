@testset "Generic.MapWithInverse.constructors..." begin
   f = map_from_func(ZZ, ZZ, x -> x + 1)
   g = map_from_func(ZZ, ZZ, x -> x - 1)

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

   u = map_with_section_from_func(ZZ, ZZ, x -> x + 1, x -> x - 1)

   @test u(ZZ(1)) == 2

   v = map_with_retraction_from_func(ZZ, ZZ, x -> x + 1, x -> x - 1)

   @test u(ZZ(1)) == 2
end

@testset "Generic.MapWithInverse.composition..." begin
   f = map_from_func(ZZ, ZZ, x -> x + 1)
   g = map_from_func(ZZ, ZZ, x -> x - 1)
   h = map_from_func(ZZ, ZZ, x -> x + 2)
   k = map_from_func(ZZ, ZZ, x -> x - 2)

   s = map_with_section(f, g)
   t = map_with_section(h, k)

   u = compose(s, t)

   @test u(ZZ(1)) == 4

   s = map_with_retraction(f, g)
   t = map_with_retraction(h, k)

   u = compose(s, t)

   @test u(ZZ(1)) == 4
end

@testset "Generic.MapWithInverse.inv..." begin
   f = map_from_func(ZZ, ZZ, x -> x + 1)
   g = map_from_func(ZZ, ZZ, x -> x - 1)
   h = map_from_func(ZZ, ZZ, x -> x + 2)
   k = map_from_func(ZZ, ZZ, x -> x - 2)

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
