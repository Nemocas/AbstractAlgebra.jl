@testset "Fac.constructors" begin
   f = Fac(-1, Dict{Int, Int}(2 => 3, 3 => 1))
   ff = Fac(-1, [2 => 3, 3 => 1])

   @test -24 == evaluate(f)
   @test -24 == evaluate(ff)
   @test -24 == unit(f)*prod([p^e for (p, e) in f])
   @test -24 == unit(ff)*prod([p^e for (p, e) in ff])

   @test collect(ff) == [2 => 3, 3 => 1]

   ff = Fac(-1, [2 => 3, 3 => 1])
   @test ff[ZZ(2)] == 3
end

@testset "Fac.printing" begin
   f = Fac(-1, Dict{Int, Int}(2 => 3, 3 => 1))
   ff = Fac(-1, [(2, 3), (3, 1)])

   @test string(f) == "-1 * 2^3 * 3" || string(f) == "-1 * 3 * 2^3"
   @test string(ff) == "-1 * 2^3 * 3"

   @test string(Fac{BigInt}()) isa String

   R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
   f = Fac(x, Dict(x*y => 1, x + y => 1))

   @test evaluate(f) == x * (x + y) * (x*y)
   @test string(f) == "x * (x + y) * (x*y)" ||
         string(f) == "x * (x*y) * (x + y)"
end

@testset "Fac.equality" begin
  # 2025-10-20  equality test on factorizations always gives error (except if the args are ===)
   f = Fac(-1, Dict{Int, Int}(2 => 3, 3 => 1))
   ff = Fac(-1, [2 => 3, 3 => 1])
   fzz = Fac(ZZ(-1), Dict{ZZRingElem, Int}(ZZ(2) => 3, ZZ(3) => 1))

   @test  f == f
   @test_throws  MethodError  (f == ff) 
   @test_throws  MethodError  (f == fzz)
end
