@testset "Fac.constructors..." begin
   f = Fac(-1, Dict{Int, Int}(2 => 3, 3 => 1))

   @test -24 == unit(f)*prod([p^e for (p, e) in f])
end

@testset "Fac.printing..." begin
   f = Fac(-1, Dict{Int, Int}(2 => 3, 3 => 1))

   @test string(f) == "-1 * 2^3 * 3" || string(f) == "-1 * 3 * 2^3"

   @test string(Fac{BigInt}()) isa String

   R, (x, y) = PolynomialRing(ZZ, ["x", "y"])
   f = Fac(x, Dict(x*y => 1, x + y => 1))

   @test string(f) == "x * (x + y) * (x*y)" ||
         string(f) == "x * (x*y) * (x + y)"
end
