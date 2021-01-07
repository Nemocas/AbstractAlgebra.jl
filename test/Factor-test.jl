@testset "Fac.constructors..." begin
   f = Fac(-1, Dict{Int, Int}(2 => 3, 3 => 1))

   @test -24 == unit(f)*prod([p^e for (p, e) in f])
end

@testset "Fac.printing..." begin
   f = Fac(-1, Dict{Int, Int}(2 => 3, 3 => 1))

   @test string(f) == "-2^3*3" || string(f) == "-3*2^3"

   @test string(Fac{BigInt}()) isa String
end
