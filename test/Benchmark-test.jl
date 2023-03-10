@testset "Benchmark.fateman" begin
   R, x = polynomial_ring(ZZ, "x")
   S, y = polynomial_ring(R, "y")
   T, z = polynomial_ring(S, "z")
   U, t = polynomial_ring(T, "t")

   p = (x + y + z + t + 1)^10

   q = p*(p + 1)

   @test length(q) == 21
end

@testset "Benchmark.pearce" begin
   R, x = polynomial_ring(ZZ, "x")
   S, y = polynomial_ring(R, "y")
   T, z = polynomial_ring(S, "z")
   U, t = polynomial_ring(T, "t")
   V, u = polynomial_ring(U, "u")

   f = (x + y + 2z^2 + 3t^3 + 5u^5 + 1)^7
   g = (u + t + 2z^2 + 3y^3 + 5x^5 + 1)^7

   q = f*g

   @test length(q) == 43
end
