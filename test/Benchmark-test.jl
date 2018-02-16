function test_benchmark_fateman()
   print("Benchmark.fateman...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")
   U, t = PolynomialRing(T, "t")

   p = (x + y + z + t + 1)^10
   
   q = p*(p + 1)

   @test length(q) == 21

   println("PASS")
end

function test_benchmark_pearce()
   print("Benchmark.pearce...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")
   U, t = PolynomialRing(T, "t")
   V, u = PolynomialRing(U, "u")

   f = (x + y + 2z^2 + 3t^3 + 5u^5 + 1)^7
   g = (u + t + 2z^2 + 3y^3 + 5x^5 + 1)^7
   
   q = f*g

   @test length(q) == 43

   println("PASS")
end

function test_benchmarks()
   test_benchmark_fateman()
   test_benchmark_pearce()

   println("")
end

