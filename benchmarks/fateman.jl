function benchmark_fateman()
   print("benchmark_fateman ... ")
   R, x = PolynomialRing(FlintZZ, "x")
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")
   U, t = PolynomialRing(T, "t")

   p = (x + y + z + t + 1)^20
   
   tt = @elapsed p*(p + 1)
   println("$tt")
end

