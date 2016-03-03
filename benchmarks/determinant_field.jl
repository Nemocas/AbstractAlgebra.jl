function benchmark_nf_det()
   print("benchmark_nf_det ... ")
   QQx, x = PolynomialRing(FlintQQ, "x")
   K, a = AnticNumberField(x^3 + 3*x + 1, "a")
   M = MatrixSpace(K, 80, 80)()

   for i in 1:80
     for j in 1:80
       for k in 0:2
         M[i, j] = M[i, j] + a^k * (rand(-100:100))
       end
     end
   end

   tt = @elapsed determinant(M)
   println("$tt")
end

