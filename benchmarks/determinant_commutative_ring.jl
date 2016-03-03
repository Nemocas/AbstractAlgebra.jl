function benchmark_znz_det()
   print("benchmark_znz_det ... ")
   n = 2003 * 1009
   Zn = ResidueRing(FlintZZ, n)
   R, x = PolynomialRing(Zn, "x")

   M = MatrixSpace(R, 80, 80)()

   for i in 1:80
     for j in 1:80
       for k in 0:5
         M[i, j] = M[i, j] + x^k * Zn(rand(-100:100))
       end
     end
   end

   tt = @elapsed determinant(M)
   println("$tt")
end

