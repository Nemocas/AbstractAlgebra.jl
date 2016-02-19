function benchmark_fateman()
   R, x = PolynomialRing(FlintZZ, "x")
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")
   U, t = PolynomialRing(T, "t")

   p = (x + y + z + t + 1)^20
   
   q = p*(p + 1)
end

function benchmark_pearce()
   R, x = PolynomialRing(FlintZZ, "x")
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")
   U, t = PolynomialRing(T, "t")
   V, u = PolynomialRing(U, "u")

   f = (x + y + 2z^2 + 3t^3 + 5u^5 + 1)^10
   g = (u + t + 2z^2 + 3y^3 + 5x^5 + 1)^10
   
   f*g
end

function benchmark_resultant()
   R, x = FiniteField(7, 11, "x")
   S, y = PolynomialRing(R, "y")
   T = ResidueRing(S, y^3 + 3x*y + 1)
   U, z = PolynomialRing(T, "z")

   f = (3y^2 + y + x)*z^2 + ((x + 2)*y^2 + x + 1)*z + 4x*y + 3
   g = (7y^2 - y + 2x + 7)*z^2 + (3y^2 + 4x + 1)*z + (2x + 1)*y + 1

   s = f^12
   t = (s + g)^12
   
   resultant(s, t)
end

function benchmark_poly_nf_elem()
   R, x = CyclotomicField(20, "x")
   S, y = PolynomialRing(R, "y")

   f = (3x^7 + x^4 - 3x + 1)*y^3 + (2x^6-x^5+4x^4-x^3+x^2-1)*y +(-3x^7+2x^6-x^5+3x^3-2x^2+x)

   f^300
end

function benchmark_solve_poly()
   R, x = PolynomialRing(FlintZZ, "x")
   S, y = PolynomialRing(R, "y")
   M = MatrixSpace(S, 20, 20)()

   for i in 1:20
     for j in 1:20
       for k in 0:2
         M[i, j] = M[i,j] + y^k * (rand(-20:20) + rand(-20:20)*x + rand(-20:20)*x^2)
       end
     end
   end

   b = MatrixSpace(S, 20, 1)()

   for j in 1:20
     for k in 0:2
       b[j, 1] = b[j,1] + y^k * (rand(-20:20) + rand(-20:20)*x + rand(-20:20)*x^2)
     end
   end

   solve(M, b)
end

function benchmark_znz_det()
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

   determinant(M)
end


function benchmark_nf_det()
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

   determinant(M)
end

function benchmark_bernoulli()
   R, x = FlintQQ["x"]
   S, t = PowerSeriesRing(R, 1000, "t")

   u = t + O(t^1000)

   divexact((u*exp(x*u)), (exp(u)-1));
end

function benchmark_all()
   benchmark_fateman()
   benchmark_pearce()
   benchmark_resultant()
   benchmark_poly_nf_elem()
   benchmark_solve_poly()
   benchmark_znz_det()
   benchmark_nf_det()
   benchmark_bernoulli()
end

