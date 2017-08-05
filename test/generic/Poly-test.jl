function test_gen_poly_constructors()
   print("GenPoly.constructors...")
 
   R, x = ZZ["x"]
   S, y = R["y"]

   @test elem_type(S) == GenPoly{elem_type(R)}
   @test elem_type(GenPolyRing{elem_type(R)}) == GenPoly{elem_type(R)}
   @test parent_type(GenPoly{elem_type(R)}) == GenPolyRing{elem_type(R)}

   @test typeof(R) <: Nemo.Ring
   @test typeof(S) <: GenPolyRing

   @test isa(y, PolyElem)

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   @test typeof(S) <: GenPolyRing

   @test isa(y, PolyElem)

   T, z = PolynomialRing(S, "z")

   @test typeof(T) <: GenPolyRing

   @test isa(z, PolyElem)

   f = x^2 + y^3 + z + 1

   @test isa(f, PolyElem)

   g = S(2)

   @test isa(g, PolyElem)

   h = S(x^2 + 2x + 1)

   @test isa(h, PolyElem)

   j = T(x + 2)

   @test isa(j, PolyElem)

   k = S([x, x + 2, x^2 + 3x + 1])

   @test isa(k, PolyElem)

   l = S(k)

   @test isa(l, PolyElem)

   m = S([1, 2, 3])

   @test isa(m, PolyElem)

   n = S([ZZ(1), ZZ(2), ZZ(3)])

   @test isa(n, PolyElem)

   println("PASS")
end

function test_gen_poly_manipulation()
   print("GenPoly.manipulation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   @test iszero(zero(S))
   
   @test isone(one(S))

   @test isgen(gen(S))
   
   @test isunit(one(S))

   f = 2x*y + x^2 + 1

   @test lead(f) == 2x

   @test trail(2x*y + x^2) == x^2

   @test degree(f) == 1

   h = x*y^2 + (x + 1)*y + 3

   @test coeff(h, 2) == x

   @test length(h) == 3

   @test canonical_unit(-x*y + x + 1) == -1

   @test deepcopy(h) == h

   @test ismonomial(2*x*y^2)

   @test !ismonomial(2*x*y^2 + y + 1)

   println("PASS")
end

function test_gen_poly_binary_ops()
   print("GenPoly.binary_ops...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

#  Remi Imbach 31/07/17: begin
   maxIter = 100
   degmax  = 50
   bitsize = 20
   for iter = 1:maxIter
      f = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      g = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      h = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      @test (f - h) + (g + h) == (f + g)
      @test (f + g)*(f - g) == f*f - g*g
   end
   #the same over z/6Z
   T = ResidueRing(ZZ, 6)
   R,x = T["x"]
   S,y = R["y"]
   for iter = 1:maxIter
      f = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      g = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      h = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      @test (f - h) + (g + h) == (f + g)
      @test (f + g)*(f - g) == f*f - g*g
   end
#  Remi Imbach 31/07/17: end
   println("PASS")
end

function test_gen_poly_adhoc_binary()
   print("GenPoly.adhoc_binary...")
#  Remi Imbach 31/07/17: begin
   R,x = ZZ["x"]
   S,y = R["y"]
   U,z = S["z"]
   maxIter = 100
   degmax = 20
   bitsize = 20
   for iter = 1:maxIter
      f = rand(U, 0:degmax, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      c1 = rand(ZZ, -2^bitsize:2^bitsize)
      c2 = rand(ZZ, -2^bitsize:2^bitsize)
      @test c1*f - c2*f == (c1 - c2)*f
   end
   #the same over z/6Z
   T = ResidueRing(ZZ, 6)
   R,x = T["x"]
   S,y = R["y"]
   U,z = S["z"]
   for iter = 1:maxIter
      f = rand(U, 0:degmax, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      c1 = rand(T, -2^bitsize:2^bitsize)
      c2 = rand(T, -2^bitsize:2^bitsize)
      @test c1*f - c2*f == (c1 - c2)*f
   end
#  Remi Imbach 31/07/17: end
   println("PASS")
end

function test_gen_poly_comparison()
   print("GenPoly.comparison...")
#  Remi Imbach 31/07/17: begin
   R,x = ZZ["x"]
   S,y = R["y"]
   maxIter = 100
   degmax = 100
   bitsize = 32
   for iter = 1:maxIter
      f = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      g = f
      @test f==g
      @test isequal(f, g)
      g = f+(rand(ZZ, -2^bitsize:2^bitsize) +1)
      @test f!=g
   end
   #the same over z/6Z
   T = ResidueRing(ZZ, 6)
   R,x = T["x"]
   S,y = R["y"]
   for iter = 1:maxIter
      f = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      g = f
      @test f==g
      @test isequal(f, g)
      g = f+(rand(T, 0:4)+1)
      @test f!=g
   end
#  Remi Imbach 31/07/17: end
   println("PASS")
end

function test_gen_poly_adhoc_comparison()
   print("GenPoly.adhoc_comparison...")

#  Remi Imbach 31/07/17: begin
   R,x = ZZ["x"]
   S,y = R["y"]
   U,z = S["z"]
   maxIter = 100
   bitsize = 32
   for iter = 1:maxIter
      c1 = rand(ZZ, -2^bitsize:2^bitsize)
      @test U(c1) == c1 
      @test U(c1) != y + rand(ZZ, -2^bitsize:2^bitsize)*z
   end
   #the same over z/6Z
   T = ResidueRing(ZZ, 6)
   R,x = T["x"]
   S,y = R["y"]
   U,z = S["z"]
   for iter = 1:maxIter
      c1 = rand(T, -2^bitsize:2^bitsize)
      @test U(c1) == c1
      @test U(c1) != y + rand(T, 0:100)*z
   end
#  Remi Imbach 31/07/17: end
   println("PASS")
end

function test_gen_poly_unary_ops()
   print("GenPoly.unary_ops...")
#  Remi Imbach 31/07/17: begin
   R,x = ZZ["x"]
   S,y = R["y"]
   maxIter = 100
   degmax = 100
   bitsize = 32
   for iter = 1:maxIter
      f = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      c1 = rand(ZZ, -2^bitsize:2^bitsize)
      mf = -f
      @test -(S(c1)) == -c1
      @test degree(mf) == degree(f)
      for co = 1:(degree(f) + 1)
        @test coeff(mf, co) == -coeff(f,co)
      end
   end
   #the same over z/6Z
   T = ResidueRing(ZZ, 6)
   R,x = T["x"]
   S,y = R["y"]
   for iter = 1:maxIter
      f = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      c1 = rand(T, -2^bitsize:2^bitsize)
      mf = -f
      @test -(S(c1)) == -c1
      @test degree(mf) == degree(f)
      for co = 1:(degree(f) + 1)
        @test coeff(mf, co) == -coeff(f,co)
      end
   end
#  Remi Imbach 31/07/17: end
   println("PASS")
end

function test_gen_poly_truncation()
   print("GenPoly.truncation...")
#  Remi Imbach 31/07/17: begin
   R,x = ZZ["x"]
   S,y = R["y"]
   maxIter = 100
   degmax = 100
   bitsize = 32
   for iter = 1:maxIter
      f = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      t = rand(1:degmax)
      h = truncate(f, t)
      @test (degree(h) + 1) <= min(degree(f) + 1,t) #verify degree
      for co = 1:(degree(h) + 1) #verify each coeff
         @test coeff(h, co - 1)==coeff(f, co - 1)
      end
   end
#  the same over z/6Z
   T = ResidueRing(ZZ, 6)
   R,x = T["x"]
   S,y = R["y"]
   for iter = 1:maxIter
      f = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      t = rand(1:degmax)
      h = truncate(f, t)
      @test (degree(h) + 1) <= min(degree(f)+1,t) #verify degree
      for co = 1:(degree(h) + 1) #verify each coeff
         @test coeff(h, co - 1)==coeff(f, co - 1)
      end
   end
#  Remi Imbach 31/07/17: end
   println("PASS")
end

function test_gen_poly_reverse()
   print("GenPoly.reverse...")
   R,x = ZZ["x"]
   S,y = R["y"]
   maxIter = 100
   degmax = 50
   bitsize = 32
   for iter = 1:maxIter
      f = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      len = rand(1:(degmax + 1))
      f_rev = reverse(f, len)
      if f == 0
        @test f_rev == 0
      else #count the number of trailing 0 coeffs
        it = 0
        while coeff(f,it) == 0
            it = it+1
        end
        @test length(f_rev) == (len - it)
        for co = 1:length(f_rev)
          @test coeff(f_rev, co - 1) == coeff(f, len - co)
        end
      end
   end
#  the same over z/6Z
   T = ResidueRing(ZZ, 6)
   R,x = T["x"]
   S,y = R["y"]
   for iter = 1:maxIter
      f = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      len = rand(1:(degmax + 1))
      f_rev = reverse(f, len)
      if f == 0
        @test f_rev == 0
      else #count the number of trailing 0 coeffs
        it = 0
        while coeff(f,it) == 0
            it = it+1
        end
        @test length(f_rev) == (len - it)
        for co = 1:length(f_rev)
          @test coeff(f_rev, co - 1) == coeff(f, len - co)
        end
      end
   end
#  Remi Imbach 31/07/17: end
   println("PASS")
end

function test_gen_poly_shift()
   print("GenPoly.shift...")
#  Remi Imbach 31/07/17: begin
   R,x = ZZ["x"]
   S,y = R["y"]
   maxIter = 100
   degmax = 50
   bitsize = 32
   for iter = 1:maxIter
      f = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      s = rand(0:degmax)
      fsr = shift_right(f, s)
      fsl = shift_left(f, s)
      @test fsl == (y^s)*f
      f_rev = reverse(f, degree(f) + 1)
      f_rev_trun = truncate( f_rev, max(0, degree(f) + 1 - s) )
      f_rev_trun_rev = reverse( f_rev_trun, max(0, degree(f) + 1 - s) )
      @test fsr == f_rev_trun_rev
   end
#  the same over z/6Z
   T = ResidueRing(ZZ, 6)
   R,x = T["x"]
   S,y = R["y"]
   for iter = 1:maxIter
      f = rand(S, 0:degmax, 0:degmax, -2^bitsize:2^bitsize)
      s = rand(0:degmax)
      fsr = shift_right(f, s)
      fsl = shift_left(f, s)
      @test fsl == (y^s)*f
      f_rev = reverse(f, degree(f) + 1)
      f_rev_trun = truncate( f_rev, max(0, degree(f) + 1 - s) )
      f_rev_trun_rev = reverse( f_rev_trun, max(0, degree(f) + 1 - s) )
      @test fsr == f_rev_trun_rev
   end
#  Remi Imbach 31/07/17: end
   println("PASS")
end

function test_gen_poly_powering()
   print("GenPoly.powering...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   for iter = 1:10
      f = rand(S, 0:10, 0:10, -10:10)

      r2 = S(1)

      for expn = 0:10
         r1 = f^expn

         @test (f == 0 && expn == 0 && r1 == 0) || r1 == r2        

         r2 *= f
      end
   end

   for iter = 1:10
      n = rand(2:26)

      Zn = ResidueRing(ZZ, n)

      R, x = PolynomialRing(Zn, "x")
      S, y = PolynomialRing(R, "y")

      f = rand(S, 0:10, 0:10, 0:n - 1)

      r2 = S(1)

      for expn = 0:10
         r1 = f^expn

         @test (f == 0 && expn == 0 && r1 == 0) || r1 == r2         

         r2 *= f
      end
   end
   

   println("PASS")
end

function test_gen_poly_modular_arithmetic()
   print("GenPoly.modular_arithmetic...")

   for pow2n = 1:3
      R, x = PolynomialRing(ZZ, "x")
      S, y = PolynomialRing(QQ, "y")
      p = swinnerton_dyer(pow2n, x) # irreducible
      T = ResidueRing(S, S(p))
      U, z = PolynomialRing(T, "z")

      for iter = 1:10
         f = rand(U, 0:5, 0:2^pow2n, -10:10)
         g = rand(U, 0:5, 0:2^pow2n, -10:10)
         h = rand(U, 0:5, 0:2^pow2n, -10:10)
         k = U()
         while k == 0
            k = rand(U, 0:5, 0:2^pow2n, -10:10)
         end

         @test mulmod(mulmod(f, g, k), h, k) == mulmod(f, mulmod(g, h, k), k)
      end

      for iter = 1:10
         f = U()
         g = U()
         while f == 0 || g == 0 || gcd(f, g) != 1
            f = rand(U, 0:5, 0:2^pow2n, -10:10)
            g = rand(U, 0:5, 0:2^pow2n, -10:10)
         end

         @test mulmod(invmod(f, g), f, g) == mod(U(1), g)
      end

      for iter = 1:10
         f = rand(U, 0:5, 0:2^pow2n, -10:10)
         g = U()
         while g == 0
            g = rand(U, 0:5, 0:2^pow2n, -10:10)
         end
         p = mod(U(1), g)

         for expn = 0:5
            r = powmod(f, expn, g)

            @test (f == 0 && expn == 0 && r == 0) || r == p

            p = mulmod(p, f, g)
         end
      end
   end

   R, x = PolynomialRing(ZZ, "x")

   S = FractionField(R)

   T, y = PolynomialRing(S, "y")

   for iter = 1:10
      f = rand(T, 0:5, 0:3, -10:10)
      g = rand(T, 0:5, 0:3, -10:10)
      h = rand(T, 0:5, 0:3, -10:10)
      k = T()
      while k == 0
         k = rand(T, 0:5, 0:3, -10:10)
      end

      @test mulmod(mulmod(f, g, k), h, k) == mulmod(f, mulmod(g, h, k), k)
   end

   for iter = 1:10
      f = T()
      g = T()
      while f == 0 || g == 0 || gcd(f, g) != 1
         f = rand(T, 0:5, 0:3, -10:10)
         g = rand(T, 0:5, 0:3, -10:10)
      end

      @test mulmod(invmod(f, g), f, g) == mod(T(1), g)
   end

   for iter = 1:10
      f = rand(T, 0:5, 0:3, -10:10)
      g = T()
      while g == 0
         g = rand(T, 0:5, 0:3, -10:10)
      end
      p = mod(T(1), g)

      for expn = 0:5
         r = powmod(f, expn, g)

         @test (f == 0 && expn == 0 && r == 0) || r == p

         p = mulmod(p, f, g)
      end
   end

   println("PASS")
end

function test_gen_poly_exact_division()
   print("GenPoly.exact_division...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   for iter = 1:100
      f = rand(S, 0:10, 0:0, -100:100)
      g = S()
      while g == 0
         g = rand(S, 0:10, 0:0, -100:100)
      end

      @test divexact(f*g, g) == f
   end

   n = 23
   Zn = ResidueRing(ZZ, n)
   R, x = PolynomialRing(Zn, "x")
   S, y = PolynomialRing(R, "y")

   for iter = 1:100
      f = rand(S, 0:10, 0:10, 0:n - 1)
      g = S()
      while g == 0
         g = rand(S, 0:10, 0:10, 0:n - 1)
      end

      @test divexact(f*g, g) == f
   end

   println("PASS")
end

function test_gen_poly_adhoc_exact_division()
   print("GenPoly.adhoc_exact_division...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   for iter = 1:100
      f = rand(S, 0:10, 0:10, -100:100)
      g = R()
      while g == 0
         g = rand(R, 0:10, -100:100)
      end

      @test divexact(f*g, g) == f
   end

   n = 23
   Zn = ResidueRing(ZZ, n)
   R, x = PolynomialRing(Zn, "x")
   S, y = PolynomialRing(R, "y")

   for iter = 1:100
      f = rand(S, 0:10, 0:10, -100:100)
      g = R()
      while g == 0
         g = rand(R, 0:10, -100:100)
      end

      @test divexact(f*g, g) == f
   end

   println("PASS")
end

function test_gen_poly_euclidean_division()
   print("GenPoly.euclidean_division...")

   for pow2n = 1:3
      R, x = PolynomialRing(ZZ, "x")
      S, y = PolynomialRing(QQ, "y")
      p = swinnerton_dyer(pow2n, x) # irreducible
      T = ResidueRing(S, S(p))
      U, z = PolynomialRing(T, "z")

      for iter = 1:10
         f = rand(U, 0:5, 0:2^pow2n, -10:10)
         g = rand(U, 0:5, 0:2^pow2n, -10:10)
         h = U()
         while h == 0
            h = rand(U, 0:5, 0:2^pow2n, -10:10)
         end

         @test mod(f + g, h) == mod(f, h) + mod(g, h)
      end

      for iter = 1:10
         f = rand(U, 0:5, 0:2^pow2n, -10:10)
         g = U()
         while g == 0
            g = rand(U, 0:5, 0:2^pow2n, -10:10)
         end

         q, r = divrem(f, g)
         @test q*g + r == f

         @test mod(f, g) == r
      end
   end
 
   println("PASS")
end

function test_gen_poly_pseudodivision()
   print("GenPoly.pseudodivision...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   for iter = 1:100
      f = rand(S, 0:5, 0:0, -10:10)
      g = S()
      while g == 0
         g = rand(S, 0:5, 0:0, -10:10)
      end

      q, r = pseudodivrem(f, g)
      
      if length(f) < length(g)
         @test f == r && q == 0
      else
         @test q*g + r == f*lead(g)^(length(f) - length(g) + 1)
      end

      @test pseudorem(f, g) == r
   end

   println("PASS")
end

function test_gen_poly_content_primpart_gcd()
   print("GenPoly.content_primpart_gcd...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   for iter = 1:100
      f = rand(S, 0:10, 0:10, -100:100)

      g = R()
      while g == 0
         g = rand(R, 0:10, -100:100)
      end

      @test content(f*g) == divexact(g, canonical_unit(lead(g)))*content(f)

      @test primpart(f*g) == canonical_unit(lead(g))*primpart(f)
   end

   for iter = 1:20
      f = rand(S, 0:10, 0:10, -100:100)
      g = rand(S, 0:10, 0:10, -100:100)
      h = S()
      while h == 0
         h = rand(S, 0:10, 0:10, -100:100)
      end

      @test gcd(f*h, g*h) == divexact(h, canonical_unit(lead(h)))*gcd(f, g)
   end

   for pow2n = 1:3
      R, x = PolynomialRing(ZZ, "x")
      S, y = PolynomialRing(QQ, "y")
      p = swinnerton_dyer(pow2n, x) # irreducible
      T = ResidueRing(S, S(p))
      U, z = PolynomialRing(T, "z")

      for iter = 1:10
         f = U()
         g = U()
         while f == 0 || g == 0 || gcd(f, g) != 1
            f = rand(U, 0:5, 0:2^pow2n, -10:10)
            g = rand(U, 0:5, 0:2^pow2n, -10:10)
         end

         d, inv = gcdinv(f, g)

         @test d == gcd(f, g)

         @test mod(f*inv, g) == mod(U(1), g) 
      end
   end

   R, x = PolynomialRing(ZZ, "x")

   S = FractionField(R)

   T, y = PolynomialRing(S, "y")

   for iter = 1:100
      f = rand(T, 0:5, 0:3, -10:10)

      g = R()
      while g == 0
         g = rand(R, 0:3, -10:10)
      end

      @test content(f*g) == content(f)

      @test primpart(f*g) == primpart(f)*g
   end

   for iter = 1:20
      f = rand(T, 0:5, 0:3, -10:10)
      g = rand(T, 0:5, 0:3, -10:10)
      h = T()
      while h == 0
         h = rand(T, 0:5, 0:3, -10:10)
      end

      @test gcd(f*h, g*h) == inv(lead(h))*h*gcd(f, g)
   end

   for iter = 1:10
      f = T()
      g = T()
      while f == 0 || g == 0 || gcd(f, g) != 1
         f = rand(T, 0:5, 0:3, -10:10)
         g = rand(T, 0:5, 0:3, -10:10)
      end

      d, inv = gcdinv(f, g)

      @test d == gcd(f, g)

      @test mod(f*inv, g) == mod(T(1), g) 
   end

   println("PASS")
end

function test_gen_poly_evaluation()
   print("GenPoly.evaluation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x*y^2 + (x + 1)*y + 3

   @test evaluate(g, 3) == 12x + 6

   @test evaluate(g, fmpz(3)) == 12x + 6

   @test evaluate(g, f) == x^5+4*x^4+7*x^3+7*x^2+4*x+4

   for _ in 1:10
      f = rand(R, 0:4, -100:100)

      @test evaluate(g, f) == x*f^2 + (x + 1)*f + 3
      
      h = rand(S, 0:2, 0:2, -100:100)
      
      @test evaluate(h, f)^2 == evaluate(h^2, f)
   end

   println("PASS")
end

function test_gen_poly_composition()
   print("GenPoly.composition...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3
   g = (x + 1)*y + (x^3 + 2x + 2)

   @test compose(f, g) == (x^3+2*x^2+x)*y^2+(2*x^5+2*x^4+4*x^3+9*x^2+6*x+1)*y+(x^7+4*x^5+5*x^4+5*x^3+10*x^2+8*x+5)

   for d in 1:10
      g = rand(S, 0:d, 0:2, -100:100)

      @test compose(f, g) == x*g^2 + (x + 1)*g + 3

      h = rand(S, 0:d, 0:2, -100:100)

      @test compose(h, g)^2 == compose(h^2, g)

      k = rand(S, 0:d, 0:2, -100:100)

      @test compose(g, compose(h, k)) == compose(compose(g, h), k)
   end

   println("PASS")
end

function test_gen_poly_derivative()
   print("GenPoly.derivative...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   h = x*y^2 + (x + 1)*y + 3

   @test derivative(h) == 2x*y + x + 1

   for _ in 1:10
      f = rand(S, 0:4, 0:4, -100:100)
      g = rand(S, 0:4, 0:4, -100:100)

      @test derivative(f + g) == derivative(g) + derivative(f)

      @test derivative(g*f) == derivative(g)*f + derivative(f)*g
   end

   println("PASS")
end

function test_gen_poly_integral()
   print("GenPoly.integral...")

   R, x = PolynomialRing(QQ, "x")
   S = ResidueRing(R, x^3 + 3x + 1)
   T, y = PolynomialRing(S, "y")

   f = (x^2 + 2x + 1)*y^2 + (x + 1)*y - 2x + 4

   @test integral(f) == (fmpz(1)//3*x^2 + fmpz(2)//3*x + fmpz(1)//3)*y^3+(fmpz(1)//2*x+fmpz(1)//2)*y^2+(-2*x+4)*y

   for _ in 1:10
      f = rand(T, 0:2, 0:2, -100:100)

      @test derivative(integral(f)) == f

      g = rand(T, 0:2, 0:2, -100:100)

      @test integral(f + g) == integral(g) + integral(f)
      @test integral(f)*integral(g) == integral(integral(f)*g + integral(g)*f)
   end

   println("PASS")
end

function test_gen_poly_resultant()
   print("GenPoly.resultant...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = 3x*y^2 + (x + 1)*y + 3
   g = 6(x + 1)*y + (x^3 + 2x + 2)

   @test resultant(f, g) == 3*x^7+6*x^5-6*x^3+96*x^2+192*x+96

   for _ in 1:10
      f = rand(S, 0:2, 0:2, -100:100)
      g = rand(S, 0:2, 0:2, -100:100)
      h = rand(S, 0:2, 0:2, -100:100)

      @test resultant(f*g, h) == resultant(f, h) * resultant(g, h)
      @test resultant(f, g*h) == resultant(f, g) * resultant(f, h)
   end

   println("PASS")
end

function test_gen_poly_discriminant()
   print("GenPoly.discriminant...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3

   @test discriminant(f) == x^2-10*x+1

   println("PASS")
end

function test_gen_poly_gcdx()
   print("GenPoly.gcdx...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = 3x*y^2 + (x + 1)*y + 3
   g = 6(x + 1)*y + (x^3 + 2x + 2)

   @test gcdx(f, g) == (3*x^7+6*x^5-6*x^3+96*x^2+192*x+96, (36*x^2+72*x+36), (-18*x^2-18*x)*y+(3*x^4-6*x-6))

   for _ in 1:100
      f = S()
      g = S()
      while length(f) < 2 || length(g) < 2
         f = rand(S, 2:2, 0:2, 1:100)
         g = rand(S, 2:2, 0:2, -100:-1)
      end
      r, u, v = gcdx(f, g)

      @test u*f + v*g == r

      h = S()
      while h == 0
         h = rand(S, 0:2, 0:2, -100:100)
      end
      r, u, v = gcdx(f*h, g*h)

      @test (u*f + v*g)*h == r
   end

   R, x = PolynomialRing(QQ, "x")
   T = ResidueRing(R, (x^3 + 3x - 2)*(x - 1))
   S, z = PolynomialRing(T, "z")

   for _ in 1:10
      f = S()
      g = S()
      while length(f) < 2 || length(g) < 2
         f = rand(S, 2:2, 0:2, 1:100)
         g = rand(S, 2:2, 0:2, -100:-1)
      end
      r, u, v = gcdx(f, g)

      @test u*f + v*g == r

      h = S()
      while h == 0
         h = rand(S, 0:2, 0:2, -100:100)
      end
      r, u, v = gcdx(f*h, g*h)

      @test (u*f + v*g)*h == r
   end

   println("PASS")
end

function test_gen_poly_newton_representation()
   print("GenPoly.newton_representation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   for _ in 1:10
      f = rand(S, 2:2, 0:2, 1:100)

      g = deepcopy(f)
      roots = [R(1), R(2), R(3)]
      monomial_to_newton!(g.coeffs, roots)
      newton_to_monomial!(g.coeffs, roots)

      @test f == g
   end

   println("PASS")
end

function test_gen_poly_interpolation()
   print("GenPoly.interpolation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   xs = [R(1), R(2), R(3), R(4)]
   ys = [R(1), R(4), R(9), R(16)]

   f = interpolate(S, xs, ys)

   @test f == y^2

   for _ in 1:10
      ys = [rand(R, 0:2, -100:100) for i in 1:4]
      f = interpolate(S, xs, ys)

      @test all(evaluate(f, a) == b for (a,b) in zip(xs, ys))
   end

   println("PASS")
end

function test_gen_poly_special()
   print("GenPoly.special...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   @test chebyshev_t(20, y) == 524288*y^20-2621440*y^18+5570560*y^16-6553600*y^14+4659200*y^12-2050048*y^10+549120*y^8-84480*y^6+6600*y^4-200*y^2+1

   @test chebyshev_u(15, y) == 32768*y^15-114688*y^13+159744*y^11-112640*y^9+42240*y^7-8064*y^5+672*y^3-16*y

   for n in 10:20
      T = chebyshev_t(n, y)
      dT = derivative(T)
      ddT = derivative(dT)

      @test (1 - y^2)*ddT + n^2*T == y*dT

      U = chebyshev_u(n - 1, y)
      dU = derivative(U)
      ddU = derivative(dU)

      @test (1 - y^2)*ddU + (n-1)*(n+1)*U == 3*y*dU

      @test T^2 == 1 + (y^2 - 1)*U^2
   end

   println("PASS")
end

function test_gen_poly_mul_karatsuba()
   print("GenPoly.mul_karatsuba...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")
   
   f = x + y + 2z^2 + 1
   
   @test mul_karatsuba(f^10, f^10) == mul_classical(f^10, f^10)
   @test mul_karatsuba(f^10, f^30) == mul_classical(f^10, f^30)

   println("PASS")
end

function test_gen_poly_mul_ks()
   print("GenPoly.mul_ks...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")
   
   f = x + y + 2z^2 + 1
   
   @test mul_ks(f^10, f^10) == mul_classical(f^10, f^10)
   @test mul_ks(f^10, f^30) == mul_classical(f^10, f^30)

   println("PASS")
end

function test_gen_poly_remove_valuation()
   print("GenPoly.remove_valuation...")

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3-2, "a");
   S, y = PolynomialRing(K, "y")

   f = 3a*y^2 + (a + 1)*y + 3
   g = 6(a + 1)*y + (a^3 + 2x + 2)

   v, q = remove(f^3*g^4, g)

   @test valuation(f^3*g^4, g) == 4
   @test q == f^3
   @test v == 4

   v, q = divides(f*g, f)

   @test v
   @test q == g

   v, q = divides(f*g + 1, f)

   @test !v
   
   println("PASS")
end

function test_gen_poly_generic_eval()
   print("GenPoly.generic_eval...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = 3x*y^2 + (x + 1)*y + 3
   g = 6(x + 1)*y + (x^3 + 2x + 2)

   @test f(g) == (108*x^3+216*x^2+108*x)*y^2+(36*x^5+36*x^4+72*x^3+150*x^2+84*x+6)*y+(3*x^7+12*x^5+13*x^4+13*x^3+26*x^2+16*x+5)

   @test f(x + 1) == 3*x^3+7*x^2+5*x+4

   @test f(123) == 45510*x+126

   @test f(fmpz(123)) == 45510*x + 126

   for _ in 1:10
      f = rand(S, 0:2, 0:2, -100:100)
      g = rand(S, 0:2, 0:2, -100:100)
      h = rand(S, 0:2, 0:2, -100:100)

      @test f(g(h)) == f(g)(h)
   end

   R, x = PolynomialRing(ZZ, "x")
   T, y = FiniteField(103, 1, "y")

   f = x^5 + 3x^3 + 2x^2 + x + 1

   @test f(T(13)) == 20

   f = x
   b = a = T(13)
   for i in 1:5
      g = x^2 + rand(R, 0:1, -1:1)
      f = g(f)
      b = g(b)

      @test b == f(a)
   end

   println("PASS")
end

function test_gen_poly()
   test_gen_poly_constructors()
   test_gen_poly_manipulation()
   test_gen_poly_binary_ops()
   test_gen_poly_adhoc_binary()
   test_gen_poly_comparison()
   test_gen_poly_adhoc_comparison()
   test_gen_poly_unary_ops()
   test_gen_poly_truncation()
   test_gen_poly_reverse()
   test_gen_poly_shift()
   test_gen_poly_powering()
   test_gen_poly_modular_arithmetic()
   test_gen_poly_exact_division()
   test_gen_poly_adhoc_exact_division()
   test_gen_poly_euclidean_division()
   test_gen_poly_pseudodivision()
   test_gen_poly_content_primpart_gcd()
   test_gen_poly_evaluation()
   test_gen_poly_composition()
   test_gen_poly_derivative()
   test_gen_poly_integral()
   test_gen_poly_resultant()
   test_gen_poly_discriminant()
   test_gen_poly_gcdx()
   test_gen_poly_newton_representation()
   test_gen_poly_interpolation()
   test_gen_poly_special()
   test_gen_poly_mul_karatsuba()
   test_gen_poly_mul_ks()
   if VERSION >= v"0.5.0-dev+3171"
      test_gen_poly_generic_eval()
   end
   test_gen_poly_remove_valuation()
   

   println("")
end
