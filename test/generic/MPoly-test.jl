function test_gen_mpoly_constructors()
   print("GenMPoly.constructors...")

   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      @test elem_type(S) == GenMPoly{elem_type(R)}
      @test elem_type(GenMPolyRing{elem_type(R)}) == GenMPoly{elem_type(R)}
      @test parent_type(GenMPoly{elem_type(R)}) == GenMPolyRing{elem_type(R)}

      @test typeof(S) <: GenMPolyRing

      isa(vars(S), Array{Symbol, 1})

      for j = 1:num_vars
         @test isa(varlist[j], PolyElem)
         @test isa(gens(S)[j], PolyElem)
      end

      f =  rand(S, 0:5, 0:100, 0:0, -100:100)

      @test isa(f, PolyElem)

      @test isa(S(2), PolyElem)

      @test isa(S(R(2)), PolyElem)

      @test isa(S(f), PolyElem)

      V = [R(rand(-100:100)) for i in 1:5]

      W0 = UInt[rand(0:100) for i in 1:5*num_vars]
      W = reshape(W0, num_vars, 5)

      @test isa(S(V, W), PolyElem)
   end

   println("PASS")
end   

function test_gen_mpoly_manipulation()
   print("GenMPoly.manipulation...")

   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)
      g = gens(S)

      @test !isgen(S(1))

      for i = 1:num_vars
         @test isgen(varlist[i])
         @test isgen(g[i])
         @test !isgen(g[i] + 1)
      end

      f = rand(S, 0:5, 0:100, 0:0, -100:100)

      @test f == deepcopy(f)

      if length(f) > 0
         @test isa(coeff(f, rand(1:length(f)) - 1), elem_type(R))
      end

      max_degs, biggest = max_degrees(f)
      deg = isdegree(ordering(S))
      rev = isreverse(ordering(S))

      if deg
         @test max_degs[1] <= 100*num_vars
         @test max_degs[1] == biggest
         @test S.N == num_vars + 1
      else
         @test S.N == num_vars
      end

      if rev
         @test S.N == num_vars + 1
      end

      for j = Int(deg) + 1:num_vars + Int(deg)
         @test max_degs[j] <= 100
         @test max_degs[j] <= biggest
      end

      @test ord == ordering(S)

      @test nvars(f) == num_vars

      @test isone(one(S))

      @test iszero(zero(S))

      @test isconstant(S(rand(-100:100)))
      @test isconstant(S(zero(S)))

      g = rand(S, 1:1, 0:100, 0:0, 1:100)
      h = rand(S, 1:1, 0:100, 0:0, 1:1)

      @test isterm(g)
      @test ismonomial(h)
   end

   println("PASS")
end   

function test_gen_mpoly_unary_ops()
   print("GenMPoly.unary_ops...")

   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      for iter = 1:10
         f = rand(S, 0:5, 0:100, 0:0, -100:100)

         @test f == -(-f)
      end
   end

   println("PASS")
end   

function test_gen_mpoly_binary_ops()
   print("GenMPoly.binary_ops...")

   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      for iter = 1:10
         f = rand(S, 0:5, 0:100, 0:0, -100:100)
         g = rand(S, 0:5, 0:100, 0:0, -100:100)
         h = rand(S, 0:5, 0:100, 0:0, -100:100)

         @test f + g == g + f
         @test f - g == -(g - f)
         @test f*g == g*f
         @test f*g + f*h == f*(g + h)
         @test f*g - f*h == f*(g - h)

         @test f*g == Nemo.mul_classical(f, g)
      end
   end

   println("PASS")
end   

function test_gen_mpoly_adhoc_binary()
   print("GenMPoly.adhoc_binary...")

   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      for iter = 1:100
         f = rand(S, 0:5, 0:100, 0:0, -100:100)

         d1 = rand(-20:20)
         d2 = rand(-20:20)
         g1 = rand(R, 0:2, -10:10)
         g2 = rand(R, 0:2, -10:10)

         @test f*d1 + f*d2 == (d1 + d2)*f
         @test f*ZZ(d1) + f*ZZ(d2) == (ZZ(d1) + ZZ(d2))*f
         @test f*g1 + f*g2 == (g1 + g2)*f
         
         @test f + d1 + d2 == d1 + d2 + f
         @test f + ZZ(d1) + ZZ(d2) == ZZ(d1) + ZZ(d2) + f
         @test f + g1 + g2 == g1 + g2 + f
         
         @test f - d1 - d2 == -((d1 + d2) - f)
         @test f - ZZ(d1) - ZZ(d2) == -((ZZ(d1) + ZZ(d2)) - f)
         @test f - g1 - g2 == -((g1 + g2) - f)

         @test f + d1 - d1 == f
         @test f + ZZ(d1) - ZZ(d1) == f
         @test f + g1 - g1 == f
      end
   end

   println("PASS")
end   

function test_gen_mpoly_adhoc_comparison()
   print("GenMPoly.adhoc_comparison...")

   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      for iter = 1:100
         d = rand(-100:100)
         g = rand(R, 0:2, -10:10)
         
         @test S(d) == d
         @test d == S(d)
         @test S(d) == ZZ(d)
         @test ZZ(d) == S(d)
         @test S(g) == g
         @test g == S(g)
      end
   end

   println("PASS")
end   

function test_gen_mpoly_powering()
   print("GenMPoly.powering...")

   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      for iter = 1:10
         f = rand(S, 0:5, 0:100, 0:0, -100:100)

         expn = rand(0:10)

         r = S(1)
         for i = 1:expn
            r *= f
         end

         @test (f == 0 && expn == 0 && f^expn == 0) || f^expn == r 
      end
   end

   println("PASS")
end   

function test_gen_mpoly_divides()
   print("GenMPoly.divides...")

   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      for iter = 1:10
         f = S(0)
         while iszero(f)
            f = rand(S, 0:5, 0:100, 0:0, -100:100)
         end
         g = rand(S, 0:5, 0:100, 0:0, -100:100)

         p = f*g

         flag, q = divides(p, f)

         @test flag == true
         @test q == g

         q = divexact(p, f)

         @test q == g
      end
   end

   println("PASS")
end   

function test_gen_mpoly()
   test_gen_mpoly_constructors()
   test_gen_mpoly_manipulation()
   test_gen_mpoly_unary_ops()
   test_gen_mpoly_binary_ops()
   test_gen_mpoly_adhoc_binary()
   test_gen_mpoly_adhoc_comparison()
   test_gen_mpoly_powering()
   test_gen_mpoly_divides()

   println("")
end
