function test_gen_mpoly_constructors()
   print("Generic.MPoly.constructors...")

   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      @test elem_type(S) == Generic.MPoly{elem_type(R)}
      @test elem_type(Generic.MPolyRing{elem_type(R)}) == Generic.MPoly{elem_type(R)}
      @test parent_type(Generic.MPoly{elem_type(R)}) == Generic.MPolyRing{elem_type(R)}

      @test typeof(S) <: Generic.MPolyRing

      isa(symbols(S), Array{Symbol, 1})

      for j = 1:num_vars
         @test isa(varlist[j], MPolyElem)
         @test isa(gens(S)[j], MPolyElem)
      end

      f =  rand(S, 0:5, 0:100, 0:0, -100:100)

      @test isa(f, MPolyElem)

      @test isa(S(2), MPolyElem)

      @test isa(S(R(2)), MPolyElem)

      @test isa(S(f), MPolyElem)

      V = [R(rand(-100:100)) for i in 1:5]

      W0 = UInt[rand(0:100) for i in 1:5*num_vars]
      W = reshape(W0, num_vars, 5)

      @test isa(S(V, W), MPolyElem)

      W1 = [[rand(0:100) for i in 1:num_vars] for j in 1:5]

      f1 = S(V, W1)

      @test isa(f1, MPolyElem)

      f2 = S()
      fit!(f2, 5)

      for i = 1:5
         f2 = set_exponent_vector!(f2, i, W1[i])
         f2 = setcoeff!(f2, i, V[i])
      end
      f2 = sort_terms!(f2)
      f2 = combine_like_terms!(f2)

      @test f1 == f2

      C = MPolyBuildCtx(S)

      for i = 1:5
         push_term!(C, V[i], W1[i])
      end
      f3 = finish(C)

      @test f1 == f3

      _, varlist = PolynomialRing(QQ, var_names)
      y = varlist[1]
      @test x in [x, y]
      @test x in [y, x]
      @test !(x in [y])
      @test x in keys(Dict(x => 1))
      @test !(y in keys(Dict(x => 1)))
   end

   println("PASS")
end

function test_gen_mpoly_geobuckets()
   print("Generic.MPoly.geobuckets...")

   R, x = ZZ["y"]
   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)
      g = gens(S)

      C = Generic.geobucket(S)
      for j in 1:num_vars
        push!(C, sum(g[j]^i for i in 1:64))
      end
      @test finish(C) == sum(sum(g[j]^i for i in 1:64) for j in 1:num_vars)
   end
   println("PASS")
end

function test_gen_mpoly_manipulation()
   print("Generic.MPoly.manipulation...")

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
         @test gen(S, i) == g[i]
         @test var_index(gen(S, i)) == i
      end

      nv = rand(1:num_vars)
      f = S(2)

      @test length(vars(f)) == 0

      f = 3*g[nv] + 2

      @test length(vars(f)) == 1 && vars(f)[1] == g[nv]

      f = 2*g[1]*g[num_vars] + 12

      if num_vars == 1
         @test length(vars(f)) == 1 && vars(f)[1] == g[1]
      else
         @test length(vars(f)) == 2 && vars(f)[1] == g[1] && vars(f)[2] == g[num_vars]
      end

      f = rand(S, 0:5, 0:100, 0:0, -100:100)

      @test f == deepcopy(f)

      @test hash(f) == hash(deepcopy(f))

      if length(f) > 0
         @test isa(coeff(f, rand(1:length(f))), elem_type(R))
      end

      @test length(f) == length(coeffs(f))
      @test length(f) == length(monomials(f))
      @test length(f) == length(terms(f))

      m = one(S)
      r = zero(S)
      r2 = zero(S)
      for i = 1:length(f)
         m = monomial!(m, f, i)
         @test m == monomial(f, i)
         @test term(f, i) == coeff(f, i)*monomial(f, i)
         r += coeff(f, i)*monomial(f, i)
         r2 += coeff(f, monomial(f, i))*monomial(f, i)
      end
      @test r == f
      @test r2 == f

      for i = 1:length(f)
         i1 = rand(1:length(f))
         i2 = rand(1:length(f))
         @test (i1 < i2) == (monomial(f, i1) > monomial(f, i2))
         @test (i1 > i2) == (monomial(f, i1) < monomial(f, i2))
         @test (i1 == i2) == (monomial(f, i1) == monomial(f, i2))
      end

      max_degs, biggest = max_fields(f)
      deg = isdegree(ordering(S))
      rev = isreverse(ordering(S))

      if deg
         @test max_degs[num_vars + 1] <= 100*num_vars
         @test max_degs[num_vars + 1] == biggest
         @test S.N == num_vars + 1
      else
         @test S.N == num_vars
      end

      if rev
         @test S.N == num_vars + 1
      end

      degs = degrees(f)

      for j = 1:num_vars
         @test degs[j] == degree(f, j)
         if rev
            @test max_degs[j] == degs[j] || (degs[j] == -1 && max_degs[j] == 0)
         else
            @test max_degs[j] == degs[num_vars - j + 1] ||
                  (degs[num_vars - j + 1] == -1 && max_degs[j] == 0) 
         end

         @test max_degs[j] <= biggest
      end

      @test ord == ordering(S)

      @test nvars(parent(f)) == num_vars

      @test isone(one(S))

      @test iszero(zero(S))

      @test isunit(S(1))
      @test !isunit(S(0))
      @test !isunit(gen(S, 1))

      @test isconstant(S(rand(-100:100)))
      @test isconstant(S(zero(S)))

      g = rand(S, 1:1, 0:100, 0:0, 1:100)
      h = rand(S, 1:1, 0:100, 0:0, 1:1)

      @test isterm(g)
      @test ismonomial(h)
   end

   println("PASS")
end

function test_gen_mpoly_multivariate_coeff()
   print("Generic.MPoly.multivariate_coeff...")

   R = ZZ

   for iter = 1:5
      ord = rand_ordering()

      S, (x, y, z) = PolynomialRing(R, ["x", "y", "z"]; ordering=ord)

      f = -8*x^5*y^3*z^5+9*x^5*y^2*z^3-8*x^4*y^5*z^4-10*x^4*y^3*z^2+8*x^3*y^2*z-10*x*y^3*
z^4-4*x*y-10*x*z^2+8*y^2*z^5-9*y^2*z^3

      @test coeff(f, [1], [1]) == -10*y^3*z^4-4*y-10*z^2
      @test coeff(f, [2, 3], [3, 2]) == -10*x^4
      @test coeff(f, [1, 3], [4, 5]) == 0

      @test coeff(f, [x], [1]) == -10*y^3*z^4-4*y-10*z^2
      @test coeff(f, [y, z], [3, 2]) == -10*x^4
      @test coeff(f, [x, z], [4, 5]) == 0
   end

   println("PASS")
end

function test_gen_mpoly_total_degree()
   print("Generic.MPoly.total_degree...")
   max = 50
   for nvars = 1:10
      var_names = ["x$j" for j in 1:nvars]
      for nterms = 1:10
         exps = Array{Int,2}(round.(rand(nvars, nterms) .* max))
         degrees = []
         for nord = 1:20
            ord = rand_ordering()
            S, varlist = PolynomialRing(ZZ, var_names, ordering = ord)
            p = zero(S)
            for j = 1:nterms
               p += prod(varlist[i]^exps[i,j] for i = 1:nvars)
            end
            push!(degrees, total_degree(p))
         end
         @test length(Set(degrees)) == 1
      end
   end
   println("PASS")
end

function test_gen_mpoly_unary_ops()
   print("Generic.MPoly.unary_ops...")

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
   print("Generic.MPoly.binary_ops...")

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

         @test f*g == AbstractAlgebra.mul_classical(f, g)
      end
   end

   println("PASS")
end

function test_gen_mpoly_adhoc_binary()
   print("Generic.MPoly.adhoc_binary...")

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
         @test f*BigInt(d1) + f*BigInt(d2) == (BigInt(d1) + BigInt(d2))*f
         @test f*g1 + f*g2 == (g1 + g2)*f

         @test f + d1 + d2 == d1 + d2 + f
         @test f + BigInt(d1) + BigInt(d2) == BigInt(d1) + BigInt(d2) + f
         @test f + g1 + g2 == g1 + g2 + f

         @test f - d1 - d2 == -((d1 + d2) - f)
         @test f - BigInt(d1) - BigInt(d2) == -((BigInt(d1) + BigInt(d2)) - f)
         @test f - g1 - g2 == -((g1 + g2) - f)

         @test f + d1 - d1 == f
         @test f + BigInt(d1) - BigInt(d1) == f
         @test f + g1 - g1 == f

         if !iszero(d1)
           @test divexact(d1 * f, d1) == f
           @test divexact(d1 * f, BigInt(d1)) == f
         end
      end
   end

   println("PASS")
end

function test_gen_mpoly_adhoc_comparison()
   print("Generic.MPoly.adhoc_comparison...")

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
         @test S(d) == BigInt(d)
         @test BigInt(d) == S(d)
         @test S(g) == g
         @test g == S(g)
      end
   end

   println("PASS")
end

function test_gen_mpoly_powering()
   print("Generic.MPoly.powering...")

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
   print("Generic.MPoly.divides...")

   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      for iter = 1:10
         f = S(0)
         f = rand(S, 0:5, 0:100, 0:0, -100:100)
         g = rand(S, 0:5, 0:100, 0:0, -100:100)

         p = f*g

         flag, q = divides(p, f)
         flag2, q2 = divides(f, p)

         @test flag == true

         @test q * f == p

         q1 = divexact(p, f)

         @test q1 * f == p

         if !iszero(p)
           @test q1 == g
         end
      end
   end

   println("PASS")
end

function test_gen_mpoly_euclidean_division()
   print("Generic.MPoly.euclidean_division...")

   R, x = QQ["y"]

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

         q1, r = divrem(p, f)
         q2 = div(p, f)

         @test q1 == g
         @test q2 == g
         @test f*q1 + r == p

         q3, r3 = divrem(g, f)
         q4 = div(g, f)
         flag, q5 = divides(g, f)

         @test q3*f + r3 == g
         @test q3 == q4
         @test (r3 == 0 && flag == true && q5 == q3) || (r3 != 0 && flag == false)
         
      end
      
      S, varlist = PolynomialRing(QQ, var_names, ordering = ord)
      v = varlist[1+Int(round(rand() * (num_vars-1)))]
      @test divrem(v, 2*v) == (1//2, 0)
   end

   println("PASS")
end

function test_gen_mpoly_ideal_reduction()
   print("Generic.MPoly.ideal_reduction...")

   R, x = QQ["y"]

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

         q1, r = divrem(p, [f])

         @test q1[1] == g
         @test r == 0
      end

      for iter = 1:10
         num = rand(1:5)

         V = Array{elem_type(S)}(undef, num)

         for i = 1:num
            V[i] = S(0)
            while iszero(V[i])
               V[i] = rand(S, 0:5, 0:100, 0:0, -100:100)
            end
         end
         g = rand(S, 0:5, 0:100, 0:0, -100:100)

         q, r = divrem(g, V)

         p = r
         for i = 1:num
            p += q[i]*V[i]
         end

         @test p == g
      end
   end

   println("PASS")
end

function test_gen_mpoly_deflation()
   print("Generic.MPoly.deflation...")

   for num_vars = 1:4
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()
      S, varlist = PolynomialRing(ZZ, var_names, ordering = ord)

      for iter = 1:10
         f = rand(S, 0:4, 0:5, -10:10)
         shift = [rand(0:10) for i in 1:num_vars]
         defl = [rand(1:10) for i in 1:num_vars]
         f = inflate(f, shift, defl)

         s, d = deflation(f)
         g = deflate(f, s, d)
         h = inflate(g, s, d)

         @test h == f
      end
   end

   println("PASS")
end

function test_gen_mpoly_gcd()
   print("Generic.MPoly.gcd...")

   for num_vars = 1:4
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()
      S, varlist = PolynomialRing(ZZ, var_names, ordering = ord)

      for iter = 1:10
         f = rand(S, 0:4, 0:5, -10:10)
         g = rand(S, 0:4, 0:5, -10:10)
         h = rand(S, 0:4, 0:5, -10:10)

         g1 = gcd(f, g)
         g2 = gcd(f*h, g*h)

         @test g2 == g1*h || g2 == -g1*h
      end
   end

   println("PASS")
end

function test_gen_mpoly_lcm()
   print("Generic.MPoly.lcm...")

   for num_vars = 1:4
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()
      S, varlist = PolynomialRing(ZZ, var_names, ordering = ord)

      for iter = 1:10
         f = rand(S, 0:4, 0:5, -10:10)
         g = rand(S, 0:4, 0:5, -10:10)
         h = rand(S, 0:4, 0:5, -10:10)

         l1 = lcm(f, g)
         l2 = lcm(f*h, g*h)

         @test l2 == l1*h || l2 == -l1*h
         @test divides(l1, f)[1]
         @test divides(l1, g)[1]
         @test divides(l2, f)[1]
         @test divides(l2, g)[1]
         @test divides(l2, h)[1]
      end
   end

   println("PASS")
end

function test_gen_mpoly_evaluation()
   print("Generic.MPoly.evaluation...")

   R, x = ZZ["x"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      for iter = 1:50
         f = rand(S, 0:5, 0:100, 0:0, -100:100)
         g = rand(S, 0:5, 0:100, 0:0, -100:100)

         V1 = [rand(-10:10) for i in 1:num_vars]

         r1 = evaluate(f, V1)
         r2 = evaluate(g, V1)
         r3 = evaluate(f + g, V1)

         @test r3 == r1 + r2

         V2 = [BigInt(rand(-10:10)) for i in 1:num_vars]

         r1 = evaluate(f, V2)
         r2 = evaluate(g, V2)
         r3 = evaluate(f + g, V2)

         @test r3 == r1 + r2

         V3 = [R(rand(-10:10)) for i in 1:num_vars]

         r1 = evaluate(f, V3)
         r2 = evaluate(g, V3)
         r3 = evaluate(f + g, V3)

         @test r3 == r1 + r2
      end
   end

   R1, z = ZZ["z"]
   R, x = R1["x"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      for iter = 1:50
         f = rand(S, 0:5, 0:100, 0:0, 0:0, -100:100)
         g = rand(S, 0:5, 0:100, 0:0, 0:0, -100:100)

         V1 = [rand(R1, 0:0, -10:10) for i in 1:num_vars]

         r1 = evaluate(f, V1)
         r2 = evaluate(g, V1)
         r3 = evaluate(f + g, V1)

         @test r3 == r1 + r2
      end
   end

   R, x = ZZ["x"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      for iter = 1:10
         eval_num = rand(0:num_vars)
         V = Int[] # random list of variable indices
         Vflag = [false for i in 1:num_vars] # whether each variable is in V
         Vval = Int[] # value substituted for variables in V
         Vals = [0 for i in 1:num_vars] # value subst. for each variable of pol

         for i = 1:eval_num
            v = rand(1:num_vars)
            while Vflag[v]
               v = rand(1:num_vars)
            end
            push!(V, v)
            Vflag[v] = true
            c = rand(-10:10)
            push!(Vval, c)
            Vals[v] = c
         end

         W = Int[] # remaining variables
         Wval = Int[] # values for those variables

         for v = 1:num_vars
            if !Vflag[v] # no value for this var yet
               push!(W, v)
               c = rand(-10:10)
               push!(Wval, c)
               Vals[v] = c
            end
         end

         f = rand(S, 0:5, 0:100, 0:0, -100:100)
      
         f1 = evaluate(f, V, Vval)
         f2 = evaluate(f1, W, Wval)
         
         r = evaluate(f, Vals)

         @test (length(f2) == 0 && r == 0) ||
               (length(f2) == 1 && r == coeff(f2, 1))
      end

      for iter = 1:10
         eval_num = rand(0:num_vars)
         V = Int[] # random list of variable indices
         Vflag = [false for i in 1:num_vars] # whether each variable is in V
         Vval = Array{elem_type(R), 1}(undef, 0) # value substituted for variables in V
         Vals = [R(0) for i in 1:num_vars] # value subst. for each variable of pol

         for i = 1:eval_num
            v = rand(1:num_vars)
            while Vflag[v]
               v = rand(1:num_vars)
            end
            push!(V, v)
            Vflag[v] = true
            c = R(rand(-10:10))
            push!(Vval, c)
            Vals[v] = c
         end

         W = Int[] # remaining variables
         Wval = Array{elem_type(R), 1}(undef, 0) # values for those variables

         for v = 1:num_vars
            if !Vflag[v] # no value for this var yet
               push!(W, v)
               c = R(rand(-10:10))
               push!(Wval, c)
               Vals[v] = c
            end
         end

         f = rand(S, 0:5, 0:100, 0:0, -100:100)

         f1 = evaluate(f, V, Vval)
         f2 = evaluate(f1, W, Wval)

         r = evaluate(f, Vals)

         @test (length(f2) == 0 && r == 0) ||
               (length(f2) == 1 && r == coeff(f2, 1))
      end
   end

   R = ZZ
   T = MatrixAlgebra(R, 2)

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      for iter = 1:50
         f = rand(S, 0:5, 0:100, -100:100)
         g = rand(S, 0:5, 0:100, -100:100)

         V1 = [rand(-10:10) for i in 1:num_vars]

         r1 = f(V1...)
         r2 = g(V1...)
         r3 = (f + g)(V1...)

         @test r3 == r1 + r2
      end

      for iter = 1:50
         f = rand(S, 0:5, 0:100, -100:100)
         g = rand(S, 0:5, 0:100, -100:100)

         V1 = [T(rand(-10:10)) for i in 1:num_vars]

         r1 = f(V1...)
         r2 = g(V1...)
         r3 = (f + g)(V1...)

         @test r3 == r1 + r2
      end

      for iter = 1:50
         f = rand(S, 0:5, 0:100, -100:100)
         g = rand(S, 0:5, 0:100, -100:100)

         V4 = [QQ(rand(-10:10)) for i in 1:num_vars]

         r1 = evaluate(f, V4)
         r2 = evaluate(g, V4)
         r3 = evaluate(f + g, V4)

         @test r3 == r1 + r2
      end
   end

   # Test ordering is correct, see issue #184
   for iter = 1:10
      ord = rand_ordering()
      R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"], ordering = ord)

      f = x*y^2*z^3

      @test evaluate(f, [2, 3, 5]) == 2*9*125
      @test evaluate(f, BigInt[2, 3, 5]) == 2*9*125
   end

   # Individual tests

   R, (x, y) = PolynomialRing(ZZ, ["x", "y"])

   f = 2x^2*y^2 + 3x + y + 1

   @test evaluate(f, BigInt[1, 2]) == ZZ(14)
   @test evaluate(f, [QQ(1), QQ(2)]) == 14//1
   @test evaluate(f, [1, 2]) == 14
   @test f(1, 2) == 14
   @test f(ZZ(1), ZZ(2)) == ZZ(14)
   @test f(QQ(1), QQ(2)) == 14//1

   @test evaluate(f, [x + y, 2y - x]) ==
               2*x^4 - 4*x^3*y - 6*x^2*y^2 + 8*x*y^3 + 2*x + 8*y^4 + 5*y + 1
   @test f(x + y, 2y - x) ==
               2*x^4 - 4*x^3*y - 6*x^2*y^2 + 8*x*y^3 + 2*x + 8*y^4 + 5*y + 1

   S, z = PolynomialRing(R, "z")

   @test evaluate(f, [z + 1, z - 1]) == 2*z^4 - 4*z^2 + 4*z + 5
   @test f(z + 1, z - 1) == 2*z^4 - 4*z^2 + 4*z + 5

   R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])

   f = x^2*y^2 + 2x*z + 3y*z + z + 1

   @test evaluate(f, [1, 3], [3, 4]) == 9*y^2 + 12*y + 29
   @test evaluate(f, [x, z], [3, 4]) == 9*y^2 + 12*y + 29

   @test evaluate(f, [1, 2], [x + z, x - z]) ==
                  x^4 - 2*x^2*z^2 + 5*x*z + z^4 - z^2 + z + 1
   @test evaluate(f, [x, y], [x + z, x - z]) ==
                  x^4 - 2*x^2*z^2 + 5*x*z + z^4 - z^2 + z + 1

   S, t = PolynomialRing(R, "t")
   T, (x1, y1, z1) = PolynomialRing(QQ, ["x", "y", "z"])

   @test evaluate(f, [2, 3], [t + 1, t - 1]) ==
                 (x^2 + 3)*t^2 + (2*x^2 + 2*x + 1)*t + (x^2 - 2*x - 3)
   @test evaluate(f, [y, z], [t + 1, t - 1]) ==
                 (x^2 + 3)*t^2 + (2*x^2 + 2*x + 1)*t + (x^2 - 2*x - 3)

   @test evaluate(f, [2, 4, 6], QQ) == 167//1
   @test evaluate(f, [1, 3], [2, 4], QQ) == 4*y1^2 + 12*y1 + 21
   @test evaluate(f, [x, z], [2, 4], QQ) == 4*y1^2 + 12*y1 + 21

   S = MatrixAlgebra(ZZ, 2)

   M1 = S([1 2; 3 4])
   M2 = S([2 3; 1 -1])
   M3 = S([-1 1; 1 1])

   @test evaluate(f, [M1, M2, M3]) == S([64 83; 124 149])
   @test f(M1, M2, M3) == S([64 83; 124 149])

   @test f(M1, ZZ(2), M3) == S([24 53; 69 110])
   @test f(M1, ZZ(2), 3) == S([56 52; 78 134])

   K = RealField
   R, (x, y) = PolynomialRing(K, ["x", "y"])
   @test evaluate(x + y, [K(1), K(1)]) isa BigFloat

   println("PASS")
end

function test_gen_mpoly_valuation()
   print("Generic.MPoly.valuation...")

   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = PolynomialRing(R, var_names, ordering = ord)

      for iter = 1:100
         f = S()
         g = S()
         while f == 0 || g == 0 || isconstant(g)
            f = rand(S, 0:5, 0:100, 0:0, -100:100)
            g = rand(S, 0:5, 0:100, 0:0, -100:100)
         end

         d1 = valuation(f, g)

         expn = rand(1:5)

         d2 = valuation(f*g^expn, g)

         @test d2 == d1 + expn

         d3, q3 = remove(f, g)

         @test d3 == d1
         @test f == q3*g^d3

         d4, q4 = remove(q3*g^expn, g)

         @test d4 == expn
         @test q4 == q3
      end
   end

   println("PASS")
end

function test_gen_mpoly_derivative()
   print("Generic.MPoly.derivative...")

   for num_vars=1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      R, vars = PolynomialRing(ZZ, var_names; ordering=ord)

      j = 1
      for v in vars
         for iter in 1:10
            f = rand(R, 5:10, 1:10, -100:100)
            g = rand(R, 5:10, 1:10, -100:100)

            @test derivative(f + g, v) == derivative(g, v) + derivative(f, v)
            @test derivative(g*f, v) == derivative(g, v)*f + derivative(f, v)*g
            @test derivative(f, j) == derivative(f, v)
         end
         @test derivative(one(R), v) == zero(R)
         @test derivative(zero(R), v) == zero(R)
         @test derivative(v, v) == one(R)
         j += 1
      end
   end

   println("PASS")
end

function test_gen_mpoly_change_base_ring()
   print("Generic.MPoly.change_base_ring...")
   
   F2 = ResidueRing(ZZ, 2)
   R, varsR = PolynomialRing(F2, ["x"])
   S, varsS = PolynomialRing(R, ["y"])
   f = x -> x^2
   change_base_ring(varsR[1] * varsS[1], f) == f(varsR[1]) * varsS[1]

   for num_vars=1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      R, vars = PolynomialRing(ZZ, var_names; ordering=ord)

      @test typeof(AbstractAlgebra.Generic.change_base_ring(R(1), a->ZZ(a))) == AbstractAlgebra.Generic.MPoly{typeof(ZZ(1))}
      @test typeof(AbstractAlgebra.Generic.change_base_ring(R(0), a->ZZ(a))) == AbstractAlgebra.Generic.MPoly{typeof(ZZ(0))}

      for iter in 1:10
         f = rand(R, 5:10, 1:10, -100:100)
         @test evaluate(change_base_ring(f, a -> R(a)), [one(R) for i=1:num_vars]) == sum(f.coeffs[i] for i=1:f.length)
         @test evaluate(change_base_ring(f, a -> R(a)), vars) == f
         @test parent(change_base_ring(f, a -> R(a))).ord == parent(f).ord
      end
   end

   println("PASS")
end

function test_gen_mpoly_vars()
   print("Generic.MPoly.vars...")

   for num_vars=1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      R, vars_R = PolynomialRing(ZZ, var_names; ordering=ord)

      for iter in 1:10
         f = rand(R, 5:10, 1:10, -100:100)
         @test length(vars(R(evaluate(change_base_ring(f, a -> R(a)), [one(R) for i=1:num_vars])))) == 0
         # @test issubset(vars(f), Set(gens(parent(f)))) # This fails, probably due to https://github.com/Nemocas/AbstractAlgebra.jl/issues/111
      end
   end

   println("PASS")
end

function test_gen_mpoly_combine_like_terms()
   print("Generic.MPoly.combine_like_terms...")

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      R, vars_R = PolynomialRing(ZZ, var_names; ordering=ord)

      for iter in 1:10
         f = R()
         while f == 0
            f = rand(R, 5:10, 1:10, -100:100)
         end

         lenf = length(f)
         f = setcoeff!(f, rand(1:lenf), 0)
         f = combine_like_terms!(f)

         @test length(f) == lenf - 1

         while length(f) < 2
            f = rand(R, 5:10, 1:10, -100:100)
         end

         lenf = length(f)
         nrand = rand(1:lenf - 1)
         v = exponent_vector(f, nrand)
         f = set_exponent_vector!(f, nrand + 1, v)
         terms_cancel = coeff(f, nrand) == -coeff(f, nrand + 1)
         f = combine_like_terms!(f)
         @test length(f) == lenf - 1 - terms_cancel
      end
   end

   println("PASS")
end

function test_gen_mpoly_exponents()
   print("Generic.MPoly.exponents...")

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      R, vars_R = PolynomialRing(ZZ, var_names; ordering=ord)

      for iter in 1:10
         f = R()
         while f == 0
            f = rand(R, 5:10, 1:10, -100:100)
         end

         nrand = rand(1:length(f))
         v = exponent_vector(f, nrand)
         c = coeff(f, v)

         @test c == coeff(f, nrand)
         @test v == [AbstractAlgebra.Generic.exponent(f, nrand, i) for i = 1:length(vars_R)]
      end

      for iter in 1:10
         num_vars = nvars(R)

         f = R()
         rand_len = rand(0:10)

         fit!(f, rand_len)

         for i = 1:rand_len
            expi = [rand(0:10) for j in 1:num_vars]
            ci = rand(ZZ, -10:10)

            f = set_exponent_vector!(f, i, expi)
            f = setcoeff!(f, i, ci)
         end

         f = sort_terms!(f)
         f = combine_like_terms!(f)

         for i = 1:length(f) - 1
            @test AbstractAlgebra.Generic.monomial_cmp(f.exps, i, f.exps, i + 1, R.N, R, UInt(0)) > 0
         end

         f = R()
         g = R()
         h = R()

         exp_arr = unique([[rand(0:10) for j in 1:num_vars] for k in 1:rand_len])
         for i = 1:length(exp_arr)
            expi = exp_arr[i]
            ci = rand(ZZ, -10:10)

            f = set_exponent_vector!(f, i, expi)
            f = setcoeff!(f, i, ci)

            if ci != 0
               g = setcoeff!(g, expi, ci)
               h = setcoeff!(h, expi, Int(ci))
            end
         end

         f = sort_terms!(f)
         f = combine_like_terms!(f)

         @test f == g
         @test f == h
      end
   end

   println("PASS")
end

function test_gen_mpoly_to_univariate()
   print("Generic.MPoly.to_univariate...")

   for num_vars=1:10
      ord = rand_ordering()
      
      var_names = ["x$j" for j in 1:num_vars]

      R, vars_R = PolynomialRing(ZZ, var_names; ordering=ord)
      x = rand(vars_R)

      R_univ, x_univ = PolynomialRing(ZZ, "x")

      @test zero(R_univ) == to_univariate(R_univ, zero(R))
      @test one(R_univ) == to_univariate(R_univ, one(R))

      for iter in 1:10
         f = zero(R)
         f_univ = zero(R_univ)
         coeffs = rand(Int, 100)
         for i in 1:100
            f = f + coeffs[i] * x^i
            f_univ = f_univ + coeffs[i] * x_univ^i
         end
         @test to_univariate(R_univ, f) == f_univ
      end
   end

   println("PASS")
end

function test_gen_mpoly_coefficients_of_univariate()
   print("Generic.MPoly.coefficients_of_univariate...")

   for num_vars=1:10
      ord = rand_ordering()
      var_names = ["x$j" for j in 1:num_vars]

      R, vars_R = PolynomialRing(ZZ, var_names; ordering=ord)

      @test length(AbstractAlgebra.Generic.coefficients_of_univariate(zero(R), true)) == 0
      @test length(AbstractAlgebra.Generic.coefficients_of_univariate(zero(R), false)) == 0
      @test AbstractAlgebra.Generic.coefficients_of_univariate(one(R), true) == [ one(base_ring(R)) ]
      x = rand(vars_R)
      @test AbstractAlgebra.Generic.coefficients_of_univariate(one(R), false) == [ one(base_ring(R)) ]
      @test AbstractAlgebra.Generic.coefficients_of_univariate(x, true) == [ zero(base_ring(R)), one(base_ring(R)) ]
      x = rand(vars_R)
      @test AbstractAlgebra.Generic.coefficients_of_univariate(x, false) == [ zero(base_ring(R)), one(base_ring(R)) ]

      for iter in 1:10
         f = zero(R)
         l = rand(1:1000)
         coeffs = rand(Int, l)
         for i in 1:l
            f = f + coeffs[i] * x^(i-1)
         end
         @test AbstractAlgebra.Generic.coefficients_of_univariate(f, true) == coeffs
         @test AbstractAlgebra.Generic.coefficients_of_univariate(f, false) == coeffs
      end
   end

   println("PASS")
end

function test_gen_mpoly_isless()
   print("Generic.MPoly.isless...")
   
   n_mpolys = 100
   maxval = 10
   maxdeg = 20
   
   # :deglex ordering
   R, (x,y,z) = AbstractAlgebra.Generic.PolynomialRing(AbstractAlgebra.Generic.ZZ, ["x","y","z"], ordering=:lex)
   # Monomials of degree 2
   @test isless(z^2, y*z) == true
   @test isless(y*z, y^2) == true
   @test isless(y^2, x*z) == true
   @test isless(x*z, x*y) == true
   @test isless(x*y, x^2) == true

   for n_vars = 1:maxdeg
      A = unique(sortslices(reshape(map(Int,map(round, rand(n_vars * n_mpolys) * maxval)), (n_mpolys, n_vars)), dims=1),dims=1)
      var_names = ["x$j" for j in 1:n_vars]
      R, varsR = AbstractAlgebra.Generic.PolynomialRing(AbstractAlgebra.Generic.ZZ, var_names, ordering=:lex)
      for i in 1:size(A)[1]-1
         f = R([base_ring(R)(1)], [A[i,:]])
         g = R([base_ring(R)(1)], [A[i+1,:]])
         @test isless(f,g)
      end
   end

   # :deglex ordering
   R, (x,y,z) = AbstractAlgebra.Generic.PolynomialRing(AbstractAlgebra.Generic.ZZ, ["x","y","z"], ordering=:deglex)
   
   @test isless(z^2, y*z) == true
   @test isless(y*z, x*z) == true
   @test isless(y^2, x*z) == true
   @test isless(y^2, x*y) == true
   @test isless(x*y, x^2) == true
   
   for n_vars=1:maxdeg
      A = reshape(map(Int,map(round, rand(n_vars * n_mpolys) * maxval)), (n_mpolys, n_vars))
      var_names = ["x$j" for j in 1:n_vars]
      R, varsR = AbstractAlgebra.Generic.PolynomialRing(AbstractAlgebra.Generic.ZZ, var_names, ordering=:deglex)
      
      for i in 1:size(A)[1]-1
         f = R([base_ring(R)(1)], [A[i,:]])
         g = R([base_ring(R)(1)], [A[i+1,:]])
         if total_degree(f) < total_degree(g)
            @test isless(f,g)
         elseif total_degree(g) < total_degree(f)
            @test isless(g,f)
         else
            for j = 1:n_vars
               if A[i, j] < A[i+1, j]
                  @test isless(f,g)
                  break
               elseif A[i,j] > A[i+1,j]
                  @test isless(g,f)
                  break
               end
            end
         end
      end
   end
   
   # :degrevlex ordering
   R, (x,y,z) = AbstractAlgebra.Generic.PolynomialRing(AbstractAlgebra.Generic.ZZ, ["x","y","z"], ordering=:degrevlex)
   # Monomials of degree 2
   @test isless(z^2, y*z) == true
   @test isless(y*z, x*z) == true
   @test isless(x*z, y^2) == true
   @test isless(y^2, x*y) == true
   @test isless(x*y, x^2) == true
   for n_vars = 1:maxdeg
      A = reshape(map(Int,map(round, rand(n_vars * n_mpolys) * maxval)), (n_mpolys, n_vars))
      var_names = ["x$j" for j in 1:n_vars]
      R, varsR = AbstractAlgebra.Generic.PolynomialRing(AbstractAlgebra.Generic.ZZ, var_names, ordering=:degrevlex)
      for i in 1:size(A)[1]-1
         f = R([base_ring(R)(1)], [A[i,:]])
         g = R([base_ring(R)(1)], [A[i+1,:]])
         if total_degree(f) < total_degree(g)
            @test isless(f,g)
         elseif total_degree(g) < total_degree(f)
            @test isless(g,f)
         else
            for j = n_vars:-1:1
               if A[i, j] > A[i+1, j]
                  @test isless(f,g)
                  break
               elseif A[i, j] == A[i+1, j]
                  continue
               elseif A[i,j] < A[i+1,j]
                  @test isless(g,f)
                  break
               end
            end
         end
      end
   end
   println("PASS")
end

function test_gen_mpoly_lt()
   for num_vars=1:10
      ord = rand_ordering()
      var_names = ["x$j" for j in 1:num_vars]

      R, vars_R = PolynomialRing(ZZ, var_names; ordering=ord)
      
      f = rand(R, 5:10, 1:10, -100:100)
      g = rand(R, 5:10, 1:10, -100:100)

      @test lt(f*g) == lt(f)*lt(g)
      @test lt(one(R)) == one(R)
      @test lt(zero(R)) == zero(R)

      for v in vars_R
         @test lt(v) == v
      end

      @test parent(lt(f)) == parent(f)
   end
end

function test_gen_mpoly_lc()
   print("Generic.MPoly.lc...")

   for num_vars = 1:4
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()
      S, varlist = PolynomialRing(ZZ, var_names, ordering = ord)

      for iter = 1:10
         f = rand(S, 0:4, 0:5, -10:10)
         g = rand(S, 0:4, 0:5, -10:10)

         @test lc(f*g) == lc(f)*lc(g)
         @test lc(one(S)) == one(base_ring(S))

         for v in varlist
            @test lc(v) == one(base_ring(S))
         end

         @test parent(lc(f)) == base_ring(f)
      end
   end

   println("PASS")
end

function test_gen_mpoly_lm()
   print("Generic.MPoly.lm...")

   for num_vars = 1:4
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()
      S, varlist = PolynomialRing(ZZ, var_names, ordering = ord)

      for iter = 1:10
         f = zero(S)
         while iszero(f)
            f = rand(S, 0:4, 0:5, -10:10)
         end
         g = zero(S)
         while iszero(g)
            g = rand(S, 0:4, 0:5, -10:10)
         end

         @test lm(f*g) == lm(f)*lm(g)
         @test lm(one(S)) == one(S)
         @test lm(zero(S)) == zero(S)

         for v in varlist
            @test lm(v) == v
         end

         @test parent(lm(f)) == parent(f)
      end
   end

   println("PASS")
end

function test_gen_mpoly_ishomogeneous()
   print("Generic.MPoly.ishomogenous...")

   for num_vars = 1:4
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()
      S, varlist = PolynomialRing(ZZ, var_names, ordering = ord)

      for iter = 1:10
         @test ishomogeneous(zero(S))
         @test ishomogeneous(one(S))
         for v in varlist
            @test ishomogeneous(v)
            @test !ishomogeneous(v+one(S))
         end
      end
   end

   println("PASS")
end


function test_gen_mpoly()
   test_gen_mpoly_constructors()
   test_gen_mpoly_geobuckets()
   test_gen_mpoly_manipulation()
   test_gen_mpoly_multivariate_coeff()
   test_gen_mpoly_unary_ops()
   test_gen_mpoly_binary_ops()
   test_gen_mpoly_adhoc_binary()
   test_gen_mpoly_adhoc_comparison()
   test_gen_mpoly_powering()
   test_gen_mpoly_divides()
   test_gen_mpoly_euclidean_division()
   test_gen_mpoly_ideal_reduction()
   test_gen_mpoly_deflation()
   test_gen_mpoly_gcd()
   test_gen_mpoly_evaluation()
   test_gen_mpoly_valuation()
   test_gen_mpoly_derivative()
   test_gen_mpoly_change_base_ring()
   test_gen_mpoly_total_degree()
   test_gen_mpoly_vars()
   test_gen_mpoly_combine_like_terms()
   test_gen_mpoly_exponents()
   test_gen_mpoly_to_univariate()
   test_gen_mpoly_coefficients_of_univariate()
   test_gen_mpoly_isless()
   test_gen_mpoly_lc()
   test_gen_mpoly_lm()
   test_gen_mpoly_lt()
   test_gen_mpoly_lcm()
   test_gen_mpoly_ishomogeneous()

   println("")
end
