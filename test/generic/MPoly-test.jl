@testset "Generic.MPoly.constructors" begin
   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

      @test polynomial_ring(R, var_names, internal_ordering = ord, cached = true)[1] === polynomial_ring(R, var_names, internal_ordering = ord, cached = true)[1]
      @test polynomial_ring(R, var_names, internal_ordering = ord, cached = false)[1] !== polynomial_ring(R, var_names, internal_ordering = ord, cached = true)[1]
      @test polynomial_ring(R, num_vars, "x", internal_ordering = ord, cached = true)[1] === polynomial_ring(R, var_names, internal_ordering = ord, cached = true)[1]
      @test polynomial_ring(R, num_vars, internal_ordering = ord, cached = true)[1] === polynomial_ring(R, var_names, internal_ordering = ord, cached = true)[1]

      @test elem_type(S) == Generic.MPoly{elem_type(R)}
      @test elem_type(Generic.MPolyRing{elem_type(R)}) == Generic.MPoly{elem_type(R)}
      @test parent_type(Generic.MPoly{elem_type(R)}) == Generic.MPolyRing{elem_type(R)}
      @test base_ring(S) === R
      @test coefficient_ring(S) === R
      @test coefficient_ring_type(S) === typeof(R)

      @test S isa Generic.MPolyRing

      @test isa(symbols(S), Vector{Symbol})

      for j = 1:num_vars
         @test isa(varlist[j], MPolyRingElem)
         @test isa(gens(S)[j], MPolyRingElem)
      end

      f =  rand(S, 0:5, 0:100, 0:0, -100:100)

      @test isa(f, MPolyRingElem)

      @test isa(S(2), MPolyRingElem)

      @test isa(S(R(2)), MPolyRingElem)

      @test isa(S(f), MPolyRingElem)

      V = [R(rand(-100:100)) for i in 1:5]

      W0 = UInt[rand(0:100) for i in 1:5*num_vars]
      W = reshape(W0, num_vars, 5)

      @test isa(S(V, W), MPolyRingElem)

      W1 = [[rand(0:100) for i in 1:num_vars] for j in 1:5]

      f1 = S(V, W1)

      @test isa(f1, MPolyRingElem)

      Va = [rand(-100:100) for i in 1:5]
      f1a = S(Va, W1)

      @test isa(f1a, MPolyRingElem)

      Vb = [ZZ(rand(-100:100)) for i in 1:5]
      f1b = S(Vb, W1)

      @test isa(f1b, MPolyRingElem)

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

      @test_throws ErrorException push_term!(C, one(base_ring(S)), zeros(Int, num_vars + 1))

      for i = 1:5
         push_term!(C, V[i], W1[i])
      end
      f3 = finish(C)

      @test f1 == f3

      _, varlist = polynomial_ring(QQ, var_names)
      y = varlist[1]
      @test x in [x, y]
      @test x in [y, x]
      @test !(x in [y])
      @test x in keys(Dict(x => 1))
      @test !(y in keys(Dict(x => 1)))
   end

   R1, (x, y) = QQ["x", "y"]

   B = MPolyBuildCtx(R1)

   push_term!(B, QQ(6), [6, 6])
   p = finish(B)

   @test p == 6x^6*y^6

   push_term!(B, QQ(5), [5, 5])

   q = finish(B)

   @test q == 5x^5*y^5

   # test "getindex" syntax
   S, (y, z) = R["y", "z"]
   @test S isa Generic.MPolyRing{Generic.Poly{BigInt}}
   @test y isa Generic.MPoly{Generic.Poly{BigInt}}
   @test z isa Generic.MPoly{Generic.Poly{BigInt}}

   ZZxyz_ = @inferred polynomial_ring(ZZ, 'x':'z')
   ZZxyz2, xyz2 = @inferred polynomial_ring(ZZ, VarName[:x, 'y', GenericString("z")])
   ZZxyz4_ = @inferred polynomial_ring(ZZ, Union{String,Char,Symbol}["x", 'y', :z])
   ZZxyz5_ = @inferred ZZ["x", 'y', :z]
   ZZxyz6 = @polynomial_ring(ZZ, [:x, :y, :z])

   @test ZZxyz_[1] isa Generic.MPolyRing
   @test ZZxyz_ == (ZZxyz2, xyz2)
   @test ZZxyz_ == ZZxyz4_
   @test ZZxyz_ == ZZxyz5_
   @test ZZxyz_ == (ZZxyz6, [x, y, z])

   ZZxxx0_ = @inferred polynomial_ring(ZZ, :x=>Base.OneTo(3))
   ZZxxx_ = @inferred polynomial_ring(ZZ, :x=>1:3)

   @test ZZxxx_[1] isa Generic.MPolyRing
   @test ZZxxx_ == ZZxxx0_

   QQxxx_ = @inferred polynomial_ring(QQ, "x#" => 1:3)
   QQxxx00_ = @inferred polynomial_ring(QQ, 3)
   QQxxx0 = @polynomial_ring(QQ, 3)

   @test QQxxx_[1] isa Generic.MPolyRing
   @test QQxxx_ == QQxxx00_
   @test QQxxx_ == (QQxxx0, [x1, x2, x3])

   QQxxx1 = @polynomial_ring(QQ, 3, "x")
   @test QQxxx_ == (QQxxx1, [x1, x2, x3])

   QQxxx2 = @polynomial_ring(QQ, "x#" => 1:3)
   @test QQxxx_ == (QQxxx2, [x1, x2, x3])

   QQxxx3 = @polynomial_ring(QQ, :x=>1:3)
   @test QQxxx_ == (QQxxx3, [x1, x2, x3])

   # test support of kwargs (issue #1631)
   ordering = :lex
   QQxxx4 = @polynomial_ring(QQ, "x#" => 1:3; internal_ordering=ordering)
   @test QQxxx_ == (QQxxx4, [x1, x2, x3])
   QQxxx5 = @polynomial_ring(QQ, "x#" => 1:3, internal_ordering=ordering)
   @test QQxxx_ == (QQxxx5, [x1, x2, x3])

   QQxxx_deglex_ = @inferred polynomial_ring(QQ, "x#" => 1:3; internal_ordering=:deglex)
   @test internal_ordering(QQxxx_deglex_[1]) == :deglex

   QQxxx_deglex2 = @polynomial_ring(QQ, "x#" => 1:3; internal_ordering=:deglex)
   @test internal_ordering(QQxxx_deglex2) == :deglex
   @test QQxxx_deglex_ == (QQxxx_deglex2, [x1, x2, x3])

   ZZxy_ = @inferred polynomial_ring(ZZ, :x => (1:2, 1:2), :y => 0:3)

   ZZxy2_ = @inferred polynomial_ring(ZZ, :x => ["1, 1" "1, 2"; "2, 1" "2, 2"], :y => (0:3,))

   @test ZZxy_[1] isa Generic.MPolyRing
   @test ZZxy_ == ZZxy2_

   QQxy_ = @inferred polynomial_ring(QQ, "x#" => (1:2, 1:2), Symbol.(:y, 0:3))
   QQxy2 = @polynomial_ring(QQ, "x#" => (1:2, 1:2), Symbol.(:y, 0:3))

   @test QQxy_[1] isa Generic.MPolyRing
   @test QQxy_ == (QQxy2, [x11 x12; x21 x22], [y0, y1, y2, y3])

   QQxy3 = @polynomial_ring(QQ, :x => (1:2, 1:2), :y => 0:3)
   @test QQxy_ == (QQxy3, [x11 x12; x21 x22], [y0, y1, y2, y3])

   # Errors
   @test_throws ArgumentError polynomial_ring(QQ, "x###" => (1:2, 3:4))
   @test_logs (:warn, """The variable name "x-1" sadly is no Julia identifier. You can still access it as `var"x-1"`."""
      ) @macroexpand @polynomial_ring(QQ, :x => -1:1)
   @test_logs (:error, "Inconveniently, you may only use literals and variables from the global scope of the current module (`Main`) when using variable name constructor macros"
      ) @test_throws UndefVarError let local_name = 3
         @macroexpand @polynomial_ring(QQ, :x => 1:local_name)
      end
end

# these variables need to be in global scope
varname_demo = :a => 1:3
varname_demo0 = "a#" => 1:3 # for use in function constructor
varname_three = 3
@testset "Generic.MPoly.constructors.macroscoping" begin
   # test scoping of variables in macro constructor (issue #1630)
   S = ZZ
   Sx, x = polynomial_ring(S, :x)
   Sx2 = @polynomial_ring(S, :x)
   @test Sx == Sx2

   Svars, _ = polynomial_ring(S, varname_demo0)
   Svars1 = @polynomial_ring(S, varname_demo)
   Svars2 = @polynomial_ring(S, varname_demo0)
   Svars3 = @polynomial_ring(S, :a => 1:varname_three)
   @test Svars == Svars1
   @test Svars == Svars2
   @test Svars == Svars3

   K = QQ
   Kxy, _ = polynomial_ring(K, [:x, :y])
   Kxy2 = @polynomial_ring(K, [:x, :y])
   @test Kxy == Kxy2

   # some more tests for different ways to specify kwargs (issue #1631)
   # combined with scoping of keyword and all other args
   R = Kxy
   cached = true
   ordering = :deglex
   kwargs = (; internal_ordering = ordering, cached)
   Rvars, _ = polynomial_ring(R, varname_demo0; internal_ordering=ordering, cached)
   Rvars1 = @polynomial_ring(R, varname_demo; internal_ordering=ordering, cached)
   Rvars2 = @polynomial_ring(R, varname_demo, internal_ordering=ordering; cached)
   Rvars3 = @polynomial_ring(R, varname_demo, internal_ordering=ordering, cached=cached)
   Rvars4 = @polynomial_ring(R, varname_demo, cached=true; kwargs...)
   Rvars5 = @polynomial_ring(R, varname_demo; kwargs...)
   Rvars6 = @polynomial_ring(R, :a => 1:varname_three; kwargs...)
   @test Rvars == Rvars1
   @test Rvars == Rvars2
   @test Rvars == Rvars3
   @test Rvars == Rvars4
   @test Rvars == Rvars5
   @test Rvars == Rvars6
end

@testset "Generic.MPoly.printing" begin
   S, (x, y) = ZZ["x", "y"]

   @test string(zero(S)) == "0"
   @test string(one(S)) == "1"
   @test string(x) == "x"
   @test string(y) == "y"
   @test string(x^2 - y) == "x^2 - y"

   S, (x, y) = RealField["x", "y"]

   @test string(x) == "x"
   @test string(y) == "y"
end

@testset "Generic.MPoly.geobuckets" begin
   R, x = ZZ["y"]
   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)
      g = gens(S)

      C = Generic.geobucket(S)
      for j in 1:num_vars
        push!(C, sum(g[j]^i for i in 1:64))
      end
      @test finish(C) == sum(sum(g[j]^i for i in 1:64) for j in 1:num_vars)
   end
end

@testset "Generic.MPoly.rand" begin
   R, x = ZZ["y"]
   num_vars = 5
   var_names = ["x$j" for j in 1:num_vars]

   ord = rand_ordering()
   @test ord in [:lex, :deglex, :degrevlex]
   ord = rand_ordering(rng)
   @test ord in [:lex, :deglex, :degrevlex]

   S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

   test_rand(S, 0:5, 0:100, 0:0, -100:100)
end

@testset "Generic.MPoly.manipulation" begin
   R, x = ZZ["y"]

   @test characteristic(R) == 0

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)
      g = gens(S)

      @test !is_gen(S(1))

      for i = 1:num_vars
         @test is_gen(varlist[i])
         @test is_gen(g[i])
         @test !is_gen(g[i] + 1)
         @test gen(S, i) == g[i]
         @test gen(S, i) == S[i]
         @test var_index(gen(S, i)) == i
      end

      @test_throws ArgumentError gen(S, 0)
      @test_throws ArgumentError gen(S, num_vars + 1)

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
      @test f == copy(f)

      @test hash(f) == hash(deepcopy(f))

      if length(f) > 0
         @test isa(coeff(f, rand(1:length(f))), elem_type(R))
      end

      @test length(f) == length(coefficients(f))
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
      deg = is_degree(internal_ordering(S))
      rev = is_reverse(internal_ordering(S))

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

      @test ord == internal_ordering(S)

      @test nvars(parent(f)) == num_vars
      @test ngens(parent(f)) == num_vars

      @test isone(one(S))

      @test iszero(zero(S))

      @test is_unit(S(1))
      @test !is_unit(S(0))
      @test !is_unit(gen(S, 1))

      @test is_zero_divisor(S())
      @test !is_zero_divisor(gen(S, 1))

      @test is_constant(S(rand(-100:100)))
      @test is_constant(S(zero(S)))

      g = rand(S, 1:1, 0:100, 0:0, 1:100)
      h = rand(S, 1:1, 0:100, 0:0, 1:1)

      @test is_term(g)
      @test is_monomial(h)
   end

   for num_vars = 1:4
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()
      S, varlist = polynomial_ring(ZZ, var_names, internal_ordering = ord)

      for iter = 1:10
         @test is_homogeneous(zero(S))
         @test is_homogeneous(one(S))
         for v in varlist
            @test is_homogeneous(v)
            @test !is_homogeneous(v + one(S))
         end
      end
   end

   R, (x, ) = polynomial_ring(ZZ, ["x"])

   @test is_univariate(R)
   @test is_univariate(x)
   @test is_univariate(R())
   @test is_univariate(R(1))
   @test is_univariate(x^3 + 3x)

   R, (x, y) = polynomial_ring(ZZ, ["x", "y"])

   @test !is_univariate(R)
   @test is_univariate(x)
   @test is_univariate(y)
   @test is_univariate(R())
   @test is_univariate(R(1))
   @test is_univariate(x^3 + 3*x)
   @test is_univariate(2*y^4 + 3*y + 5)
   @test !is_univariate(x + y)
   @test !is_univariate(x*y)
   @test !is_univariate(x^3 + 3x + y + 1)
   @test !is_univariate(x^3 + 3x + y)
   @test !is_univariate(y^4 + 3x + 1)

   @test is_univariate_with_data(y) == (true, 2)
   @test is_univariate_with_data(R()) == (true, 0)
   @test is_univariate_with_data(x + y) == (false, 0)
end


@testset "Generic.MPoly.multivariate_coeff" begin
   R = ZZ

   for iter = 1:5
      ord = rand_ordering()

      S, (x, y, z) = polynomial_ring(R, ["x", "y", "z"]; internal_ordering=ord)

      f = -8*x^5*y^3*z^5+9*x^5*y^2*z^3-8*x^4*y^5*z^4-10*x^4*y^3*z^2+8*x^3*y^2*z-10*x*y^3*
z^4-4*x*y-10*x*z^2+8*y^2*z^5-9*y^2*z^3

      @test coeff(f, [1], [1]) == -10*y^3*z^4-4*y-10*z^2
      @test coeff(f, [2, 3], [3, 2]) == -10*x^4
      @test coeff(f, [1, 3], [4, 5]) == 0

      @test coeff(f, [x], [1]) == -10*y^3*z^4-4*y-10*z^2
      @test coeff(f, [y, z], [3, 2]) == -10*x^4
      @test coeff(f, [x, z], [4, 5]) == 0
   end
end

@testset "Generic.MPoly.leading_term" begin
   for num_vars=1:10
      ord = rand_ordering()
      var_names = ["x$j" for j in 1:num_vars]

      R, vars_R = polynomial_ring(ZZ, var_names; internal_ordering=ord)

      f = rand(R, 5:10, 1:10, -100:100)
      g = rand(R, 5:10, 1:10, -100:100)

      if !iszero(f) && !iszero(g)
         @test leading_term(f*g) == leading_term(f)*leading_term(g)
      else
         @test_throws ArgumentError leading_term(f)*leading_term(g)
      end
      @test leading_term(one(R)) == one(R)
      @test_throws ArgumentError leading_term(zero(R))

      for v in vars_R
         @test leading_term(v) == v
      end

      @test parent(leading_term(f)) == parent(f)
   end

   for num_vars = 1:4
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()
      S, varlist = polynomial_ring(ZZ, var_names, internal_ordering = ord)

      for iter = 1:10
         f = rand(S, 0:4, 0:5, -10:10)
         g = rand(S, 0:4, 0:5, -10:10)

         @test leading_coefficient(f*g) ==
               leading_coefficient(f)*leading_coefficient(g)
         @test leading_coefficient(one(S)) == one(base_ring(S))

         for v in varlist
            @test leading_coefficient(v) == one(base_ring(S))
         end

         @test parent(leading_coefficient(f)) == base_ring(f)
      end
   end

   for num_vars = 1:4
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()
      S, varlist = polynomial_ring(ZZ, var_names, internal_ordering = ord)

      for iter = 1:10
         f = zero(S)
         while iszero(f)
            f = rand(S, 0:4, 0:5, -10:10)
         end
         g = zero(S)
         while iszero(g)
            g = rand(S, 0:4, 0:5, -10:10)
         end

         @test leading_monomial(f*g) == leading_monomial(f)*leading_monomial(g)
         @test leading_exponent_vector(f*g) == leading_exponent_vector(f) +
                                               leading_exponent_vector(g)
         @test leading_monomial(one(S)) == one(S)
         @test_throws ArgumentError leading_monomial(zero(S))
         @test_throws ArgumentError leading_exponent_vector(zero(S))

         for v in varlist
            @test leading_monomial(v) == v
         end

         @test parent(leading_monomial(f)) == parent(f)
      end
   end

   R, (x, y) = polynomial_ring(ZZ, ["x", "y"])

   @test constant_coefficient(R()) == 0
   @test constant_coefficient(2x + 1) == 1
   @test constant_coefficient(2x) == 0
   @test constant_coefficient(2x^2 + 3y^3 + 4) == 4

   @test trailing_coefficient(x^2*y + 7x*y + 3x + 2y + 5) == 5
   @test trailing_coefficient(x^2*y + 7x*y + 3x + 2y) == 2
   @test trailing_coefficient(R(2)) == 2
   @test trailing_coefficient(R()) == 0

   @test tail(2x^2 + 2x*y + 3) == 2x*y + 3
   @test tail(R(1)) == 0
   @test tail(R()) == 0
end

@testset "Generic.MPoly.total_degree" begin
   max = 50
   for nvars = 1:10
      var_names = ["x$j" for j in 1:nvars]
      for nterms = 1:10
         exps = Matrix{Int}(round.(rand(nvars, nterms) .* max))
         degrees = []
         for nord = 1:20
            ord = rand_ordering()
            S, varlist = polynomial_ring(ZZ, var_names, internal_ordering = ord)
            p = zero(S)
            for j = 1:nterms
               p += prod(varlist[i]^exps[i,j] for i = 1:nvars)
            end
            push!(degrees, total_degree(p))
         end
         @test length(Set(degrees)) == 1
      end
   end
end

@testset "Generic.MPoly.unary_ops" begin
   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

      for iter = 1:10
         f = rand(S, 0:5, 0:100, 0:0, -100:100)

         @test f == -(-f)
      end
   end
end

@testset "Generic.MPoly.binary_ops" begin
   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

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
end

@testset "Generic.MPoly.adhoc_binary" begin
   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

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
end

@testset "Generic.MPoly.adhoc_comparison" begin
   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

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
end

@testset "Generic.MPoly.powering" begin
   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

      for iter = 1:10
         f = rand(S, 0:5, 0:100, 0:0, -100:100)

         expn = rand(0:10)

         r = S(1)
         for i = 1:expn
            r *= f
         end

         @test f^expn == r
      end

      @test_throws DomainError rand(varlist)^-1
      @test_throws DomainError rand(varlist)^-rand(2:100)
   end

   # Over field of nonzero characteristic
   R, (x, y, z, t) = polynomial_ring(GF(2), ["x", "y", "z", "t"])
   f = 1 + x + y + z + t

   @test zero(R)^0 == one(R)

   for i = 1:5
      @test f^i == f*f^(i - 1)
   end

   # Over non-exact field
   F = RealField
   R, u = polynomial_ring(F, 3)
   x = F(0.77)*u[1] + F(0.77)*one(F)*u[2] + F(0.63)*u[3]
   @test x^5 == x*x*x*x*x
end

@testset "Generic.MPoly.divides" begin
   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

      for iter = 1:10
         f = S(0)
         f = rand(S, 0:5, 0:100, 0:0, -100:100)
         g = rand(S, 0:5, 0:100, 0:0, -100:100)

         p = f*g

         flag, q = divides(p, f)
         flag2, q2 = divides(f, p)

         @test flag == true
         @test is_divisible_by(p, f)
         @test q * f == p

         q1 = divexact(p, f)

         @test q1 * f == p

         if !iszero(p)
           @test q1 == g
         end
      end
   end
end

@testset "Generic.MPoly.square_root" begin
   for R in [ZZ, QQ]
      for num_vars = 1:10
         var_names = ["x$j" for j in 1:num_vars]
         ord = rand_ordering()

         S, varlist = polynomial_ring(ZZ, var_names, internal_ordering = ord)

         for iter = 1:10
            f = rand(S, 0:5, 0:100, -100:100)

            p = f^2

            @test is_square(p)

            q = sqrt(f^2)

            @test q^2 == f^2

            q = sqrt(f^2; check=false)

            @test q^2 == f^2

            if f != 0
               x = varlist[rand(1:num_vars)]
               @test_throws ErrorException sqrt(f^2*(x^2 - x))
            end

            f1, s1 = is_square_with_sqrt(f)

            @test !f1 || s1^2 == f
         end
      end
   end

   # Field of characteristic p
   for p in [2, 7, 13, 65537, ZZ(2), ZZ(7), ZZ(37), ZZ(65537)]
      R, = residue_field(ZZ, p)
      for num_vars = 1:10
         var_names = ["x$j" for j in 1:num_vars]
         ord = rand_ordering()

         S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

         for iter = 1:10
            f = rand(S, 0:5, 0:100, 0:Int(p))

            s = f^2

            @test is_square(s)

            q = sqrt(f^2)

            @test q^2 == f^2

            q = sqrt(f^2)

            @test q^2 == f^2

            if f != 0
               x = varlist[rand(1:num_vars)]
               @test_throws ErrorException sqrt(f^2*(x^2 - x))
            end

            f1, s1 = is_square_with_sqrt(f)

            @test !f1 || s1^2 == f
         end
      end
   end
end

@testset "Generic.MPoly.euclidean_division" begin
   R, x = QQ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

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

      S, varlist = polynomial_ring(QQ, var_names, internal_ordering = ord)
      v = varlist[1+Int(round(rand() * (num_vars-1)))]
      @test divrem(v, 2*v) == (1//2, 0)
   end
end

@testset "Generic.MPoly.ideal_reduction" begin
   R, x = QQ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

      @test divrem(varlist[1], elem_type(S)[]) == (elem_type(S)[], varlist[1])

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

         V = Vector{elem_type(S)}(undef, num)

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
end

@testset "Generic.MPoly.deflation" begin
   for num_vars = 1:4
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()
      S, varlist = polynomial_ring(ZZ, var_names, internal_ordering = ord)

      for iter = 1:10
         p = rand(S, 0:4, 0:5, -10:10)
         shift = [rand(0:10) for i in 1:num_vars]
         defl = [rand(1:10) for i in 1:num_vars]
         f = inflate(p, shift, defl)

         s, d = deflation(f)
         g = deflate(f, s, d)
         h = inflate(g, s, d)

         @test h == f

         @test deflate(inflate(f, d), d) == f

         g = inflate(f, defl)
         h, defl = deflate(g)
         @test g == inflate(h, defl)

         vars = unique([rand(varlist) for _ in varlist])
         shift = [rand(0:10) for _ in vars]
         defl = [rand(1:10) for _ in vars]

         @test p == deflate(inflate(p, vars, shift, defl), vars, shift, defl)

         x = rand(varlist)
         f = p + evaluate(p, [x], [1])  # otherwise f0 is usually zero
         f0 = coeff(f, [x], [0])
         f1 = deflate(f - f0, [x], [1], [1])

         @test f == f0 + x * f1
         @test f == f0 + inflate(f1, [x], [1], [1])
      end
   end
end

@testset "Generic.MPoly.gcd" begin
   for num_vars = 1:4
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()
      S, varlist = polynomial_ring(ZZ, var_names, internal_ordering = ord)

      for iter = 1:10
         f = rand(S, 0:4, 0:5, -10:10)
         g = rand(S, 0:4, 0:5, -10:10)
         h = rand(S, 0:4, 0:5, -10:10)

         g1 = gcd(f, g)
         g2 = gcd(f*h, g*h)

         @test g2 == g1*h || g2 == -g1*h
      end
   end
end

@testset "Generic.MPoly.lcm" begin
   for num_vars = 1:4
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()
      S, varlist = polynomial_ring(ZZ, var_names, internal_ordering = ord)

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
end

@testset "Generic.MPoly.evaluation" begin
   R, x = ZZ["x"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

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

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

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

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

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
         Vval = Vector{elem_type(R)}(undef, 0) # value substituted for variables in V
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
         Wval = Vector{elem_type(R)}(undef, 0) # values for those variables

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
   T = matrix_ring(R, 2)

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

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
      R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"], internal_ordering = ord)

      f = x*y^2*z^3

      @test evaluate(f, [2, 3, 5]) == 2*9*125
      @test evaluate(f, BigInt[2, 3, 5]) == 2*9*125
   end

   # Individual tests

   R, (x, y) = polynomial_ring(ZZ, ["x", "y"])

   f = 2x^2*y^2 + 3x + y + 1

   @test evaluate(f, [0*x, 0*y]) == 1
   @test evaluate(f, Any[0*x, 0*y]) == 1
   @test_throws ArgumentError evaluate(f, [])

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

   @test evaluate(f, [x, 0]) == 3x + 1 # see https://github.com/oscar-system/Oscar.jl/issues/2331

   @test f(x=1, y=2) == 14
   @test f(y=2, x=1) == 14
   @test f(x=x+y, y=2y-x) ==
               2*x^4 - 4*x^3*y - 6*x^2*y^2 + 8*x*y^3 + 2*x + 8*y^4 + 5*y + 1
   @test f(y=2y-x, x=x+y) ==
               2*x^4 - 4*x^3*y - 6*x^2*y^2 + 8*x*y^3 + 2*x + 8*y^4 + 5*y + 1

   S, z = polynomial_ring(R, "z")

   @test evaluate(f, [z + 1, z - 1]) == 2*z^4 - 4*z^2 + 4*z + 5
   @test f(z + 1, z - 1) == 2*z^4 - 4*z^2 + 4*z + 5

   R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"])

   f = x^2*y^2 + 2x*z + 3y*z + z + 1

   @test evaluate(f, [1, 3], [3, 4]) == 9*y^2 + 12*y + 29
   @test evaluate(f, [x, z], [3, 4]) == 9*y^2 + 12*y + 29

   @test evaluate(f, [1, 2], [x + z, x - z]) ==
                  x^4 - 2*x^2*z^2 + 5*x*z + z^4 - z^2 + z + 1
   @test evaluate(f, [x, y], [x + z, x - z]) ==
                  x^4 - 2*x^2*z^2 + 5*x*z + z^4 - z^2 + z + 1

   S, t = polynomial_ring(R, "t")
   T, (x1, y1, z1) = polynomial_ring(QQ, ["x", "y", "z"])
   f1 = x1^2*y1^2 + 2x1*z1 + 3y1*z1 + z1 + 1

   @test evaluate(f, [2, 3], [t + 1, t - 1]) ==
                 (x^2 + 3)*t^2 + (2*x^2 + 2*x + 1)*t + (x^2 - 2*x - 3)
   @test evaluate(f, [y, z], [t + 1, t - 1]) ==
                 (x^2 + 3)*t^2 + (2*x^2 + 2*x + 1)*t + (x^2 - 2*x - 3)

   @test evaluate(change_base_ring(QQ, f1), [2, 4, 6]) == 167//1
   @test evaluate(f1, [1, 3], [2, 4]) == 4*y1^2 + 12*y1 + 21
   @test evaluate(f1, [x1, z1], [2, 4]) == 4*y1^2 + 12*y1 + 21

   S = matrix_ring(ZZ, 2)

   M1 = S([1 2; 3 4])
   M2 = S([2 3; 1 -1])
   M3 = S([-1 1; 1 1])

   @test evaluate(f, [M1, M2, M3]) == S([64 83; 124 149])
   @test f(M1, M2, M3) == S([64 83; 124 149])

   @test f(M1, ZZ(2), M3) == S([24 53; 69 110])
   @test f(M1, ZZ(2), 3) == S([56 52; 78 134])

   K = RealField
   R, (x, y) = polynomial_ring(K, ["x", "y"])
   @test evaluate(x + y, [K(1), K(1)]) isa BigFloat

   # Issue oscar-system/Oscar.jl#4762
   F,t = rational_function_field(QQ, :t)
   P,(x,y) = polynomial_ring(F, [:x, :y])
   @test x(t,y) == t
   @test x == gen(P, 1) # evaluation used to modify the polynomial

   # Issue #1219
   Qx, (x, y) = QQ["x", "y"];
   @test typeof(zero(Qx)(x, y)) == typeof(one(Qx)(x, y)) == typeof((x+y)(x, y))
end

@testset "Generic.MPoly.valuation" begin
   R, x = ZZ["y"]

   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      S, varlist = polynomial_ring(R, var_names, internal_ordering = ord)

      for iter = 1:100
         f = S()
         g = S()
         while f == 0 || g == 0 || is_constant(g)
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
   Qxy, (x, y) = QQ[:x, :y]
   @test_throws ArgumentError remove(0*x, y)
   @test_throws ArgumentError remove(x, y^0)
end

@testset "Generic.MPoly.derivative" begin
   for num_vars=1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      R, vars = polynomial_ring(ZZ, var_names; internal_ordering=ord)

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
end

@testset "Generic.MPoly.change_base_ring" begin
   F2, = residue_ring(ZZ, 2)
   R, varsR = polynomial_ring(F2, ["x"])
   S, varsS = polynomial_ring(R, ["y"])
   f = x -> x^2
   @test map_coefficients(f, varsR[1] * varsS[1]) == f(varsR[1]) * varsS[1]

   for num_vars=1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      R, vars = polynomial_ring(ZZ, var_names; internal_ordering=ord)

      F2x, varss = polynomial_ring(F2, var_names; internal_ordering = ord)

      @test typeof(AbstractAlgebra.Generic.change_base_ring(ZZ, R(1))) == AbstractAlgebra.Generic.MPoly{typeof(ZZ(1))}
      @test typeof(AbstractAlgebra.Generic.change_base_ring(ZZ, R(0))) == AbstractAlgebra.Generic.MPoly{typeof(ZZ(0))}

      for iter in 1:10
         f = rand(R, 5:10, 1:10, -100:100)
         @test evaluate(change_base_ring(R, f), [one(R) for i=1:num_vars]) == sum(f.coeffs[i] for i=1:f.length)
         @test evaluate(change_base_ring(R, f), vars) == f
         @test internal_ordering(parent(change_base_ring(R, f))) == internal_ordering(parent(f))

         g = change_base_ring(F2, f, parent = F2x)
         @test base_ring(g) === F2
         @test parent(g) === F2x

         g = map_coefficients(z -> z + 1, f, parent = R)
         @test parent(g) === R
      end
   end
end

@testset "Generic.MPoly.vars" begin
   for num_vars=1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      R, vars_R = polynomial_ring(ZZ, var_names; internal_ordering=ord)

      for iter in 1:10
         f = rand(R, 5:10, 1:10, -100:100)
         @test length(vars(R(evaluate(f, [one(R) for i=1:num_vars])))) == 0
      end
   end
end

@testset "Generic.MPoly.combine_like_terms" begin
   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      R, vars_R = polynomial_ring(ZZ, var_names; internal_ordering=ord)

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
end

@testset "Generic.MPoly.exponents" begin
   for num_vars = 1:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = rand_ordering()

      R, vars_R = polynomial_ring(ZZ, var_names; internal_ordering=ord)

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
end

@testset "Generic.MPoly.to_univariate" begin
   for num_vars=1:10
      ord = rand_ordering()

      var_names = ["x$j" for j in 1:num_vars]

      R, vars_R = polynomial_ring(ZZ, var_names; internal_ordering=ord)
      x = rand(vars_R)

      R_univ, x_univ = polynomial_ring(ZZ, "x")

      @test zero(R_univ) == to_univariate(R_univ, zero(R))
      @test one(R_univ) == to_univariate(R_univ, one(R))

      p = to_univariate(vars_R[1])
      Rp = parent(p)

      @test string(symbols(Rp)[1]) == var_names[1]

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
end

@testset "Generic.MPoly.coefficients_of_univariate" begin
   for num_vars=1:10
      ord = rand_ordering()
      var_names = ["x$j" for j in 1:num_vars]

      R, vars_R = polynomial_ring(ZZ, var_names; internal_ordering=ord)

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
end

@testset "Generic.MPoly.ordering" begin
   n_mpolys = 100
   maxval = 10
   maxdeg = 20

   # :deglex ordering
   R, (x,y,z) = AbstractAlgebra.Generic.polynomial_ring(AbstractAlgebra.Generic.ZZ, ["x", "y", "z"], internal_ordering=:lex)
   # Monomials of degree 2
   @test isless(z^2, y*z) == true
   @test isless(y*z, y^2) == true
   @test isless(y^2, x*z) == true
   @test isless(x*z, x*y) == true
   @test isless(x*y, x^2) == true

   for n_vars = 1:maxdeg
      A = unique(sortslices(reshape(map(Int,map(round, rand(n_vars * n_mpolys) * maxval)), (n_mpolys, n_vars)), dims=1),dims=1)
      var_names = ["x$j" for j in 1:n_vars]
      R, varsR = AbstractAlgebra.Generic.polynomial_ring(AbstractAlgebra.Generic.ZZ, var_names, internal_ordering=:lex)
      for i in 1:size(A)[1]-1
         f = R([base_ring(R)(1)], [A[i,:]])
         g = R([base_ring(R)(1)], [A[i+1,:]])
         @test isless(f,g)
      end
   end

   # :deglex ordering
   R, (x,y,z) = AbstractAlgebra.Generic.polynomial_ring(AbstractAlgebra.Generic.ZZ, ["x", "y", "z"], internal_ordering=:deglex)

   @test isless(z^2, y*z) == true
   @test isless(y*z, x*z) == true
   @test isless(y^2, x*z) == true
   @test isless(y^2, x*y) == true
   @test isless(x*y, x^2) == true

   for n_vars=1:maxdeg
      A = reshape(map(Int,map(round, rand(n_vars * n_mpolys) * maxval)), (n_mpolys, n_vars))
      var_names = ["x$j" for j in 1:n_vars]
      R, varsR = AbstractAlgebra.Generic.polynomial_ring(AbstractAlgebra.Generic.ZZ, var_names, internal_ordering=:deglex)

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
   R, (x,y,z) = AbstractAlgebra.Generic.polynomial_ring(AbstractAlgebra.Generic.ZZ, ["x", "y", "z"], internal_ordering=:degrevlex)
   # Monomials of degree 2
   @test isless(z^2, y*z) == true
   @test isless(y*z, x*z) == true
   @test isless(x*z, y^2) == true
   @test isless(y^2, x*y) == true
   @test isless(x*y, x^2) == true
   for n_vars = 1:maxdeg
      A = reshape(map(Int,map(round, rand(n_vars * n_mpolys) * maxval)), (n_mpolys, n_vars))
      var_names = ["x$j" for j in 1:n_vars]
      R, varsR = AbstractAlgebra.Generic.polynomial_ring(AbstractAlgebra.Generic.ZZ, var_names, internal_ordering=:degrevlex)
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
end

@testset "Generic.MPoly.equality" begin
    R, _ = polynomial_ring(ZZ, Symbol[])
    @test R(1) != R(2)
end

@testset "Generic.MPoly.exceptions" begin
   @test_throws MethodError polynomial_ring(Char, [:x])
   @test_throws Exception mpoly_ring_type(Char)
   @test_throws ArgumentError mpoly_type(Char)
end

@testset "Generic.MPoly.as_univariate" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])
   f = x * y + 2 * x^2 * y + x * y^2
   @test coefficients(f, 1) == typeof(f)[R(), y + y^2, 2 * y]
   @test leading_coefficient(f, 1) == coefficients(f, 1)[end]
   @test content(f, 1) == y
end

@testset "Generic.MPoly.Ring_interface" begin
  S, = polynomial_ring(QQ, 0)
  ConformanceTests.test_Ring_interface_recursive(S)

  S, = polynomial_ring(QQ, 1)
  ConformanceTests.test_Ring_interface_recursive(S)

  S, = polynomial_ring(ZZ, 2)
  ConformanceTests.test_Ring_interface_recursive(S)

  R, = QQ[:x]
  S, = polynomial_ring(R, :z => 1:3)
  ConformanceTests.test_Ring_interface(S) # _recursive needs too many ressources
end

@testset "Generic.MPoly.zero_rings" begin
  R, = residue_ring(ZZ, 1)
  S, = polynomial_ring(R, 2)
  @test is_zero(gen(S, 1)) && is_one(gen(S, 1))
  @test is_zero(one(S))
  ConformanceTests.test_Ring_interface_recursive(S)
end

# -------------------------------------------------------

# Coeff rings needed for the tests below
ZeroRing,_ = residue_ring(ZZ,1);
ZZmod720,_ = residue_ring(ZZ, 720);

## MPoly over ZeroRing
@testset "Nilpotent/unit for ZeroRing[x,y]" begin
  P,(x,y) = polynomial_ring(ZeroRing, ["x", "y"]);
  @test is_nilpotent(P(0))
  @test is_nilpotent(P(1))
  @test is_nilpotent(x)
  @test is_nilpotent(-x)
  @test is_nilpotent(x+y)
  @test is_nilpotent(x-y)
  @test is_nilpotent(x*y)

  @test is_unit(P(0))
  @test is_unit(P(1))
  @test is_unit(x)
  @test is_unit(-x)
  @test is_unit(x+y)
  @test is_unit(x-y)
  @test is_unit(x*y)
end

## MPoly over ZZ
@testset "Nilpotent/unit for ZZ[x,y]" begin
  P,(x,y) = polynomial_ring(ZZ, ["x", "y"]);
  @test is_nilpotent(P(0))
  @test !is_nilpotent(P(1))
  @test !is_nilpotent(x)
  @test !is_nilpotent(-x)
  @test !is_nilpotent(x+y)
  @test !is_nilpotent(x-y)
  @test !is_nilpotent(x*y)

  @test !is_unit(P(0))
  @test is_unit(P(1))
  @test is_unit(P(-1))
  @test !is_unit(P(-2))
  @test !is_unit(P(-2))
  @test !is_unit(x)
  @test !is_unit(-x)
  @test !is_unit(x+1)
  @test !is_unit(x-1)
  @test !is_unit(x+y)
  @test !is_unit(x-y)
  @test !is_unit(x*y)
end

## MPoly over QQ
@testset "Nilpotent/unit for QQ[x,y]" begin
  P,(x,y) = polynomial_ring(QQ, ["x", "y"]);
  @test is_nilpotent(P(0))
  @test !is_nilpotent(P(1))
  @test !is_nilpotent(x)
  @test !is_nilpotent(-x)
  @test !is_nilpotent(x+y)
  @test !is_nilpotent(x-y)
  @test !is_nilpotent(x*y)

  @test !is_unit(P(0))
  @test is_unit(P(1))
  @test is_unit(P(-1))
  @test is_unit(P(-2))
  @test is_unit(P(-2))
  @test !is_unit(x)
  @test !is_unit(-x)
  @test !is_unit(x+1)
  @test !is_unit(x-1)
  @test !is_unit(x+y)
  @test !is_unit(x-y)
  @test !is_unit(x*y)
end

## MPoly over ZZ/720
@testset "Nilpotent/unit for ZZ/(720)[x,y]" begin
  FactorsOf30 = [2, 3, 5, 6, 10, 15]; # non-nilpotent zero-divisors
  P,(x,y) = polynomial_ring(ZZmod720, ["x", "y"]);
  @test is_nilpotent(P(0))
  for ZeroDiv in FactorsOf30
    @test !is_nilpotent(P(ZeroDiv))
  end
  @test is_nilpotent(P(30))
  @test !is_nilpotent(x)
  @test !is_nilpotent(-x)
  for ZeroDiv in FactorsOf30
    @test !is_nilpotent(ZeroDiv*x)
  end
  @test is_nilpotent(30*x)
  @test is_nilpotent(30*x+120*y)
  @test is_nilpotent(30*x-120*y)
  @test !is_nilpotent(x*y)

  @test !is_unit(P(0))
  @test is_unit(P(1))
  @test is_unit(P(-1))
  @test !is_unit(P(2))
  @test !is_unit(P(-2))
  @test is_unit(P(7))
  @test is_unit(P(-7))
  @test !is_unit(x)
  @test !is_unit(-x)
  @test !is_unit(x+1)
  for ZeroDiv in FactorsOf30
    @test !is_unit(x+ZeroDiv)
  end
  @test !is_unit(x+30)
  @test !is_unit(x-30)
  @test is_unit(1+30*x)
  for ZeroDiv in FactorsOf30
    @test !is_unit(1+x*ZeroDiv)
  end
  @test is_unit(7+60*x)
  @test is_unit(7-60*x)
  @test is_unit(1+30*(x+y))
  for ZeroDiv in FactorsOf30
    @test !is_unit(1+ZeroDiv*(x+y))
  end
  @test is_unit(7+60*x*y)
  @test is_unit(7-60*x*y)
end

@testset "Generic.MPoly.Issue#2010" begin
  R1, (y1, z1) = polynomial_ring(QQ, [:y, :z])
  R2, (x2, y2) = polynomial_ring(QQ, [:x, :y])
  @test_throws ErrorException z1 + y2
end
