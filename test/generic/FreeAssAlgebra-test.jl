@testset "Generic.FreeAssAlgebra.constructors" begin
   R, x = ZZ["y"]

   for num_vars = 1:5
      var_names = ["x$j" for j in 1:num_vars]

      S, varlist = free_associative_algebra(R, var_names)

      @test free_associative_algebra(R, var_names, cached = true)[1] === free_associative_algebra(R, var_names, cached = true)[1]
      @test free_associative_algebra(R, var_names, cached = false)[1] !== free_associative_algebra(R, var_names, cached = true)[1]
      @test free_associative_algebra(R, num_vars, "x", cached = true)[1] === free_associative_algebra(R, var_names, cached = true)[1]
      @test free_associative_algebra(R, num_vars, cached = true)[1] === free_associative_algebra(R, var_names, cached = true)[1]

      @test elem_type(S) == Generic.FreeAssAlgElem{elem_type(R)}
      @test elem_type(Generic.FreeAssAlgebra{elem_type(R)}) == Generic.FreeAssAlgElem{elem_type(R)}
      @test parent_type(Generic.FreeAssAlgElem{elem_type(R)}) == Generic.FreeAssAlgebra{elem_type(R)}
      @test base_ring(S) === R
      @test coefficient_ring(S) === R
      @test ngens(S) == length(gens(S))

      @test typeof(S) <: Generic.FreeAssAlgebra

      isa(symbols(S), Vector{Symbol})

      for j = 1:num_vars
         @test coefficient_ring(varlist[j]) === R
         @test isa(varlist[j], FreeAssAlgElem)
         @test isa(gens(S)[j], FreeAssAlgElem)
      end

      f = rand(S, 0:5, 0:10, 0:0, -100:100)

      @test isa(f, FreeAssAlgElem)

      @test isa(S(2), FreeAssAlgElem)

      @test isa(S(R(2)), FreeAssAlgElem)

      @test isa(S(f), FreeAssAlgElem)

      V = [R(rand(-100:100)) for i in 1:5]
      W = [[rand(1:num_vars) for i in 1:rand(0:9)] for j in 1:5]
      f1 = S(V, W)
      @test isa(f1, FreeAssAlgElem)

      f2 = S()
      fit!(f2, 5)

      for i in 1:5
         if iseven(num_vars)
            f2 = set_exponent_word!(f2, i, W[i])
            f2 = setcoeff!(f2, i, V[i])
         else
            f2 = setcoeff!(f2, i, V[i])
            f2 = set_exponent_word!(f2, i, W[i])
         end
      end
      f2 = sort_terms!(f2)
      f2 = combine_like_terms!(f2)

      @test f1 == f2

      C = MPolyBuildCtx(S)

      @test string(C) isa String

      for i in 1:5
         push_term!(C, V[i], W[i])
      end
      f3 = finish(C)

      @test f1 == f3

      f4 = zero(S)
      for t in terms(f1)
         f4 = addeq!(f4, t)
      end
      @test f1 == f4

      f4 = zero(S)
      for t in zip(coefficients(f1), monomials(f1))
         f4 = addmul!(f4, S(t[1]), t[2], S())
      end
      @test f1 == f4

      if !iszero(f1)
         @test leading_term(f1) == leading_coefficient(f1)*leading_monomial(f1)
         @test total_degree(f1) >= total_degree(f1 - leading_term(f1))
      else
         @test_throws ArgumentError leading_term(f1)
         @test_throws ArgumentError leading_monomial(f1)
         @test_throws ArgumentError leading_exponent_word(f1)
      end

      @test canonical_unit(f1) == canonical_unit(leading_coefficient(f1))

      @test !is_gen(zero(S))
      @test !is_gen(one(S))
      for i in 1:num_vars
         g = gen(S, i)
         @test g == S[i]
         @test is_gen(g)
         @test !is_gen(g + 1)
         @test leading_exponent_word(g) == [i]
      end

      @test_throws ArgumentError gen(S, 0)
      @test_throws ArgumentError gen(S, num_vars + 1)

      @test_throws ArgumentError leading_term(zero(S))
      @test_throws ArgumentError leading_monomial(zero(S))
      @test_throws ArgumentError leading_exponent_word(zero(S))

      @test collect(exponent_words(varlist[1] + 1)) == [Int[1], Int[]]
      @test isone(varlist[1]^0)

      _, varlist = polynomial_ring(QQ, var_names)
      y = varlist[1]
      @test x in [x, y]
      @test x in [y, x]
      @test !(x in [y])
      @test x in keys(Dict(x => 1))
      @test !(y in keys(Dict(x => 1)))
   end
end

@testset "Generic.FreeAssAlgebra.printing" begin
   R, x = ZZ["y"]
   S = free_associative_algebra(R, 5)
   @test !occursin("\n", sprint(show, R))
end

function test_elem(R::Generic.FreeAssAlgebra{elem_type(ZZ)})
   return rand(R, 0:4, 0:5, -10:10)
end

@testset "Generic.FreeAssAlgebra.change_base_ring" begin
   F5, = residue_ring(ZZ, 5)
   R, varsR = polynomial_ring(F5, ["x"])
   S, varsS = free_associative_algebra(R, ["y"])
   f = x -> x^2 + F5(3)
   @test map_coefficients(f, varsR[1] * varsS[1]) == f(varsR[1]) * varsS[1]

   for num_vars = 1:5
      var_names = ["x$j" for j in 1:num_vars]

      R, t = ZZ["t"]
      Rx, varsRx = free_associative_algebra(R, var_names)
      S, _ = polynomial_ring(R, ["y", "z"])
      Sx, varsSx = free_associative_algebra(S, var_names)

      @test typeof(change_base_ring(R, Rx(0))) == typeof(Rx(0))
      @test typeof(change_base_ring(R, Rx(1))) == typeof(Rx(1))
      @test typeof(change_base_ring(S, Rx(0))) == typeof(Sx(0))
      @test typeof(change_base_ring(S, Rx(1))) == typeof(Sx(1))

      @test change_base_ring(R, Rx(0)) == Rx(0)
      @test change_base_ring(R, Rx(1)) == Rx(1)
      @test change_base_ring(S, Rx(0)) == Sx(0)
      @test change_base_ring(S, Rx(1)) == Sx(1)

      # some pseudo-random non-trivial polynomial
      f = sum((5 + i) * x^i for (i, x) in enumerate(varsRx)) + t^2 * prod(varsRx) - t + 42
      g = sum((5 + i) * x^i for (i, x) in enumerate(varsSx)) + S(t)^2 * prod(varsSx) - S(t) + 42
      @test change_base_ring(S, f, parent=Sx) == g

      for _ in 1:10
         f1 = rand(Rx, 0:5, 0:10, 0:0, -100:100)
         f2 = rand(Rx, 0:5, 0:10, 0:0, -100:100)

         g = change_base_ring(R, f1)
         @test base_ring(g) === R
         @test g == f1

         g = change_base_ring(S, f1, parent=Sx)
         @test base_ring(g) === S
         @test parent(g) === Sx

         @test change_base_ring(S, f1 + f2) == change_base_ring(S, f1) + change_base_ring(S, f2)
         @test change_base_ring(S, f1 * f2) == change_base_ring(S, f1) * change_base_ring(S, f2)
      end
   end
end

@testset "Generic.FreeAssAlgebra.divexact" begin
   R, (x, y, z) = free_associative_algebra(ZZ, [:x, :x, :x])
   a = x + y + z
   @test divexact(2*a, 2) == a
   @test divexact_left(2*a, 2) == a
   @test divexact_right(a*2, 2) == a
   @test_throws ArgumentError divexact(2*a + 3, 2)
   @test_throws ArgumentError divexact_left(2*a + 3, R(2))
   @test_throws ArgumentError divexact_right(3 - 2*a, R(2))
   @test_throws ArgumentError divexact_left(a*x + 1, x)
   @test_throws ArgumentError divexact_right(1 + a*x, x)
   @test_throws ArgumentError divexact_left(x*a - 1, x)
   @test_throws ArgumentError divexact_right(1 - x*a, x)
end

@testset "Generic.FreeAssAlgebra.is_unit" begin
   R, (x,) = free_associative_algebra(residue_ring(ZZ, 4)[1], ["x"])

   @test !is_unit(x)
   @test !is_unit(2*x)
   try
      res = is_unit(1 + 2*x)
      @test res
   catch e
      @test e isa NotImplementedError
   end
end

@testset "Generic.FreeAssAlgebra.NCRing_interface" begin
   test_NCRing_interface(free_associative_algebra(ZZ, 3)[1])
end

