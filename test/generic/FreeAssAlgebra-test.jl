@testset "Generic.FreeAssAlgebra.constructors" begin
   R, x = ZZ["y"]

   for num_vars = 1:5
      var_names = ["x$j" for j in 1:num_vars]

      S, varlist = FreeAssociativeAlgebra(R, var_names)

      @test FreeAssociativeAlgebra(R, var_names, cached = true)[1] === FreeAssociativeAlgebra(R, var_names, cached = true)[1]
      @test FreeAssociativeAlgebra(R, var_names, cached = false)[1] !== FreeAssociativeAlgebra(R, var_names, cached = true)[1]
      @test FreeAssociativeAlgebra(R, num_vars, "x", cached = true)[1] === FreeAssociativeAlgebra(R, var_names, cached = true)[1]
      @test FreeAssociativeAlgebra(R, num_vars, cached = true)[1] === FreeAssociativeAlgebra(R, var_names, cached = true)[1]

      @test elem_type(S) == Generic.FreeAssAlgElem{elem_type(R)}
      @test elem_type(Generic.FreeAssAlgebra{elem_type(R)}) == Generic.FreeAssAlgElem{elem_type(R)}
      @test parent_type(Generic.FreeAssAlgElem{elem_type(R)}) == Generic.FreeAssAlgebra{elem_type(R)}
      @test base_ring(S) === R
      @test coefficient_ring(S) === R

      @test typeof(S) <: Generic.FreeAssAlgebra

      isa(symbols(S), Vector{Symbol})

      for j = 1:num_vars
         @test isa(varlist[j], FreeAssAlgElem)
         @test isa(gens(S)[j], FreeAssAlgElem)
      end

      f =  rand(S, 0:5, 0:10, 0:0, -100:100)

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

      @test leading_term(f1) == leading_coefficient(f1)*leading_monomial(f1)
      @test total_degree(f1) >= total_degree(f1 - leading_term(f1))

      _, varlist = PolynomialRing(QQ, var_names)
      y = varlist[1]
      @test x in [x, y]
      @test x in [y, x]
      @test !(x in [y])
      @test x in keys(Dict(x => 1))
      @test !(y in keys(Dict(x => 1)))
   end
end

function test_elem(R::Generic.FreeAssAlgebra{elem_type(ZZ)})
   return rand(R, 0:4, 0:5, -10:10)
end

@testset "Generic.FreeAssAlgebra.divexact" begin
   R, (x, y, z) = FreeAssociativeAlgebra(ZZ, [:x, :x, :x])
   a = x + y + z
   @test divexact(2*a, 2) == a
   @test_throws ArgumentError divexact(2*a + 3, 2)
   @test_throws ArgumentError divexact_left(2*a + 3, R(2))
   @test_throws ArgumentError divexact_right(2*a + 3, R(2))
   @test_throws ArgumentError divexact_left(a*x + 1, x)
   @test_throws ArgumentError divexact_right(a*x + 1, x)
   @test_throws ArgumentError divexact_left(x*a + 1, x)
   @test_throws ArgumentError divexact_right(x*a + 1, x)
end

@testset "Generic.FreeAssAlgebra.NCRing_interface" begin
   test_NCRing_interface(FreeAssociativeAlgebra(ZZ, 3)[1])
end

