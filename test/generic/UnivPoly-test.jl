@testset "Generic.UnivPoly.constructors" begin
   for R in [ZZ, QQ]
      for iters = 1:5
         ord = rand_ordering()

         @test UniversalPolynomialRing(R; ordering=ord, cached = true) === UniversalPolynomialRing(R; ordering=ord, cached = true)
         UniversalPolynomialRing(R; ordering=ord, cached = true) !== UniversalPolynomialRing(R; ordering=ord, cached = true)

         S = UniversalPolynomialRing(R; ordering=ord)

         x = gen(S, "x")

         @test isa(x, UnivPolyElem)

         y, z = gens(S, ["y", "z"])

         @test elem_type(S) == Generic.UnivPoly{elem_type(R), Generic.MPoly{elem_type(R)}}
         @test elem_type(Generic.UnivPolyRing{elem_type(R), Generic.MPoly{elem_type(R)}}) == Generic.UnivPoly{elem_type(R), Generic.MPoly{elem_type(R)}}
         @test parent_type(Generic.UnivPoly{elem_type(R), Generic.MPoly{elem_type(R)}}) == Generic.UnivPolyRing{elem_type(R), Generic.MPoly{elem_type(R)}}

         @test base_ring(S) === R
         @test coefficient_ring(S) === R

         @test typeof(S) <: Generic.UnivPolyRing

         @test isa(x, UnivPolyElem)
         @test isa(y, UnivPolyElem)
         @test isa(z, UnivPolyElem)

         f = rand(S, 0:5, 0:10, -10:10)

         @test isa(f, UnivPolyElem)
         @test isa(S(2), UnivPolyElem)
         @test isa(S(R(2)), UnivPolyElem)
         @test isa(S(f), UnivPolyElem)

         V = [R(rand(-10:10)) for i in 1:5]
         nused = rand(0:3)
         W = [[rand(0:100) for i in 1:nused] for i = 1:5]
         
         f1 = S(V, W)

         @test isa(f1, UnivPolyElem)

         f2 = S()
         fit!(f2, 5)

         for i = 1:5
            f2 = set_exponent_vector!(f2, i, W[i])
            f2 = setcoeff!(f2, i, V[i])
         end
         f2 = sort_terms!(f2)
         f2 = combine_like_terms!(f2)

         @test f1 == f2

         C = MPolyBuildCtx(S)

         for i = 1:5
            push_term!(C, V[i], W[i])
         end
         f3 = finish(C)

         @test f1 == f3

         @test x in [x, y]
         @test x in [y, x]
         @test y in [y, z]
         @test !(x in [y])
         @test !(y in [x])
         @test x in keys(Dict(x => 1))
         @test !(y in keys(Dict(x => 1)))
         @test !(y in keys(Dict(z => 1)))
      end
   end
end

@testset "Generic.UnivPoly.parent_type" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         ord = rand_ordering()

         S = UniversalPolynomialRing(R; ordering=ord, cached=false)

         @test symbols(S) == Symbol[]
         @test nvars(S) == 0

         f1 = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         f2 = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         f3 = rand(S, 0:5, 0:10, -10:10)

         @test isdomain_type(typeof(f1))
         @test isdomain_type(typeof(f2))
         @test isdomain_type(typeof(f3))

         @test isexact_type(typeof(f1))
         @test isexact_type(typeof(f2))
         @test isexact_type(typeof(f3))

         @test parent(f1) === S
         @test parent(f2) === S
         @test parent(f3) === S

         @test base_ring(S) === R
         @test base_ring(f1) === R
         @test base_ring(f2) === R
         @test base_ring(f3) === R

         @test nvars(S) == 3

         @test symbols(S) == [:x, :y, :z]

         @test length(vars(S())) == 0
         @test vars(x^2 + 2x + 1) == [x]
         @test vars(x*y + 1) == [x, y]
         @test vars(x*y + z) == [x, y, z]

         @test ordering(S) == ord
      end
   end
end


