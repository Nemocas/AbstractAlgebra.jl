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

function test_gen_mpoly()
   test_gen_mpoly_constructors()

   println("")
end
