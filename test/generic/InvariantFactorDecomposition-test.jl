function test_invariant_factors_constructors()
   print("Generic.InvariantFactorDecomposition.constructors...")

   for R in [ZZ, QQ]
      for iter = 1:100
         M = rand_module(R, -10:10)

         I, f = snf(M)

         @test isa(I, Generic.SNFModule)
      end
   end

   F = FreeModule(ZZ, 0)
   D, f = snf(F)
   m = D([])

   @test isa(m, Generic.snf_module_elem)

   println("PASS")
end

function test_invariant_factors()
   print("Generic.InvariantFactorDecomposition.invariant_factors...")

   for R in [ZZ, QQ]
      for iter = 1:100
         M = rand_module(R, -10:10)

         I, f = snf(M)
         invs = invariant_factors(M)

         @test invariant_factors(I) == invs
      end
   end

   println("PASS")
end

function test_invariant_factors_isomorphism()
   print("Generic.InvariantFactorDecomposition.isomorphism...")

   for R in [ZZ, QQ]
      for iter = 1:100
         M = rand_module(R, -10:10)

         I, f = snf(M)

         K, g = kernel(f)

         @test length(invariant_factors(K)) == 0

         m = rand(I, -10:10)

         @test m == inv(f)(f(m))
      end
   end
   
   println("PASS")
end

function test_invariant_factor_decomposition()
   test_invariant_factors_constructors()
   test_invariant_factors()
   test_invariant_factors_isomorphism()

   println("")
end
