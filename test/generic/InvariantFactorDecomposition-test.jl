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

function test_isomorphism()
   print("Generic.InvariantFactorDecomposition.isomorphism...")

   for R in [ZZ, QQ]
      for iter = 1:100
         M = rand_module(R, -10:10)

         I, f = snf(M)

         K, g = kernel(f)

         @test length(invariant_factors(K)) == 0
      end
   end
   
   println("PASS")
end

function test_invariant_factor_decomposition()
   test_invariant_factors()
   test_isomorphism()

   println("")
end
