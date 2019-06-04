function test_invariant_factors()
   print("Generic.InvariantFactorDecomposition.invariant_factors...")

   for R in [ZZ, QQ]
      for iter = 1:100
         M = rand_module(R, -10:10)

         I, f = InvariantFactorDecomposition(M)
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

         I, f = InvariantFactorDecomposition(M)

         K, g = kernel(f)

         @test ncols(invariant_factors(K)) == 0
      end
   end
   
   println("PASS")
end

function test_invariant_factor_decomposition()
   test_invariant_factors()
   test_isomorphism()

   println("")
end
