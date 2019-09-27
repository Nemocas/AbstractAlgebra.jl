@testset "Generic.InvariantFactorDecomposition.constructors..." begin
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

   @test isa(m, Generic.SNFModuleElem)
end

@testset "Generic.InvariantFactorDecomposition.invariant_factors..." begin
   for R in [ZZ, QQ]
      for iter = 1:100
         M = rand_module(R, -10:10)

         I, f = snf(M)
         invs = invariant_factors(M)

         @test invariant_factors(I) == invs
      end
   end
end

@testset "Generic.InvariantFactorDecomposition.isomorphism..." begin
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
end
