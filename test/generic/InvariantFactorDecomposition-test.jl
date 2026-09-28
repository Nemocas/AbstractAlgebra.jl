@testset "Generic.InvariantFactorDecomposition.constructors" begin
   for R in [ZZ, QQ]
      for iter = 1:100
         M = rand_module(R, -10:10)

         I, f = snf(M)

         @test isa(I, Generic.SNFModule)
      end
   end

   F = free_module(ZZ, 0)
   D, f = snf(F)
   m = D([])

   @test isa(m, Generic.SNFModuleElem)
end

@testset "Generic.InvariantFactorDecomposition.invariant_factors" begin
   for R in [ZZ, QQ]
      for iter = 1:100
         M = rand_module(R, -10:10)

         I, f = snf(M)
         invs = invariant_factors(M)

         @test invariant_factors(I) == invs
         if R == QQ
            @test rank(I) == ngens(I)
            @test dim(I) == ngens(I)
            @test vector_space_dim(I) == ngens(I)
         end
      end
   end

   F = free_module(ZZ, 1)
   S, _ = sub(F, [2 * gen(F, 1)])
   Q, _ = quo(F, S)
   I, _ = snf(Q)
   @test invariant_factors(I) == BigInt[2]
   @test_throws MethodError rank(I)
   @test_throws MethodError dim(I)
end

@testset "Generic.InvariantFactorDecomposition.isomorphism" begin
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
