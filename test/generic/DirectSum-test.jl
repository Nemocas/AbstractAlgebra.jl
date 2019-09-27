using Random

@testset "Generic.DirectSum.constructors..." begin
   for R in [ZZ, QQ]
      for iter = 1:100
         num = rand(1:5)
         M = [rand_module(R, -10:10) for i in 1:num]
         D, f, g = DirectSum(M)

         @test isa(D, Generic.DirectSumModule)
      end
   end

   F = FreeModule(ZZ, 2)
   m = F(BigInt[1, 0])
   n = F(BigInt[0, 1])
   S, f = sub(F, [m, n])
   Q, g = quo(F, S)
   D, h = DirectSum(Q)
   m = D([])

   @test isa(m, Generic.DirectSumModuleElem)
end

@testset "Generic.DirectSum.basic_manipulation..." begin

   for R in [ZZ, QQ]
      for iter = 1:100
         num = rand(1:5)
         M = [rand_module(R, -10:10) for i in 1:num]

         D, f, g = DirectSum(M...)

         @test summands(D) == M
      end
   end
end

@testset "Generic.DirectSum.maps..." begin
   for R in [ZZ, QQ]
      for iter = 1:100
         num = rand(1:5)
         M = [rand_module(R, -10:10) for i in 1:num]

         D, f, g = DirectSum(M)

         m = [rand(M[i], -10:10) for i in 1:num]

         for i = 1:num
            @test g[i](f[i](m[i])) == m[i]
         end
      end
   end
end

@testset "Generic.DirectSum.isomorphism..." begin
   for R in [ZZ, QQ]
      for iter = 1:100
         num = rand(1:5)
         M = [rand_module(R, -10:10) for i in 1:num]

         D1, f, g = DirectSum(M)

         D2, f, g = DirectSum(shuffle(M))

         @test isisomorphic(D1, D2)
      end
   end
end
