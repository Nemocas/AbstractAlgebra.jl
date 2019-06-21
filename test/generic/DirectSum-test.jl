using Random

function test_direct_sum_constructors()
   print("Generic.DirectSum.constructors...")

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
   S, f = Submodule(F, [m, n])
   Q, g = QuotientModule(F, S)
   D, h = DirectSum(Q)
   m = D([])

   @test isa(m, Generic.direct_sum_module_elem)

   println("PASS")
end

function test_direct_sum_basic_manipulation()
   print("Generic.DirectSum.basic_manipulation...")

   for R in [ZZ, QQ]
      for iter = 1:100
         num = rand(1:5)
         M = [rand_module(R, -10:10) for i in 1:num]

         D, f, g = DirectSum(M...)

         @test summands(D) == M
      end
   end

   println("PASS")
end

function test_direct_sum_maps()
   print("Generic.DirectSum.maps...")

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

   println("PASS")
end

function test_direct_sum_isomorphism()
   print("Generic.DirectSum.isomorphism...")

   for R in [ZZ, QQ]
      for iter = 1:100
         num = rand(1:5)
         M = [rand_module(R, -10:10) for i in 1:num]

         D1, f, g = DirectSum(M)

         D2, f, g = DirectSum(shuffle(M))

         @test isisomorphic(D1, D2)
      end
   end
   
   println("PASS")
end

function test_direct_sum()
   test_direct_sum_constructors()
   test_direct_sum_basic_manipulation()
   test_direct_sum_maps()
   test_direct_sum_isomorphism()

   println("")
end
