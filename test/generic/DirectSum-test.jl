using Random

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
   test_direct_sum_basic_manipulation()
   test_direct_sum_maps()
   test_direct_sum_isomorphism()

   println("")
end
