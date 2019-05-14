function rand_module(R::AbstractAlgebra.Ring, vals...)
   rk = rand(1:5)
   M = FreeModule(R, rk)
   levels = rand(0:1)
   for i = 1:levels
      if ngens(M) == 0
         break
      end
      G = [rand(M, vals...) for i in 1:rand(1:ngens(M))]
      S, f = Submodule(M, G)
      if rand(1:2) == 1
         M, f = QuotientModule(M, S)
      else
         M = S
      end
   end
   return M
end

function test_module_intersection()
   print("Generic.Module.intersection...")

   for R in [ZZ, QQ]
      for iter = 1:100
         M = rand_module(R, -10:10)

         ngens1 = rand(1:5)
         gens1 = [rand(M, -10:10) for j in 1:ngens1]
         M1, f1 = Submodule(M, gens1)

         ngens2 = rand(1:5)
         gens2 = [rand(M, -10:10) for j in 1:ngens2]
         M2, f2 = Submodule(M, gens2)

         ngens3 = rand(1:5)
         gens3 = [rand(M, -10:10) for j in 1:ngens3]
         M3, f3 = Submodule(M, gens3)

         I1, g1 = intersect(M1, M2)
         I2, g2 = intersect(I1, M3)

         J1, h1 = intersect(M2, M3)
         J2, h2 = intersect(M1, J1)

         @test I2 == J2
      end
   end
   
   println("PASS")
end

function test_module()
   test_module_intersection()

   println("")
end
