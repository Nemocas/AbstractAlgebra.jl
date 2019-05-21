function rand_module(R::AbstractAlgebra.Ring, vals...)
   rk = rand(1:5)
   M = FreeModule(R, rk)
   levels = rand(0:3)
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

function test_module_manipulation()
   print("Generic.Module.manipulation...")

   for R in [ZZ, QQ]
      for iter = 1:100
         F = FreeModule(R, 3)
         M = rand_module(R, -10:10)

         ngens1 = rand(1:5)
         gens1 = [rand(M, -10:10) for j in 1:ngens1]
         M1, f1 = Submodule(M, gens1)

         ngens2 = rand(1:5)
         gens2 = [rand(M1, -10:10) for j in 1:ngens2]
         M2, f2 = Submodule(M1, gens2)

         Q, g = QuotientModule(M, M1)

         @test issubmodule(M, M)
         @test issubmodule(M, M1)
         @test issubmodule(M, M2)
         @test !issubmodule(M, Q)
         @test !issubmodule(M1, F)
         @test !issubmodule(M1, M)

         ngens3 = rand(1:5)
         gens3 = [rand(M, -10:10) for j in 1:ngens3]
         M3, f3 = Submodule(M, gens3)

         ngens4 = rand(1:5)
         gens4 = [rand(M3, -10:10) for j in 1:ngens4]
         M4, f4 = Submodule(M3, gens4)

         flag, U = iscompatible(M1, M3)
         @test flag == true && U == M

         flag, U = iscompatible(M1, M4)
         @test flag == true && U == M

         flag, U = iscompatible(M2, M3)
         @test flag == true && U == M

         flag, U = iscompatible(M2, M4)
         @test flag == true && U == M

         flag, U = iscompatible(M, M)
         @test flag == true && U == M

         flag, U = iscompatible(M1, Q) 
         @test flag == false

         flag, U = iscompatible(M2, Q)
         @test flag == false

         flag, U = iscompatible(Q, M1)
         @test flag == false

         flag, U = iscompatible(Q, M2)
         @test flag == false
      end
   end

   # Regression test: triggers crash if Submodule
   # doesn't doesn't reduce generators mod relations
   # and put in hnf

   F = FreeModule(ZZ, 3)

   m = F(BigInt[-10, -5, -1])
   n = F(BigInt[0, -2, -9])
   S, f = Submodule(F, [m, n])

   Q, g = QuotientModule(F, S)

   b = Q(BigInt[2, 0, 4])
   c = Q(BigInt[7, 1, -59])
   d = Q(BigInt[4, 0, 23])
   T, h = Submodule(Q, [b, c, d])

   U, k = QuotientModule(Q, T)

   println("PASS")
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

         @test I2 == J2 # test associativity

         I3, g3 = intersect(M2, M1)

         @test I1 == I3 # test commutativity

         I4, g4 = intersect(M1, M1)

         @test I4 == M1 # test reflexivity

         I5, g5 = intersect(I1, M2)
         I6, g6 = intersect(I1, M1)

         @test I1 == I5 # test absorption
         @test I1 == I6
      end
   end
   
   println("PASS")
end

function test_module()
   test_module_manipulation()
   test_module_intersection()

   println("")
end
