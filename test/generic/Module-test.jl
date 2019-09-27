function rand_homomorphism(M::AbstractAlgebra.FPModule{T}, vals...) where T <: RingElement
   rk = rand(1:5)
   m = ngens(M)
   R = base_ring(M)
   F = FreeModule(R, rk)
   S = MatrixSpace(R, rk, m)
   mat = rand(S, vals...)
   f = ModuleHomomorphism(F, M, mat)
   ngens1 = rand(1:3)
   gens1 = [rand(F, vals...) for j in 1:ngens1]
   S, g = sub(F, gens1)
   hom1 = compose(g, f)
   return S, hom1
end

@testset "Generic.Module.rand..." begin
   F = FreeModule(ZZ, 3)
   f = rand(F, 1:9)
   @test f isa Generic.FreeModuleElem
   f = rand(rng, F, 1:9)
   @test f isa Generic.FreeModuleElem
end

@testset "Generic.Module.manipulation..." begin
   for R in [ZZ, QQ]
      for iter = 1:100
         F = FreeModule(R, 3)
         M = rand_module(R, -10:10)

         ngens1 = rand(1:5)
         gens1 = [rand(M, -10:10) for j in 1:ngens1]
         M1, f1 = sub(M, gens1)

         ngens2 = rand(1:5)
         gens2 = [rand(M1, -10:10) for j in 1:ngens2]
         M2, f2 = sub(M1, gens2)

         Q, g = quo(M, M1)

         @test issubmodule(M, M)
         @test issubmodule(M, M1)
         @test issubmodule(M, M2)
         @test !issubmodule(M, Q)
         @test !issubmodule(M1, F)
         @test !issubmodule(M1, M)

         ngens3 = rand(1:5)
         gens3 = [rand(M, -10:10) for j in 1:ngens3]
         M3, f3 = sub(M, gens3)

         ngens4 = rand(1:5)
         gens4 = [rand(M3, -10:10) for j in 1:ngens4]
         M4, f4 = sub(M3, gens4)

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
   S, f = sub(F, [m, n])

   Q, g = quo(F, S)

   b = Q(BigInt[2, 0, 4])
   c = Q(BigInt[7, 1, -59])
   d = Q(BigInt[4, 0, 23])
   T, h = sub(Q, [b, c, d])

   U, k = quo(Q, T)
end

@testset "Generic.Module.elem_getindex..." begin

   for R in [ZZ, QQ]
      for iter = 1:100
         F = FreeModule(R, 3)
         M = rand_module(R, -10:10)

         m1 = rand(M, -10:10)

         m2 = zero(M)

         if ngens(M) != 0
            m2 = sum(m1[i]*gen(M, i) for i in 1:ngens(M))
         end

         @test m1 == m2
      end
   end
end

@testset "Generic.Module.intersection..." begin
   for R in [ZZ, QQ]
      for iter = 1:100
         M = rand_module(R, -10:10)

         ngens1 = rand(1:5)
         gens1 = [rand(M, -10:10) for j in 1:ngens1]
         M1, f1 = sub(M, gens1)

         ngens2 = rand(1:5)
         gens2 = [rand(M, -10:10) for j in 1:ngens2]
         M2, f2 = sub(M, gens2)

         ngens3 = rand(1:5)
         gens3 = [rand(M, -10:10) for j in 1:ngens3]
         M3, f3 = sub(M, gens3)

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
end

@testset "Generic.Module.is_isomorphic..." begin
   # Test the first isomorphism theorem
   for R in [ZZ, QQ]
      for iter = 1:100
         M = rand_module(R, -10:10)
         S, f = rand_homomorphism(M, -10:10)

         K, g = kernel(f)
         Q, h = quo(S, K)

         I, k = image(f)

         @test isisomorphic(Q, I)
      end
   end

   # Test submodules are isomorphic to their image
   for R in [ZZ, QQ]
      for iter = 1:100
         M = rand_module(R, -10:10)

         ngens1 = rand(1:5)
         gens1 = [rand(M, -10:10) for j in 1:ngens1]
         M1, f1 = sub(M, gens1)

         I, g = image(f1)

         @test isisomorphic(I, M1)
      end
   end
end

@testset "Generic.Module.coercions..." begin
   # Test the first isomorphism theorem
   for R in [ZZ, QQ]
      for iter = 1:20
         M = rand_module(R, -10:10)

         S2, f2 = sub(M, [rand(M, -10:10)])
         S1, f1 = sub(S2, [rand(S2, -10:10)])

         SQ1, g1 = sub(M, [rand(M, -10:10)])
         Q1, h1 = quo(M, SQ1)
         SQ2, g2 = sub(Q1, [rand(Q1, -10:10)])
         Q2, h2 = quo(Q1, SQ2)
         SQ3, g3 = sub(Q2, [rand(Q2, -10:10)])
         Q3, h3 = quo(Q2, SQ3)

         Sa, k1 = sub(Q2, [rand(Q2, -10:10)])
         Sb, k2 = sub(Sa, [rand(Sa, -10:10)])

         m1 = rand(M, -10:10)
         @test parent(M(m1)) === M
         @test parent(Q1(m1)) === Q1
         @test parent(Q2(m1)) === Q2

         m2 = rand(Q1, -10:10)
         @test parent(Q1(m2)) == Q1
         @test parent(Q2(m2)) == Q2

         m3 = rand(S2, -10:10)
         @test parent(M(m3)) === M
         @test parent(Q1(m3)) === Q1
         @test parent(Q2(m3)) === Q2

         m4 = rand(S1, -10:10)
         @test parent(M(m4)) === M
         @test parent(Q1(m4)) === Q1
         @test parent(Q2(m4)) === Q2

         m5 = rand(Sa, -10:10)
         @test parent(Sa(m5)) === Sa
         @test parent(Q2(m5)) === Q2
         @test parent(Q3(m5)) === Q3

         m6 = rand(Sb, -10:10)
         @test parent(Sa(m6)) === Sa
         @test parent(Q2(m6)) === Q2
         @test parent(Q3(m6)) === Q3

         m7 = rand(S2, -10:10)
         @test parent(M(m7)) === M
         @test parent(Q1(m7)) === Q1
         @test parent(Q2(m7)) === Q2

         m8 = rand(S1, -10:10)
         @test parent(M(m7)) === M
         @test parent(Q1(m7)) === Q1
         @test parent(Q2(m7)) === Q2

         D, v = DirectSum(S1, Sa)
         m9 = rand(D, -10:10)
         @test parent(D(m9)) === D
      end
   end
end
