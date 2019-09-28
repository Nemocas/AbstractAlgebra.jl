@testset "Generic.ModuleHomomorphism.constructors..." begin
   M = FreeModule(ZZ, 2)

   f = ModuleHomomorphism(M, M, matrix(ZZ, 2, 2, [1, 2, 3, 4]))

   @test isa(f, Generic.Map(FunctionalMap))

   m = M([ZZ(1), ZZ(2)])

   @test f(m) == M([ZZ(7), ZZ(10)])

   N, nmap = sub(M, [m])

   g = ModuleHomomorphism(N, N, matrix(ZZ, 1, 1, [3]))

   m2 = N([ZZ(4)])

   @test g(m2) == N([ZZ(12)])
end

@testset "Generic.ModuleHomomorphism.kernel..." begin
   for R in [ZZ, QQ]
      for iter = 1:100
         # test kernels of canonical injection and projection
         M = rand_module(R, -10:10)
         ngens1 = rand(1:5)
         gens1 = [rand(M, -10:10) for j in 1:ngens1]
         M1, f1 = sub(M, gens1)

         Q, g = quo(M, M1)
         k1, h1 = kernel(f1)

         @test ngens(k1) == 0

         k2, h2 = kernel(g)
         for i = 1:3
            v = rand(k2, -10:10)

            @test iszero(g(h2(v)))
         end

         @test k2 == M1
      end
   end

   for R in [ZZ, QQ]
      for iter = 1:100
         # test random homomorphisms from free modules
         M = rand_module(R, -10:10)

         m = rand(1:5)
         F = FreeModule(R, m)

         S = MatrixSpace(R, m, ngens(M))
         f = ModuleHomomorphism(F, M, rand(S, -10:10))

         k, h = kernel(f)

         for i = 1:3
            v = rand(k, -10:10)

            @test iszero(f(h(v)))
         end
      end
   end
end

@testset "Generic.ModuleHomomorphism.image..." begin
   # To make it work on julia nightlies

   R = AbstractAlgebra.JuliaZZ
   for iter = 1:100
      # test image of composition of canonical injection and projection
      M = rand_module(R, -10:10)::AbstractAlgebra.FPModule{elem_type(R)}

      ngens1 = rand(1:5)
      gens1 = [rand(M, -10:10) for j in 1:ngens1]
      M1, f1 = sub(M, gens1)

      Q, g = quo(M, M1)
      k1, h1 = kernel(g)

      k = compose(h1, g)
      I, f = image(k)
      T, t = sub(Q, elem_type(Q)[])

      @test I == T
   end

   S = AbstractAlgebra.JuliaQQ
   for iter = 1:100
      # test image of composition of canonical injection and projection
      N = rand_module(S, -10:10)::AbstractAlgebra.FPModule{elem_type(S)}

      ngens1 = rand(1:5)
      gens1 = [rand(N, -10:10) for j in 1:ngens1]
      M1, f1 = sub(N, gens1)

      Q, g = quo(N, M1)
      k1, h1 = kernel(g)

      k = compose(h1, g)
      I, f = image(k)
      T, t = sub(Q, elem_type(Q)[])

      @test I == T
   end
end

@testset "Generic.ModuleIsomorphism..." begin
   R = AbstractAlgebra.JuliaQQ
   for iter = 1:100
      # test image of composition of canonical injection and projection
      M = rand_module(R, -10:10)::AbstractAlgebra.FPModule{elem_type(R)}

      n = ngens(M)
      R = base_ring(M)
      S = MatrixSpace(R, n, n)
      N = randmat_with_rank(S, n, -10:10)
      f = ModuleIsomorphism(M, M, N)

      @test mat(f) == N
      @test N*inverse_mat(f) == 1

      m = rand(M, -10:10)

      @test inv(f)(f(m)) == m

      @test isa(image_fn(f), Function)
      @test isa(inverse_image_fn(f), Function)
   end
end
