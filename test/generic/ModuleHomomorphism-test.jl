import AbstractAlgebra.PrettyPrinting

@testset "Generic.ModuleHomomorphism.constructors" begin
   M = free_module(ZZ, 2)

   f = ModuleHomomorphism(M, M, matrix(ZZ, 2, 2, [1, 2, 3, 4]))

   @test domain(f) == M
   @test codomain(f) == M

   @test isa(f, Generic.Map(FunctionalMap))

   m = M([ZZ(1), ZZ(2)])

   @test f(m) == M([ZZ(7), ZZ(10)])

   N, nmap = sub(M, [m])

   g = ModuleHomomorphism(N, N, matrix(ZZ, 1, 1, [3]))

   m2 = N([ZZ(4)])

   @test g(m2) == N([ZZ(12)])

   @test !occursin("\n", sprint(show, g))
end

@testset "Generic.ModuleHomomorphism.kernel" begin
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

         Q, g = quo(M, M)
         @test all(x -> iszero(g(x)), gens(Q))
      end
   end

   for R in [ZZ, QQ]
      for iter = 1:100
         # test random homomorphisms from free modules
         M = rand_module(R, -10:10)

         m = rand(1:5)
         F = free_module(R, m)

         S = matrix_space(R, m, ngens(M))
         f = ModuleHomomorphism(F, M, rand(S, -10:10))

         k, h = kernel(f)

         for i = 1:3
            v = rand(k, -10:10)

            @test iszero(f(h(v)))
         end
      end
   end
end

@testset "Generic.ModuleHomomorphism.image" begin
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

      Q, g = quo(M, M)
      @test all(x -> iszero(g(x)), gens(Q))
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

      Q, g = quo(N, N)
      @test all(x -> iszero(g(x)), gens(Q))
   end
end

@testset "Generic.ModuleHomomorphism.preimage" begin
   for R in [ZZ, QQ]
      for iter = 1:5
         N = rand_module(R, -10:10)

         ngens1 = rand(1:5)
         gens1 = [rand(N, -10:10) for j in 1:ngens1]
         M, f = sub(N, gens1)

         m = rand(M, -10:10)
         n = f(m)
         fl, x = has_preimage_with_preimage(f, n)
         if fl
          @test x == m
          y = f\n
          @test y == m
         end

         for i in 0:5
            m = elem_type(M)[ rand(M, -10:10) for j in 1:i ]
            n = elem_type(N)[ f(mm) for mm in m ]

            fl, x = has_preimage_with_preimage(f, n)
            if fl
              @test x == m
              y = f\n
              @test y == m
            end
         end
      end
   end
end

@testset "Generic.ModuleHomomorphism.printing" begin
  M = free_module(ZZ, 2)
  f = ModuleHomomorphism(M, M, matrix(ZZ, 2, 2, [1, 2, 3, 4]))
  str = """
        Module homomorphism
          from free module of rank 2 over integers
          to free module of rank 2 over integers"""
  @test PrettyPrinting.repr_detailed(f) == str
  @test PrettyPrinting.repr_oneline(f) == "Hom: free module of rank 2 over integers -> free module of rank 2 over integers"
  @test PrettyPrinting.repr_terse(f) == "Module homomorphism"
end

@testset "Generic.ModuleIsomorphism" begin
  let R = QQ
    for iter in 1:100
      # test image of composition of canonical injection and projection
      M = rand_module(R, -10:10)::AbstractAlgebra.FPModule{elem_type(R)}

      n = ngens(M)
      R = base_ring(M)
      S = matrix_space(R, n, n)
      N = randmat_with_rank(S, n, -10:10)
      f = ModuleIsomorphism(M, M, N)

      @test domain(f) == M
      @test codomain(f) == M
      @test matrix(f) == N
      @test N*inverse_mat(f) == 1

      m = rand(M, -10:10)

      @test inv(f)(f(m)) == m

      @test isa(image_fn(f), Function)
      @test isa(inverse_image_fn(f), Function)
    end
  end
  let R = ZZ # for ZZ not every random matrix gives an isomorphism, thus only calculate some fixed example
    F = free_module(R, 1)
    M = quo(F, sub(F, [F([R(304)])])[1])[1]

    n = ngens(M)
    S = matrix_space(R, n, n)
    for iter in 1:10
      N = matrix(R, 1, 1, [rand(R, -20:20)])
      gcd(N[1,1], 304) == 1 || continue
      f = ModuleIsomorphism(M, M, N)
      @test mod(matrix(f)[1,1]*matrix(inv(f))[1,1], 304) == 1
    end
  end
end

@testset "Generic.ModuleIsomorphism.printing" begin
  M = free_module(ZZ, 2)
  f = ModuleIsomorphism(M, M, matrix(ZZ, 2, 2, [1, 0, 0, 1]))
  str = """
        Module isomorphism
          from free module of rank 2 over integers
          to free module of rank 2 over integers"""
  @test PrettyPrinting.repr_detailed(f) == str
  @test PrettyPrinting.repr_oneline(f) == "Hom: free module of rank 2 over integers -> free module of rank 2 over integers"
  @test PrettyPrinting.repr_terse(f) == "Module isomorphism"
end
