@testset "Generic.QuotientModule.constructors" begin
   R = ZZ
   M = free_module(R, 2)

   m = M([R(1), R(2)])

   N, f = sub(M, [m])
   Q, g = quo(M, N)

   V = gens(M)
   for v in V
      @test parent(g(v)) == Q
   end

   @test isa(Q, Generic.QuotientModule)

   @test elem_type(Q) == Generic.QuotientModuleElem{elem_type(R)}
   @test elem_type(Generic.QuotientModule{elem_type(R)}) == Generic.QuotientModuleElem{elem_type(R)}
   @test parent_type(Generic.QuotientModuleElem{elem_type(R)}) == Generic.QuotientModule{elem_type(R)}

   @test coefficient_ring(M) == R
   @test coefficient_ring_type(M) == typeof(coefficient_ring(M))

   @test isa(Q([R(2)]), Generic.QuotientModuleElem)

   R = QQ
   M = vector_space(R, 2)
   m = M([R(1), R(3)])
   N, f = sub(M, [m])
   Q, g = quo(M, N)

   V = gens(M)
   for v in V
      @test parent(g(v)) == Q
   end

   @test isa(Q, Generic.QuotientModule)

   for R in [ZZ, QQ]
      for iter = 1:40
         M = rand_module(R, -10:10)
         ngens1 = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens1]
         N, f = sub(M, S)
         ngens2 = rand(1:5)
         T = [rand(N, -10:10) for i in 1:ngens2]
         P, g = sub(N, T)
         Q1, h1 = quo(M, P)

         U = [f(v) for v in T]
         V, k = sub(M, U)

         Q2, h2 = quo(M, V)

         @test Q2 == Q1
      end
   end

   F = free_module(ZZ, 2)
   m = F(BigInt[1, 0])
   n = F(BigInt[0, 1])
   S, f = sub(F, [m, n])
   Q, g = quo(F, S)
   m = Q([])

   @test isa(m, Generic.QuotientModuleElem)

   F = free_module(QQ, 3)
   m = F([1, 0, 0])
   n = F([0, 1, 0])
   S, f = sub(F, [m, n])
   SS, ff = sub(F, [m])
   Q, g = quo(S, SS)
   @test dim(Q) == 1
end

@testset "Generic.QuotientModule.rank" begin
    F = free_module(QQ, 3)
    S1, _ = sub(F, elem_type(F)[])
    Q1, _ = quo(F, S1)
    @test rank(Q1) == 3

    S2, _ = sub(F, gens(F))
    Q2, _ = quo(F, S2)
    @test rank(Q2) == 0

    S3, _ = sub(F, [gen(F, 1)])
    Q3, _ = quo(F, S3)
    @test rank(Q3) == 2

    S4, _ = sub(F, [gen(F, 1), gen(F, 2)])
    Q4, _ = quo(F, S4)
    @test rank(Q4) == 1

    for iter = 1:40
        n = rand(1:20)
        F = free_module(QQ, n)
        k = rand(0:n)
        gens = [rand(F, -10:10) for i in 1:k]
        S, _ = sub(F, gens)
        Q, _ = quo(F, S)
        @test rank(Q) == n - dim(S)
    end
end

@testset "Generic.QuotientModule.manipulation" begin
   R = ZZ
   M = free_module(R, 2)

   m = M([R(1), R(2)])

   N, f = sub(M, [m])
   Q, g = quo(M, N)

   @test ngens(Q) == 1
   G = gens(Q)
   for i = 1:ngens(Q)
      @test gen(Q, i) == G[i]
      @test gen(Q, i) == Q[i]
   end

   @test_throws ArgumentError gen(Q, 0)
   @test_throws ArgumentError gen(Q, ngens(Q) + 1)

   @test supermodule(Q) == M

   for iter in 1:20
     M = rand_module(QQ, -10:10)
     ngens = rand(1:5)
     S = [rand(M, -10:10) for i in 1:ngens]
     N, f = sub(M, S)
     Q, g = quo(M, N)

     @test dim(Q) == dim(M) - rank(matrix(QQ, [s[ix] for s in S, ix in 1:dim(M)]))
   end
end

@testset "Generic.QuotientModule.unary_ops" begin
   for R in [ZZ, QQ]
      for iter = 1:40
         M = rand_module(R, -10:10)
         ngens = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens]
         N, f = sub(M, S)
         Q, g = quo(M, N)

         m = rand(Q, -10:10)

         @test -(-m) == m
      end
   end
end

@testset "Generic.QuotientModule.binary_ops" begin
   for R in [ZZ, QQ]
      for iter = 1:40
         M = rand_module(R, -10:10)
         ngens = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens]
         N, f = sub(M, S)
         Q, g = quo(M, N)

         m = rand(Q, -10:10)
         n = rand(Q, -10:10)

         @test m + n - n == m
         @test m - n == m + (-n)
      end
   end
end

@testset "Generic.QuotientModule.adhoc_binary" begin
   for R in [ZZ, QQ]
      for iter = 1:40
         M = rand_module(R, -10:10)
         ngens = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens]
         N, f = sub(M, S)
         Q, g = quo(M, N)

         m = rand(Q, -10:10)
         n = rand(Q, -10:10)
         c = rand(-10:10)

         @test 2*m == m + m
         @test m*c == c*m
         @test c*(m + n) == c*m + c*n
         @test c*(m - n) == c*m - c*n
      end
   end
end

@testset "Generic.QuotientModule.canonical_projection" begin
   for R in [ZZ, QQ]
      for iter = 1:40
         M = rand_module(R, -10:10)
         ngens = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens]
         N, f = sub(M, S)
         Q, h = quo(M, N)

         I, g = image(h)

         @test I == Q

         m = rand(M, -10:10)
         n = h(m)

         @test h(preimage(h, n)) == n
      end
   end
end
