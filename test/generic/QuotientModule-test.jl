function test_quotient_module_constructors()
   print("Generic.QuotientModule.constructors...")

   R = ZZ
   M = FreeModule(R, 2)

   m = M([R(1), R(2)])

   N, f = Submodule(M, [m])
   Q, g = QuotientModule(M, N)

   V = gens(M)
   for v in V
      @test parent(g(v)) == Q
   end

   @test isa(Q, Generic.QuotientModule)

   @test elem_type(Q) == Generic.quotient_module_elem{elem_type(R)}
   @test elem_type(Generic.QuotientModule{elem_type(R)}) == Generic.quotient_module_elem{elem_type(R)}
   @test parent_type(Generic.quotient_module_elem{elem_type(R)}) == Generic.QuotientModule{elem_type(R)}

   @test isa(Q([R(2)]), Generic.quotient_module_elem)

   R = QQ
   M = VectorSpace(R, 2)
   m = M([R(1), R(3)])
   N, f = Subspace(M, [m])
   Q, g = QuotientSpace(M, N)

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
         N, f = Submodule(M, S)
         ngens2 = rand(1:5)
         T = [rand(N, -10:10) for i in 1:ngens2]
         P, g = Submodule(N, T)
         Q1, h1 = QuotientModule(M, P)
         
         U = [f(v) for v in T]
         V, k = Submodule(M, U)

         Q2, h2 = QuotientModule(M, V)

         @test Q2 == Q1
      end
   end

   F = FreeModule(ZZ, 2)
   m = F(BigInt[1, 0])
   n = F(BigInt[0, 1])
   S, f = Submodule(F, [m, n])
   Q, g = QuotientModule(F, S)
   m = Q([])

   @test isa(m, Generic.quotient_module_elem)

   println("PASS")
end

function test_quotient_module_manipulation()
   print("Generic.QuotientModule.manipulation...")

   R = ZZ
   M = FreeModule(R, 2)

   m = M([R(1), R(2)])

   N, f = Submodule(M, [m])
   Q, g = QuotientModule(M, N)

   @test ngens(Q) == 1
   G = gens(Q)
   for i = 1:ngens(Q)
      @test gen(Q, i) == G[i]
   end

   @test supermodule(Q) == M

   println("PASS")
end

function test_quotient_module_unary_ops()
   print("Generic.QuotientModule.unary_ops...")

   for R in [ZZ, QQ]
      for iter = 1:40
         M = rand_module(R, -10:10)
         ngens = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens]
         N, f = Submodule(M, S)
         Q, g = QuotientModule(M, N)

         m = rand(Q, -10:10)
   
         @test -(-m) == m
      end
   end

   println("PASS")
end

function test_quotient_module_binary_ops()
   print("Generic.QuotientModule.binary_ops...")

   for R in [ZZ, QQ]
      for iter = 1:40
         M = rand_module(R, -10:10)
         ngens = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens]
         N, f = Submodule(M, S)
         Q, g = QuotientModule(M, N)

         m = rand(Q, -10:10)
         n = rand(Q, -10:10)

         @test m + n - n == m
         @test m - n == m + (-n)
      end
   end

   println("PASS")
end

function test_quotient_module_adhoc_binary()
   print("Generic.QuotientModule.adhoc_binary...")

   for R in [ZZ, QQ]
      for iter = 1:40
         M = rand_module(R, -10:10)
         ngens = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens]
         N, f = Submodule(M, S)
         Q, g = QuotientModule(M, N)

         m = rand(Q, -10:10)
         n = rand(Q, -10:10)
         c = rand(-10:10)

         @test 2*m == m + m
         @test m*c == c*m
         @test c*(m + n) == c*m + c*n
         @test c*(m - n) == c*m - c*n
      end
   end

   println("PASS")
end

function test_quotient_module_canonical_projection()
   print("Generic.QuotientModule.canonical_projection...")

   for R in [ZZ, QQ]
      for iter = 1:40
         M = rand_module(R, -10:10)
         ngens = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens]
         N, f = Submodule(M, S)
         Q, h = QuotientModule(M, N)

         I, g = image(h)

         @test I == Q
      end
   end

   println("PASS")
end

function test_quotient_module()
   test_quotient_module_constructors()
   test_quotient_module_manipulation()
   test_quotient_module_unary_ops()
   test_quotient_module_binary_ops()
   test_quotient_module_adhoc_binary()
   test_quotient_module_canonical_projection()

   println("")
end
