function test_submodule_constructors()
   print("Generic.Submodule.constructors...")

   R = ZZ
   M = FreeModule(R, 2)
   m = M([R(1), R(3)])
   n = M([R(2), R(-1)])
   N, f = sub(M, [m, n])

   V = gens(N)
   for v in V
      @test parent(f(v)) == M
   end

   @test isa(N, Generic.Submodule)

   @test elem_type(N) == Generic.submodule_elem{elem_type(R)}
   @test elem_type(Generic.Submodule{elem_type(R)}) == Generic.submodule_elem{elem_type(R)}
   @test parent_type(Generic.submodule_elem{elem_type(R)}) == Generic.Submodule{elem_type(R)}

   @test isa(N([R(2), R(7)]), Generic.submodule_elem)

   F = FreeModule(R, 5)
   nsubs = rand(0:5)
   subs = Generic.Submodule{elem_type(R)}[sub(F, [rand(F, -10:10)])[1] for i in 1:nsubs]
   N, h = sub(F, subs)

   @test isa(N, Generic.Submodule)

   R = QQ
   M = VectorSpace(R, 2)
   m = M([R(1), R(3)])
   n = M([R(2), R(-1)])
   N, f = sub(M, [m, n])

   @test isa(N, Generic.Submodule)
   
   V = gens(N)
   for v in V
      @test parent(f(v)) == M
   end

   M = VectorSpace(R, 5)
   nsubs = rand(1:5)
   subs = [sub(M, [rand(M, -10:10)])[1] for i in 1:nsubs]
   N, h = sub(M, subs)

   @test isa(N, Generic.Submodule)

   F = FreeModule(ZZ, 2)
   S, f = sub(F, [])
   m = S([])

   @test isa(S, Generic.Submodule)
   @test isa(m, Generic.submodule_elem)

   println("PASS")
end

function test_submodule_manipulation()
   print("Generic.Submodule.manipulation...")

   R = ZZ
   M = FreeModule(R, 2)
   m = M([R(1), R(3)])
   n = M([R(2), R(-1)])
   N, f = sub(M, [m, n])

   @test ngens(N) == 2
   G = gens(N)
   for i = 1:ngens(N)
      @test gen(N, i) == G[i]
   end

   @test supermodule(N) == M

   R = QQ
   M = VectorSpace(R, 2)
   m = M([R(1), R(3)])
   n = M([R(2), R(-1)])
   N, f = sub(M, [m, n])

   @test dim(N) == 2

   println("PASS")
end

function test_submodule_unary_ops()
   print("Generic.Submodule.unary_ops...")

   for R in [ZZ, QQ]
      for iter = 1:20
         M = rand_module(R, -10:10)
         ngens = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens]
         N, f = sub(M, S)

         m = rand(N, -10:10)
   
         @test -(-m) == m
      end
   end

   println("PASS")
end

function test_submodule_binary_ops()
   print("Generic.Submodule.binary_ops...")

   for R in [ZZ, QQ]
      for iter = 1:40
         M = rand_module(R, -10:10)
         ngens = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens]
         N, f = sub(M, S)

         m = rand(N, -10:10)
         n = rand(N, -10:10)

         @test m + n - n == m
         @test m - n == m + (-n)
      end
   end

   println("PASS")
end

function test_submodule_adhoc_binary()
   print("Generic.Submodule.adhoc_binary...")

   for R in [ZZ, QQ]
      for iter = 1:40
         M = rand_module(R, -10:10)
         ngens = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens]
         N, f = sub(M, S)

         m = rand(N, -10:10)
         n = rand(N, -10:10)
         c = rand(-10:10)

         @test 2*m == m + m
         @test m*c == c*m
         @test c*(m + n) == c*m + c*n
         @test c*(m - n) == c*m - c*n
      end
   end

   println("PASS")
end

function test_submodule_canonical_injection()
   print("Generic.Submodule.canonical_injection...")

   for R in [ZZ, QQ]
      for iter = 1:40
         M = rand_module(R, -10:10)
         ngens = rand(1:5)
         S = [rand(M, -10:10) for i in 1:ngens]
         N, f = sub(M, S)

         I, g = image(f)

         @test I == N

         m = rand(N, -10:10)

         n = f(m)
         pre = preimage(f, n)

         @test pre == m
      end
   end

   println("PASS")
end

function test_submodule()
   test_submodule_constructors()
   test_submodule_manipulation()
   test_submodule_unary_ops()
   test_submodule_binary_ops()
   test_submodule_adhoc_binary()
   test_submodule_canonical_injection()

   println("")
end
