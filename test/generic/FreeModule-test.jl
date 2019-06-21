function test_free_module_constructors()
   print("Generic.FreeModule.constructors...")

   R, x = PolynomialRing(ZZ, "x")
   M = FreeModule(R, 5)

   @test isa(M, Generic.FreeModule)

   @test elem_type(M) == Generic.free_module_elem{elem_type(R)}
   @test elem_type(Generic.FreeModule{elem_type(R)}) == Generic.free_module_elem{elem_type(R)}
   @test parent_type(Generic.free_module_elem{elem_type(R)}) == Generic.FreeModule{elem_type(R)}

   @test isa(M, Generic.FreeModule)

   @test isa(M([x, x, x, x, x]), Generic.free_module_elem)

   F = FreeModule(ZZ, 0)

   @test isa(F([]), Generic.free_module_elem)

   println("PASS")
end

function test_free_module_manipulation()
   print("Generic.FreeModule.manipulation...")

   R, x = PolynomialRing(ZZ, "x")
   M = FreeModule(R, 5)

   @test rank(M) == 5

   println("PASS")
end

function test_free_module_unary_ops()
   print("Generic.FreeModule.unary_ops...")

   R, x = PolynomialRing(ZZ, "x")

   for iter = 1:10
      M = FreeModule(R, 3)

      v = [rand(R, 0:2, -10:10) for i in 1:3]
      w = [-c for c in v]
   
      @test -M(v) == M(w)
   end

   R = QQ

   for iter = 1:10
      M = FreeModule(R, 3)

      v = [rand(R, -10:10) for i in 1:3]
      w = [-c for c in v]

      @test -M(v) == M(w)
   end

   println("PASS")
end

function test_free_module_binary_ops()
   print("Generic.FreeModule.binary_ops...")

   R, x = PolynomialRing(ZZ, "x")

   for iter = 1:10
      M = FreeModule(R, 3)

      m = rand(M, 0:2, -10:10)
      n = rand(M, 0:2, -10:10)

      @test m + n - n == m
   end

   R = QQ

   for iter = 1:10
      M = FreeModule(R, 3)

      m = rand(M, -10:10)
      n = rand(M, -10:10)

      @test m + n - n == m
   end
   
   println("PASS")
end

function test_free_module_adhoc_binary()
   print("Generic.FreeModule.adhoc_binary...")

   R, x = PolynomialRing(ZZ, "x")

   for iter = 1:10
      M = FreeModule(R, 3)
      c = rand(1:10)

      m = rand(M, 0:2, -10:10)

      v = [rand(R, 0:2, -10:10) for i in 1:3]
      w = [x*c for c in v]

      @test 2*m == m + m
      @test m*c == c*m
      @test x*M(v) == M(w)
      @test m*x == x*m
   end

   R = QQ

   for iter = 1:10
      M = FreeModule(R, 3)
      c = rand(1:10)

      m = rand(M, -10:10)

      @test 2*m == m + m
      @test m*c == c*m
   end

   println("PASS")
end

function test_free_module()
   test_free_module_constructors()
   test_free_module_manipulation()
   test_free_module_unary_ops()
   test_free_module_binary_ops()
   test_free_module_adhoc_binary()

   println("")
end
