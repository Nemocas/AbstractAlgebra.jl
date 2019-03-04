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
   M = FreeModule(R, 3)

   m = M([x, x, x])
   
   @test -m == M([-x, -x, -x])

   println("PASS")
end

function test_free_module_binary_ops()
   print("Generic.FreeModule.binary_ops...")

   R, x = PolynomialRing(ZZ, "x")
   M = FreeModule(R, 3)

   m = M([x, x, x])

   @test m + m - m == m

   println("PASS")
end

function test_free_module_adhoc_binary()
   print("Generic.FreeModule.adhoc_binary...")

   R, x = PolynomialRing(ZZ, "x")
   M = FreeModule(R, 3)

   m = M([x, x, x])

   @test 2*m == M([2x, 2x, 2x])
   @test m*2 == 2*m
   @test x*m == M([x^2, x^2, x^2])
   @test m*x == x*m

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
