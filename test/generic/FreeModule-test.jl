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

function test_free_module()
   test_free_module_constructors()
   test_free_module_manipulation()

   println("")
end
