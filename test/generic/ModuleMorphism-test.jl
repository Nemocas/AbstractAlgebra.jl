function test_module_morphism_constructors()
   print("Generic.ModuleMorphism.constructors...")

   M = FreeModule(ZZ, 2)
 
   f = ModuleMorphism(M, M, matrix(ZZ, 2, 2, [1, 2, 3, 4]))

   @test isa(f, Generic.Map(FunctionalMap))

   m = M([ZZ(1), ZZ(2)])

   @test f(m) == M([ZZ(7), ZZ(10)])

   println("PASS")
end

function test_module_morphism()
   test_module_morphism_constructors()

   println("")
end
