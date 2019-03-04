function test_module_homomorphism_constructors()
   print("Generic.ModuleHomomorphism.constructors...")

   M = FreeModule(ZZ, 2)
 
   f = ModuleHomomorphism(M, M, matrix(ZZ, 2, 2, [1, 2, 3, 4]))

   @test isa(f, Generic.Map(FunctionalMap))

   m = M([ZZ(1), ZZ(2)])

   @test f(m) == M([ZZ(7), ZZ(10)])

   N, nmap = Submodule(M, [m])

   g = ModuleHomomorphism(N, N, matrix(ZZ, 1, 1, [3]))

   m2 = N([ZZ(4)])

   @test g(m2) == N([ZZ(12)])

   println("PASS")
end

function test_module_homomorphism()
   test_module_homomorphism_constructors()

   println("")
end
