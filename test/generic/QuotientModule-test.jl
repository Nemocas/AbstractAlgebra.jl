function test_quotient_module_constructors()
   print("Generic.QuotientModule.constructors...")

   R = ZZ
   M = FreeModule(R, 2)

   m = M([R(1), R(2)])

   N, f = Submodule(M, [m])
   Q, g = QuotientModule(M, N)

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

   @test isa(Q, Generic.QuotientModule)

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

   R = ZZ
   M = FreeModule(R, 2)

   m = M([R(1), R(2)])

   N, f = Submodule(M, [m])
   Q, g = QuotientModule(M, N)

   m = Q([ZZ(2)])
   
   @test -m == Q([-ZZ(2)])

   println("PASS")
end

function test_quotient_module_binary_ops()
   print("Generic.QuotientModule.binary_ops...")

   R = ZZ
   M = FreeModule(R, 2)

   m = M([R(1), R(2)])

   N, f = Submodule(M, [m])
   Q, g = QuotientModule(M, N)

   m = Q([ZZ(2)])

   @test m + m - m == m

   println("PASS")
end

function test_quotient_module_adhoc_binary()
   print("Generic.QuotientModule.adhoc_binary...")

   R = ZZ
   M = FreeModule(R, 2)

   m = M([R(1), R(2)])

   N, f = Submodule(M, [m])
   Q, g = QuotientModule(M, N)

   m = Q([ZZ(2)])

   @test 2*m == Q([ZZ(4)])
   @test m*2 == 2*m

   println("PASS")
end

function test_quotient_module()
   test_quotient_module_constructors()
   test_quotient_module_manipulation()
   test_quotient_module_unary_ops()
   test_quotient_module_binary_ops()
   test_quotient_module_adhoc_binary()

   println("")
end
