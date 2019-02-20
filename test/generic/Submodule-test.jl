function test_submodule_constructors()
   print("Generic.Submodule.constructors...")

   R = ZZ
   M = FreeModule(R, 2)
   m = M([R(1), R(3)])
   n = M([R(2), R(-1)])
   N, f = Submodule(M, [m, n])

   @test isa(N, Generic.Submodule)

   @test elem_type(N) == Generic.submodule_elem{elem_type(R)}
   @test elem_type(Generic.Submodule{elem_type(R)}) == Generic.submodule_elem{elem_type(R)}
   @test parent_type(Generic.submodule_elem{elem_type(R)}) == Generic.Submodule{elem_type(R)}

   @test isa(N([R(2), R(7)]), Generic.submodule_elem)

   R = QQ
   M = VectorSpace(R, 2)
   m = M([R(1), R(3)])
   n = M([R(2), R(-1)])
   N, f = Subspace(M, [m, n])

   @test isa(N, Generic.Submodule)

   println("PASS")
end

function test_submodule_manipulation()
   print("Generic.Submodule.manipulation...")

   R = ZZ
   M = FreeModule(R, 2)
   m = M([R(1), R(3)])
   n = M([R(2), R(-1)])
   N, f = Submodule(M, [m, n])

   @test ngens(N) == 2
   G = gens(N)
   for i = 1:ngens(N)
      @test gen(N, i) == G[i]
   end

   @test supermodule(N) == M

   println("PASS")
end

function test_submodule_unary_ops()
   print("Generic.Submodule.unary_ops...")

   R = ZZ
   M = FreeModule(R, 2)
   m = M([R(1), R(3)])
   n = M([R(2), R(-1)])
   N, f = Submodule(M, [m, n])

   m = N([ZZ(2), ZZ(5)])
   
   @test -m == N([-ZZ(2), -ZZ(5)])

   println("PASS")
end

function test_submodule_binary_ops()
   print("Generic.Submodule.binary_ops...")

   R = ZZ
   M = FreeModule(R, 2)
   m = M([R(1), R(3)])
   n = M([R(2), R(-1)])
   N, f = Submodule(M, [m, n])

   m = N([ZZ(2), ZZ(5)])

   @test m + m - m == m

   println("PASS")
end

function test_submodule_adhoc_binary()
   print("Generic.Submodule.adhoc_binary...")

   R = ZZ
   M = FreeModule(R, 2)
   m = M([R(1), R(3)])
   n = M([R(2), R(-1)])
   N, f = Submodule(M, [m, n])

   m = N([ZZ(2), ZZ(5)])

   @test 2*m == N([ZZ(4), ZZ(10)])
   @test m*2 == 2*m

   println("PASS")
end

function test_submodule()
   test_submodule_constructors()
   test_submodule_manipulation()
   test_submodule_unary_ops()
   test_submodule_binary_ops()
   test_submodule_adhoc_binary()

   println("")
end
