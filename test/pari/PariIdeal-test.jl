function test_PariIdeal_constructors()
   print("PariIdeal.constructors...")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = parent(Ideal(O, 2, 3))

   @test isa(S, PariIdealCollection)

   @test isa(Ideal(O, 2, 3), PariIdeal)

   @test isa(Ideal(O, x + 1, R(3)), PariIdeal)

   @test isa(Ideal(O, a + 1, K(3)), PariIdeal)

   @test isa(Ideal(O, O(a + 1), O(3)), PariIdeal)

   println("PASS")
end

function test_PariIdeal_lists()
   print("PariIdeal.lists...")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = prime_decomposition(O, 13)
   T = bounded_ideals(O, 100)

   @test length(T) == 100

   a = S[2]
   b = T[97][1]

   @test isa(a, PariIdeal)
   @test isa(b, PariIdeal)

   println("PASS")
end

function test_PariIdeal()
   test_PariIdeal_constructors()
   test_PariIdeal_lists()

   println("")
end
