function test_pari_maximal_order_elem_constructors()
   print("pari_maximal_order_elem.constructors...")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   @test isa(O, PariMaximalOrder)

   b = O(123)

   @test isa(b, pari_maximal_order_elem)

   c = O(b)

   @test isa(c, pari_maximal_order_elem)

   d = O(fmpz(12))

   @test isa(d, pari_maximal_order_elem)

   f = O()

   @test isa(f, pari_maximal_order_elem)

   g = O(x^2 + 2x + 3)

   @test isa(g, pari_maximal_order_elem)

   h = O(a^2 + 2a + 3)

   @test isa(h, pari_maximal_order_elem)

   println("PASS")
end

function test_pari_maximal_order_elem_manipulation()
   print("pari_maximal_order_elem.manipulation...")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   B = basis(O)
   b1 = R(B[1])
   b2 = R(B[2])
   b3 = R(B[3])

   @test b1 == 1
   @test b2 == x
   @test b3 == x^2 + 2

   println("PASS")
end

function test_pari_maximal_order_elem()
   test_pari_maximal_order_elem_constructors()
   test_pari_maximal_order_elem_manipulation()

   println("")
end
