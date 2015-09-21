RR = ArbField(64)

function test_arb_constructors()
   print("arb.constructors()...")

   @test isa(RR, ArbField)
   @test isa(RR(2), FieldElem)

   println("PASS")
end

function test_arb_comparison()
   print("arb.comparison()...")

   exact3 = RR(3)
   exact4 = RR(4)
   approx3 = RR("3 +/- 0.000001")

   @test exact3 == exact3
   @test !(exact3 != exact3)
   @test exact4 > exact3
   @test exact4 >= exact3
   @test exact3 < exact4
   @test exact3 <= exact4

   @test exact4 > approx3
   @test exact4 >= approx3
   @test approx3 < exact4
   @test approx3 <= exact4

   @test !(exact3 == approx3)
   @test !(exact3 != approx3)
   @test !(exact3 > approx3)
   @test !(exact3 >= approx3)
   @test !(exact3 < approx3)
   @test !(exact3 <= approx3)

   println("PASS")
end

function test_arb()
   test_arb_constructors()
   test_arb_comparison()

   println("")
end
