if !on_windows64

RR = ArbField(64)

function test_arb_poly_constructors()
   print("arb_poly.constructors()...")

   R, x = PolynomialRing(RR, "x")

   @test typeof(R) <: ArbPolyRing

   #@test isa(x, PolyElem)

   f = x^3 + 2x^2 + x + 1

   #@test isa(f, PolyElem)

   g = R(2)

   #@test isa(g, PolyElem)

   h = R(x^2 + 2x + 1)

   #@test isa(h, PolyElem)

   k = R([RR(1), RR(0), RR(3)])

   #@test isa(k, PolyElem)

   println("PASS")
end



function test_arb_poly()
   test_arb_poly_constructors()

   println("")
end

end # on_windows64
