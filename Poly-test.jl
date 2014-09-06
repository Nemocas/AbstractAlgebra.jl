function test_constructors()
   print("Poly.constructors...")
 
   R, x = PolynomialRing(ZZ, "x")

   @test R <: Poly

   @test isa(x, Poly)

   S, y = PolynomialRing(R, "y")

   @test S <: Poly

   @test isa(y, Poly)

   T, z = PolynomialRing(S, "z")

   @test T <: Poly

   @test isa(z, Poly)

   f = x^2 + y^3 + z + 1

   @test isa(f, Poly)

   g = S(2)

   @test isa(g, Poly)

   h = S(x^2 + 2x + 1)

   @test isa(h, Poly)

   j = T(x + 2)

   @test isa(j, Poly)

   k = Poly(S, [x, x + 2, x^2 + 3x + 1])

   @test isa(k, Poly)

   l = S(k)

   @test isa(l, Poly)

   println("PASS")
end

function test_manipulation()
   print("Poly.manipulation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   @test iszero(zero(R))
   
   @test isone(one(S))

   @test isgen(gen(S))
   
   @test isunit(one(S))

   f = x^2 + 2x + 1

   @test lead(f) == 1

   h = x*y^2 + (x + 1)*y + 3

   @test coeff(h, 2) == x

   @test canonical_unit(-x + 1) == -1

   @test canonical_unit(-x*y + x + 1) == -1

   println("PASS")
end

function test_binary_ops()
   print("Poly.binary_ops...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 2

   @test f + g == x^3+x^2+5*x+3

   @test f*g == x^5+2*x^4+4*x^3+8*x^2+7*x+2

   k = x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y + (x^3 + 2x + 2)
   
   @test k - l == x*y^2+(-x^3-2*x+1)

   @test k*(f + g) == (x^4+x^3+5*x^2+3*x)*y^2+(x^4+2*x^3+6*x^2+8*x+3)*y+(3*x^3+3*x^2+15*x+9)

   @test (k - l) - g == x*y^2+-2*x^3-5*x-1

   println("PASS")
end

function test_adhoc_binary()
   print("Poly.adhoc_binary...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 2

   @test f + 12 == x^2+2*x+13

   @test ZZ(3)*g == 3*x^3+9*x+6

   k = x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y + (x^3 + 2x + 2)

   @test k - 4 == x*y^2+(x+1)*y-1
   
   @test ZZ(5) - l == (-x-1)*y+(-x^3-2*x+3)

   println("PASS")
end

function test_comparison()
   print("Poly.comparison...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 2

   @test f != g

   k = x*y^2 + (x + 1)*y + 3
   l = x*y^2 + (x + 1)*y + 3

   @test k == l

   println("PASS")
end

function test_adhoc_comparison()
   print("Poly.adhoc_comparison...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1

   @test f != 1 

   @test S(1) == ZZ(1) 

   println("PASS")
end

function test_unary_ops()
   print("Poly.unary_ops...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1

   @test -f == -x^2 - 2x - 1

   h = x*y^2 + (x + 1)*y + 3

   @test -h == -x*y^2 - (x + 1)*y - 3

   println("PASS")
end

function test_truncation()
   print("Poly.truncation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test truncate(f, 2) == 2*x+1

   @test mullow(f, g, 3) == 7*x^2+5*x+1

   k = x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y + (x^3 + 2x + 2)

   @test truncate(k, 1) == 3

   @test mullow(k, l, 4) == (x^2+x)*y^3+(x^4+3*x^2+4*x+1)*y^2+(x^4+x^3+2*x^2+7*x+5)*y+(3*x^3+6*x+6)

   println("PASS")
end

function test_poly()
   test_constructors()
   test_manipulation()
   test_binary_ops()
   test_adhoc_binary()
   test_comparison()
   test_adhoc_comparison()
   test_unary_ops()
   test_truncation()
end
