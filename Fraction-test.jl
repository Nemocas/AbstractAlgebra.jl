function test_fraction_constructors()
   print("Fraction.constructors...")

   R = FractionField(ZZ)

   @test R <: Fraction

   @test isa(R(ZZ(3), ZZ(2)), Fraction)

   S, x = PolynomialRing(ZZ, "x")
   T = FractionField(S)

   @test isa(T(3), Fraction)

   @test isa(T(ZZ(7)), Fraction)

   @test isa(T(x + 2), Fraction)

   @test isa(T(T(x + 2)), Fraction)

   @test isa(T(), Fraction)

   @test isa(ZZ(5)/ZZ(7), Fraction)

   @test isa((x + 3)/(x^2 + 2), Fraction)

   println("PASS")
end

function test_fraction_manipulation()
   print("Fraction.manipulation...")

   R = FractionField(ZZ)
   S, x = PolynomialRing(ZZ, "x")

   @test num(ZZ(4)/ZZ(6)) == 2

   @test den((x + 1)/(-x^2 + 1)) == x - 1

   @test iszero(zero(R))

   @test isone(one(S))

   @test canonical_unit((x + 1)/(-x^2 + 1)) == -1/(x-1)

   @test isunit((x + 1)/(-x^2 + 1))

   println("PASS")
end

function test_fraction_unary_ops()
   print("Fraction.unary_ops...")

   S, x = PolynomialRing(ZZ, "x")

   @test -(ZZ(4)/ZZ(6)) == ZZ(-2)/ZZ(3)

   @test -((x + 1)/(-x^2 + 1)) == 1/(x - 1)

   println("PASS")
end

function test_fraction_comparison()
   print("Fraction.comparison...")

   S, x = PolynomialRing(ZZ, "x")

   a = -(ZZ(4)/ZZ(6))
   b = -((x + 1)/(-x^2 + 1))

   @test a == -ZZ(2)/ZZ(3)
   
   @test b == 1/(x - 1)

   println("PASS")
end

function test_fraction_adhoc_comparison()
   print("Fraction.adhoc_comparison...")

   S, x = PolynomialRing(ZZ, "x")

   a = 1/(x - 1)

   @test 1/a == x - 1

   @test ZZ(4)/ZZ(2) == 2

   @test ZZ(3) == ZZ(6)/ZZ(2)

   @test 1 == one(S)

   println("PASS")
end

function test_fraction_binary_ops()
   print("Fraction.binary_ops...")

   S, x = PolynomialRing(ZZ, "x")

   a = -ZZ(2)/ZZ(3) + ZZ(3)/ZZ(7)
   b = (x + 1)/(-x^2 + 1) - x/(2x + 1)
   c = (ZZ(5)/ZZ(3))*(ZZ(7)/ZZ(2))

   @test a == -ZZ(5)/21

   @test b == (-x^2-x-1)/(2*x^2-x-1)

   @test c == ZZ(35)/6

   @test a/c == -ZZ(2)/49

   @test gcd(b, (x + 1)/(x - 1)) == 1/(2*x^2-x-1)

   println("PASS")
end

function test_fraction_adhoc_binary()
   print("Fraction.adhoc_binary...")

   S, x = PolynomialRing(ZZ, "x")

   a = -ZZ(2)/ZZ(3)
   b = (x + 1)/(-x^2 + 1)

   @test a + 2 == ZZ(4)/3

   @test ZZ(3) - a == ZZ(11)/3

   @test b*(x + 1) == (-x-1)/(x-1)

   @test a/5 == -ZZ(2)/15

   @test (x + 1)/b == -x^2+1

   println("PASS")
end

function test_fraction_powering()
   print("Fraction.powering()...")

   S, x = PolynomialRing(ZZ, "x")

   a = -ZZ(2)/ZZ(3)
   b = (x + 1)/(-x^2 + 1)

   @test a^12 == ZZ(4096)/531441

   @test b^-12 == x^12-12*x^11+66*x^10-220*x^9+495*x^8-792*x^7+924*x^6-792*x^5+495*x^4-220*x^3+66*x^2-12*x+1

   println("PASS")
end

function test_fraction_inversion()
   print("Fraction.inversion()...")

   S, x = PolynomialRing(ZZ, "x")

   a = -ZZ(2)/ZZ(3)
   b = (x + 1)/(-x^2 + 1)

   @test inv(a) == -ZZ(3)/2

   @test inv(b) == -x + 1

   println("PASS")
end

function test_fraction()
   test_fraction_constructors()
   test_fraction_manipulation()
   test_fraction_unary_ops()
   test_fraction_comparison()
   test_fraction_adhoc_comparison()
   test_fraction_binary_ops()
   test_fraction_adhoc_binary()
   test_fraction_powering()
   test_fraction_inversion()

   println("")
end
