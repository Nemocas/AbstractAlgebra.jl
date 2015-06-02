function test_fraction_constructors()
   print("Fraction.constructors...")

   R = FractionField(ZZ)

   @test isa(R, RationalField)

   @test isa(R(ZZ(3), ZZ(2)), fmpq)

   S, x = PolynomialRing(ZZ, "x")
   T = FractionField(S)

   @test isa(T(3), Fraction)

   @test isa(T(ZZ(7)), Fraction)

   @test isa(T(x + 2), Fraction)

   @test isa(T(T(x + 2)), Fraction)

   @test isa(T(), Fraction)

   @test isa(ZZ(5)//ZZ(7), fmpq)

   @test isa((x + 3)//(x^2 + 2), Fraction)

   println("PASS")
end

function test_fraction_manipulation()
   print("Fraction.manipulation...")

   R = FractionField(ZZ)
   S, x = PolynomialRing(ZZ, "x")

   @test num(ZZ(4)//ZZ(6)) == 2

   @test den((x + 1)//(-x^2 + 1)) == x - 1

   @test iszero(zero(R))

   @test isone(one(S))

   @test canonical_unit((x + 1)//(-x^2 + 1)) == -1//(x-1)

   @test isunit((x + 1)//(-x^2 + 1))

   println("PASS")
end

function test_fraction_unary_ops()
   print("Fraction.unary_ops...")

   S, x = PolynomialRing(ZZ, "x")

   @test -(ZZ(4)//ZZ(6)) == ZZ(-2)//ZZ(3)

   @test -((x + 1)//(-x^2 + 1)) == 1//(x - 1)

   println("PASS")
end

function test_fraction_comparison()
   print("Fraction.comparison...")

   S, x = PolynomialRing(ZZ, "x")

   a = -(ZZ(4)//ZZ(6))
   b = -((x + 1)//(-x^2 + 1))

   @test a == -ZZ(2)//ZZ(3)
   
   @test b == 1//(x - 1)

   println("PASS")
end

function test_fraction_adhoc_comparison()
   print("Fraction.adhoc_comparison...")

   S, x = PolynomialRing(ZZ, "x")

   a = 1//(x - 1)

   @test 1//a == x - 1

   @test ZZ(4)//ZZ(2) == 2

   @test ZZ(3) == ZZ(6)//ZZ(2)

   @test 1 == one(S)

   println("PASS")
end

function test_fraction_binary_ops()
   print("Fraction.binary_ops...")

   S, x = PolynomialRing(ZZ, "x")

   a = -ZZ(2)//ZZ(3) + ZZ(3)//ZZ(7)
   b = (x + 1)//(-x^2 + 1) - x//(2x + 1)
   c = (ZZ(5)//ZZ(3))*(ZZ(7)//ZZ(2))

   @test a == -ZZ(5)//21

   @test b == (-x^2-x-1)//(2*x^2-x-1)

   @test c == ZZ(35)//6

   @test a//c == -ZZ(2)//49

   @test gcd(b, (x + 1)//(x - 1)) == -1//(2*x^2-x-1)

   println("PASS")
end

function test_fraction_adhoc_binary()
   print("Fraction.adhoc_binary...")

   S, x = PolynomialRing(ZZ, "x")

   a = -ZZ(2)//ZZ(3)
   b = (x + 1)//(-x^2 + 1)

   @test a + 2 == ZZ(4)//3

   @test ZZ(3) - a == ZZ(11)//3

   @test b*(x + 1) == (-x-1)//(x-1)

   @test a//5 == -ZZ(2)//15

   @test (x + 1)//b == -x^2+1

   println("PASS")
end

function test_fraction_powering()
   print("Fraction.powering()...")

   S, x = PolynomialRing(ZZ, "x")

   a = -ZZ(2)//ZZ(3)
   b = (x + 1)//(-x^2 + 1)

   @test a^12 == ZZ(4096)//531441

   @test b^-12 == x^12-12*x^11+66*x^10-220*x^9+495*x^8-792*x^7+924*x^6-792*x^5+495*x^4-220*x^3+66*x^2-12*x+1

   println("PASS")
end

function test_fraction_inversion()
   print("Fraction.inversion()...")

   S, x = PolynomialRing(ZZ, "x")

   a = -ZZ(2)//ZZ(3)
   b = (x + 1)//(-x^2 + 1)

   @test inv(a) == -ZZ(3)//2

   @test inv(b) == -x + 1

   println("PASS")
end

function test_qq_constructors()
   print("Fraction.QQ.constructors()...")

   @test isa(QQ, RationalField)

   @test isa(QQ(2), fmpq)

   println("PASS")
end

function test_qq_manipulation()
   print("Fraction.QQ.manipulation()...")

   a = -ZZ(2)//ZZ(3)
   b = ZZ(123)//ZZ(234)

   @test height(a) == 3

   @test height_bits(b) == 7

   println("PASS")
end

function test_qq_comparison()
   print("Fraction.QQ.comparison()...")

   a = -ZZ(2)//ZZ(3)
   b = ZZ(1)//ZZ(2)

   @test a < b

   @test b >= a

   @test cmp(a, b) < 0

   println("PASS")
end

function test_qq_adhoc_comparison()
   print("Fraction.QQ.adhoc_comparison()...")

   a = -ZZ(2)//ZZ(3)
   b = ZZ(1)//ZZ(2)

   @test a < 1

   @test a <= 0

   @test b >= ZZ(0)
   
   println("PASS")
end

function test_qq_adhoc_binary()
   print("Fraction.QQ.adhoc_binary()...")

   a = -ZZ(2)//ZZ(3)
   b = ZZ(1)//ZZ(2)

   @test a << 3 == -ZZ(16)//3

   @test b >> 5 == ZZ(1)//64
   
   println("PASS")
end

function test_qq_modular_arithmetic()
   print("Fraction.QQ.modular_arithmetic()...")

   a = -ZZ(2)//ZZ(3)
   b = ZZ(1)//ZZ(2)

   @test mod(a, 7) == 4

   @test mod(b, ZZ(5)) == 3
   
   println("PASS")
end

function test_qq_rational_reconstruction()
   print("Fraction.QQ.rational_reconstruction()...")

   @test reconstruct(7, 13) == ZZ(1)//2
   
   @test reconstruct(ZZ(15), 31) == -ZZ(1)//2
   
   @test reconstruct(ZZ(123), ZZ(237)) == ZZ(9)//2
   
   println("PASS")
end

function test_qq_rational_enumeration()
   print("Fraction.QQ.rational_enumeration()...")

   @test next_minimal(ZZ(2)//ZZ(3)) == ZZ(3)//2

   @test next_signed_minimal(-ZZ(21)//ZZ(31)) == ZZ(31)//21

   @test next_calkin_wilf(ZZ(321)//ZZ(113)) == ZZ(113)//244

   @test next_signed_calkin_wilf(-ZZ(51)//ZZ(17)) == ZZ(1)//4
   
   println("PASS")
end

function test_qq_special_functions()
   print("Fraction.QQ.special_functions()...")

   @test harmonic(12) == ZZ(86021)//27720

   @test dedekind_sum(12, 13) == -ZZ(11)//13

   @test dedekind_sum(-120, ZZ(1305)) == -ZZ(575)//522
   
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

   test_qq_constructors()
   test_qq_manipulation()
   test_qq_comparison()
   test_qq_adhoc_comparison()
   test_qq_adhoc_binary()
   test_qq_modular_arithmetic()
   test_qq_rational_reconstruction()
   test_qq_rational_enumeration()
   test_qq_special_functions()

   println("")
end
