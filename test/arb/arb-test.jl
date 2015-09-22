RR = ArbField(64)

function test_arb_constructors()
   print("arb.constructors()...")

   @test isa(RR, ArbField)
   @test isa(RR(2), FieldElem)

   println("PASS")
end

function test_arb_basic_ops()
   print("arb.basic_ops()...")

   @test one(RR) == 1
   @test zero(RR) == 0

   a = one(RR)
   @test RR(1) == a
   @test RR(ZZ(1)) == a
   @test RR(QQ(1)) == a
   @test RR(1.0) == a
   @test RR(UInt(1.0)) == a
   @test RR(RR(1)) == a
   @test RR("1.0") == a
   @test RR("1.0 +/- 0") == a
   @test RR("+1.00000e+0") == a

   @test contains(RR((RR(1), RR(0.0001))), 1)

   println("PASS")
end

function test_arb_comparison()
   print("arb.comparison()...")

   exact3 = RR(3)
   exact4 = RR(4)
   approx3 = RR("3 +/- 0.000001")
   approx4 = RR("4 +/- 0.000001")

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

   @test strongequal(approx3, approx3)
   @test !strongequal(approx3, exact3)

   @test overlaps(approx3, exact3)
   @test overlaps(exact3, approx3)
   @test overlaps(approx3, approx3)
   @test !overlaps(approx3, approx4)

   @test contains(approx3, exact3)
   @test contains(approx3, approx3)
   @test !contains(exact3, approx3)

   @test contains(approx3, QQ(3))
   @test contains(approx3, ZZ(3))
   @test contains(approx3, 3)
   @test contains(approx3, BigFloat(3))

   @test !contains_zero(approx3)
   @test !contains_negative(approx3)
   @test contains_positive(approx3)
   @test !contains_nonpositive(approx3)
   @test contains_nonnegative(approx3)

   @test !contains_zero(-approx3)
   @test contains_negative(-approx3)
   @test !contains_positive(-approx3)
   @test contains_nonpositive(-approx3)
   @test !contains_nonnegative(-approx3)

   @test contains_zero(approx3 - 3)
   @test contains_negative(approx3 - 3)
   @test contains_positive(approx3 - 3)
   @test contains_nonpositive(approx3 - 3)
   @test contains_nonnegative(approx3 - 3)

   println("PASS")
end

function test_arb_adhoc_comparison()
   print("arb.adhoc_comparison()...")

   a = RR(3)

   @test a == 3
   @test !(a == 4)
   @test a != 4
   @test !(a != 3)
   @test a >= 3
   @test a >= 2
   @test !(a >= 4)
   @test a > 2
   @test !(a > 3)
   @test a <= 3
   @test a <= 4
   @test !(a <= 2)
   @test a < 4
   @test !(a < 3)

   @test 3 == a
   @test !(4 == a)
   @test 4 != a
   @test !(3 != a)
   @test 3 <= a
   @test 2 <= a
   @test !(4 <= a)
   @test 2 < a
   @test !(3 < a)
   @test 3 >= a
   @test 4 >= a
   @test !(2 >= a)
   @test 4 > a
   @test !(3 > a)

   @test a == ZZ(3)
   @test !(a == ZZ(4))
   @test a != ZZ(4)
   @test !(a != ZZ(3))
   @test a >= ZZ(3)
   @test a >= ZZ(2)
   @test !(a >= ZZ(4))
   @test a > ZZ(2)
   @test !(a > ZZ(3))
   @test a <= ZZ(3)
   @test a <= ZZ(4)
   @test !(a <= ZZ(2))
   @test a < ZZ(4)
   @test !(a < ZZ(3))

   @test ZZ(3) == a
   @test !(ZZ(4) == a)
   @test ZZ(4) != a
   @test !(ZZ(3) != a)
   @test ZZ(3) <= a
   @test ZZ(2) <= a
   @test !(ZZ(4) <= a)
   @test ZZ(2) < a
   @test !(ZZ(3) < a)
   @test ZZ(3) >= a
   @test ZZ(4) >= a
   @test !(ZZ(2) >= a)
   @test ZZ(4) > a
   @test !(ZZ(3) > a)

   @test a == 3.0
   @test !(a == 4.0)
   @test a != 4.0
   @test !(a != 3.0)
   @test a >= 3.0
   @test a >= 2.0
   @test !(a >= 4.0)
   @test a > 2.0
   @test !(a > 3.0)
   @test a <= 3.0
   @test a <= 4.0
   @test !(a <= 2.0)
   @test a < 4.0
   @test !(a < 3.0)

   @test 3.0 == a
   @test !(4.0 == a)
   @test 4.0 != a
   @test !(3.0 != a)
   @test 3.0 <= a
   @test 2.0 <= a
   @test !(4.0 <= a)
   @test 2.0 < a
   @test !(3.0 < a)
   @test 3.0 >= a
   @test 4.0 >= a
   @test !(2.0 >= a)
   @test 4.0 > a
   @test !(3.0 > a)

   println("PASS")
end

function test_arb_predicates()
   print("arb.predicates()...")

   @test iszero(RR(0))
   @test !iszero(RR(1))
   @test !iszero(RR("0 +/- 0.01"))

   @test !isnonzero(RR(0))
   @test isnonzero(RR(1))
   @test !isnonzero(RR("0 +/- 0.01"))

   @test isone(RR(1))
   @test !isone(RR(0))

   @test isfinite(RR(3))
   @test !isfinite(RR("0 +/- inf"))
   @test !isfinite(RR("nan"))

   @test isexact(RR(3))
   @test !isexact(RR("3 +/- 0.01"))
   @test isexact(RR(QQ(1,4)))
   @test !isexact(RR(QQ(1,3)))

   @test isint(RR(3))
   @test !isint(RR("3 +/- 0.01"))

   @test ispositive(RR(3))
   @test isnonnegative(RR(3))
   @test isnegative(RR(-3))
   @test isnonpositive(RR(-3))

   @test !ispositive(RR(0))
   @test isnonnegative(RR(0))
   @test !isnegative(RR(0))
   @test isnonpositive(RR(0))

   println("PASS")
end

function test_arb_parts()
   print("arb.parts()...")

   @test midpoint(RR(3)) == 3
   @test radius(RR(3)) == 0
   @test midpoint(RR("3 +/- 0.25")) == 3
   @test radius(RR("3 +/- 0.25")) >= 0.25

   println("PASS")
end

function test_arb_unary_ops()
   print("arb.unary_ops()...")

   @test -RR(3) == RR(-3)
   @test abs(-RR(3)) == 3
   @test abs(RR(3)) == 3
   @test inv(RR(2)) == RR(0.5)

   println("PASS")
end

function test_arb_binary_ops()
   print("arb.binary_ops()...")

   x = RR(2)
   y = RR(4)

   @test x + y == 6
   @test x - y == -2
   @test x * y == 8
   @test x / y == 0.5

   @test x + UInt(4) == 6
   @test x - UInt(4) == -2
   @test x * UInt(4) == 8
   @test x / UInt(4) == 0.5
   @test UInt(2) + y == 6
   @test UInt(2) - y == -2
   @test UInt(2) * y == 8
   @test UInt(2) / y == 0.5

   @test x + Int(4) == 6
   @test x - Int(4) == -2
   @test x * Int(4) == 8
   @test x / Int(4) == 0.5
   @test Int(2) + y == 6
   @test Int(2) - y == -2
   @test Int(2) * y == 8
   @test Int(2) / y == 0.5

   @test x + ZZ(4) == 6
   @test x - ZZ(4) == -2
   @test x * ZZ(4) == 8
   @test x / ZZ(4) == 0.5
   @test ZZ(2) + y == 6
   @test ZZ(2) - y == -2
   @test ZZ(2) * y == 8
   @test ZZ(2) / y == 0.5

   @test x ^ y == 16
   @test x ^ ZZ(4) == 16
   @test x ^ UInt(4) == 16
   @test x ^ Int(4) == 16
   @test x ^ QQ(4) == 16

   println("PASS")
end

function test_arb()
   test_arb_constructors()
   test_arb_basic_ops()
   test_arb_comparison()
   test_arb_adhoc_comparison()
   test_arb_predicates()
   test_arb_parts()
   test_arb_unary_ops()
   test_arb_binary_ops()

   println("")
end
