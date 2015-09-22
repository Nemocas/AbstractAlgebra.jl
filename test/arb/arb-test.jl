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

   println("PASS")
end

function test_arb()
   test_arb_constructors()
   test_arb_basic_ops()
   test_arb_comparison()

   println("")
end
