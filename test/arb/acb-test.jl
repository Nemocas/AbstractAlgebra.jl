if !on_windows64

RR = ArbField(64)
CC = AcbField(64)

function test_acb_constructors()
   print("acb.constructors()...")

   @test isa(CC, AcbField)
   @test isa(CC(2), FieldElem)

   @test elem_type(CC) == acb
   @test base_ring(CC) == None

   println("PASS")
end

function test_acb_basic_ops()
   print("acb.basic_ops()...")

   @test one(CC) == 1
   @test zero(CC) == 0

   a = one(CC)
   @test CC(1) == a
   @test CC(ZZ(1)) == a
   @test CC(QQ(1)) == a
   @test CC(RR(1)) == a
   @test CC(UInt(1)) == a
   @test CC(RR(1)) == a
   @test CC("1.0") == a
   @test CC("1.0 +/- 0") == a
   @test CC("+1.00000e+0") == a

   b = CC(2,3)
   @test CC("2","3") == b
   @test CC(RR(2),RR(3)) == b
   @test real(b) == 2
   @test imag(b) == 3

   println("PASS")
end

function test_acb_comparison()
   print("acb.comparison()...")

   exact3 = CC(3)
   exact4 = CC(4)
   approx3 = CC("3 +/- 0.000001")
   approx4 = CC("4 +/- 0.000001")

   @test exact3 == exact3
   @test !(exact3 != exact3)

   @test !(exact3 == approx3)
   @test !(exact3 != approx3)

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

   @test !contains_zero(approx3)
   @test !contains_zero(-approx3)
   @test contains_zero(approx3 - 3)

   println("PASS")
end


function test_acb_predicates()
   print("acb.predicates()...")

   @test iszero(CC(0))
   @test !iszero(CC(1))
   @test !iszero(CC("0 +/- 0.01"))

   # @test !isnonzero(RR(0))
   # @test isnonzero(RR(1))
   # @test !isnonzero(RR("0 +/- 0.01"))

   @test isone(CC(1))
   @test !isone(CC(0))

   @test isfinite(CC(3))
   @test !isfinite(CC("0 +/- inf"))
   @test !isfinite(CC("nan"))

   @test isexact(CC(3))
   @test !isexact(CC("3 +/- 0.01"))
   @test isexact(CC(QQ(1,4)))
   @test !isexact(CC(QQ(1,3)))

   @test isint(CC(3))
   @test !isint(CC("3 +/- 0.01"))

   println("PASS")
end

function test_acb_unary_ops()
   print("acb.unary_ops()...")

   @test -CC(3) == CC(-3)
   @test abs(-CC(3)) == 3
   @test abs(CC(3)) == 3
   @test inv(CC(2)) == CC(QQ(1,2))

   println("PASS")
end

function test_acb_binary_ops()
   print("acb.binary_ops()...")

   x = CC(2)
   y = CC(4)

   @test x + y == 6
   @test x - y == -2
   @test x * y == 8
   @test x // y == 0.5

   @test x + UInt(4) == 6
   @test x - UInt(4) == -2
   @test x * UInt(4) == 8
   @test x // UInt(4) == 0.5
   @test UInt(2) + y == 6
   @test UInt(2) - y == -2
   @test UInt(2) * y == 8
   @test UInt(2) // y == 0.5

   @test x + Int(4) == 6
   @test x - Int(4) == -2
   @test x * Int(4) == 8
   @test x // Int(4) == 0.5
   @test Int(2) + y == 6
   @test Int(2) - y == -2
   @test Int(2) * y == 8
   @test Int(2) // y == 0.5

   @test x + ZZ(4) == 6
   @test x - ZZ(4) == -2
   @test x * ZZ(4) == 8
   @test x // ZZ(4) == 0.5
   @test ZZ(2) + y == 6
   @test ZZ(2) - y == -2
   @test ZZ(2) * y == 8
   @test ZZ(2) // y == 0.5

   @test x + QQ(4) == 6
   @test x - QQ(4) == -2
   @test x * QQ(4) == 8
   @test x // QQ(4) == 0.5
   @test QQ(2) + y == 6
   @test QQ(2) - y == -2
   @test QQ(2) * y == 8
   @test QQ(2) // y == 0.5

   @test x ^ y == 16
   @test x ^ ZZ(4) == 16
   @test x ^ UInt(4) == 16
   @test x ^ Int(4) == 16
   @test x ^ QQ(4) == 16

   @test ZZ(2) ^ y == 16
   @test UInt(2) ^ y == 16
   @test Int(2) ^ y == 16
   @test QQ(2) ^ y == 16


   println("PASS")
end

function test_acb_misc_ops()
   print("acb.misc_ops()...")

   @test ldexp(CC(3), 2) == 12
   @test ldexp(CC(3), ZZ(2)) == 12
   @test contains(trim(CC("1.1 +/- 0.001")), CC("1.1"))

   @test accuracy_bits(CC(0)) == typemax(Int)
   @test accuracy_bits(CC("+/- inf")) == -typemax(Int)
   @test accuracy_bits(CC("0.1")) > prec(CC) - 4

   uniq, n = unique_integer(CC("3 +/- 0.001"))
   @test uniq
   @test n == 3

   uniq, n = unique_integer(CC("3 +/- 1.001"))
   @test !uniq

   println("PASS")
end

function test_acb_unsafe_ops()
   print("acb.unsafe_ops()...")

   z = CC(1)
   x = CC(2)
   y = CC(3)

   add!(z, x, y)
   @test z == 5

   sub!(z, x, y)
   @test z == -1

   mul!(z, x, y)
   @test z == 6

   div!(z, y, x)
   @test z == 1.5

   println("PASS")
end

function test_acb_constants()
   print("acb.constants()...")

   @test overlaps(const_pi(CC), CC("3.141592653589793238462643 +/- 4.03e-25"))

   println("PASS")
end


function test_acb()
   test_acb_constructors()
   test_acb_basic_ops()
   test_acb_comparison()
   test_acb_predicates()
   test_acb_unary_ops()
   test_acb_binary_ops()
   test_acb_misc_ops()
   test_acb_unsafe_ops()
   test_acb_constants()

   println("")
end

end # on_windows64
