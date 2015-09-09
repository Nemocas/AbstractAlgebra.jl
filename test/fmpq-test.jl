function test_fmpq_constructors()
   print("fmpq.constructors()...")

   @test isa(QQ, FlintRationalField)

   @test isa(QQ(2), FractionElem)

   println("PASS")
end

function test_fmpq_manipulation()
   print("fmpq.manipulation()...")

   a = -ZZ(2)//ZZ(3)
   b = ZZ(123)//ZZ(234)

   @test height(a) == 3

   @test height_bits(b) == 7

   println("PASS")
end

function test_fmpq_comparison()
   print("fmpq.comparison()...")

   a = -ZZ(2)//ZZ(3)
   b = ZZ(1)//ZZ(2)

   @test a < b

   @test b >= a

   @test cmp(a, b) < 0

   println("PASS")
end

function test_fmpq_adhoc_comparison()
   print("fmpq.adhoc_comparison()...")

   a = -ZZ(2)//ZZ(3)
   b = ZZ(1)//ZZ(2)

   @test a < 1

   @test a <= 0

   @test b >= ZZ(0)
   
   println("PASS")
end

function test_fmpq_adhoc_binary()
   print("fmpq.adhoc_binary()...")

   a = -ZZ(2)//ZZ(3)
   b = ZZ(1)//ZZ(2)

   @test a << 3 == -ZZ(16)//3

   @test b >> 5 == ZZ(1)//64
   
   println("PASS")
end

function test_fmpq_modular_arithmetic()
   print("fmpq.modular_arithmetic()...")

   a = -ZZ(2)//ZZ(3)
   b = ZZ(1)//ZZ(2)

   @test mod(a, 7) == 4

   @test mod(b, ZZ(5)) == 3
   
   println("PASS")
end

function test_fmpq_rational_reconstruction()
   print("fmpq.rational_reconstruction()...")

   @test reconstruct(7, 13) == ZZ(1)//2
   
   @test reconstruct(ZZ(15), 31) == -ZZ(1)//2
   
   @test reconstruct(ZZ(123), ZZ(237)) == ZZ(9)//2
   
   println("PASS")
end

function test_fmpq_rational_enumeration()
   print("fmpq.rational_enumeration()...")

   @test next_minimal(ZZ(2)//ZZ(3)) == ZZ(3)//2

   @test next_signed_minimal(-ZZ(21)//ZZ(31)) == ZZ(31)//21

   @test next_calkin_wilf(ZZ(321)//ZZ(113)) == ZZ(113)//244

   @test next_signed_calkin_wilf(-ZZ(51)//ZZ(17)) == ZZ(1)//4
   
   println("PASS")
end

function test_fmpq_special_functions()
   print("fmpq.special_functions()...")

   @test harmonic(12) == ZZ(86021)//27720

   @test dedekind_sum(12, 13) == -ZZ(11)//13

   @test dedekind_sum(-120, ZZ(1305)) == -ZZ(575)//522
   
   println("PASS")
end

function test_fmpq()
   test_fmpq_constructors()
   test_fmpq_manipulation()
   test_fmpq_comparison()
   test_fmpq_adhoc_comparison()
   test_fmpq_adhoc_binary()
   test_fmpq_modular_arithmetic()
   test_fmpq_rational_reconstruction()
   test_fmpq_rational_enumeration()
   test_fmpq_special_functions()

   println("")
end
