function test_padics_constructors()
   print("Padics.constructors()...")  

   R = PadicField(30, 7)
   
   @test isa(R, PadicField)

   S = PadicField(30, ZZ(65537))

   @test isa(S, PadicField)

   @test isa(R(), padic)

   @test isa(R(1), padic)

   @test isa(R(ZZ(123)), padic)

   @test isa(R(ZZ(1)//7^2), padic)

   @test isa(1 + 2*7 + 4*7^2 + O(R, 7^3), padic)

   @test isa(13 + 357*ZZ(65537) + O(S, ZZ(65537)^12), padic)

   @test isa(ZZ(1)//7^2 + ZZ(2)//7 + 3 + 4*7 + O(R, 7^2), padic)

   println("PASS")
end

function test_padics_manipulation()
   print("Padics.manipulation()...")  

   R = PadicField(30, 7)

   a = 1 + 2*7 + 4*7^2 + O(R, 7^3)
   b = 7^2 + 3*7^3 + O(R, 7^5)
   c = R(2)

   @test isone(one(R))

   @test iszero(zero(R))

   @test precision(a) == 3

   @test prime(R) == 7

   @test valuation(b) == 2

   println("PASS")
end

function test_padics_unary_ops()
   print("Padics.unary_ops()...")  

   R = PadicField(30, 7)

   a = 1 + 2*7 + 4*7^2 + O(R, 7^3)
   b = R(0)

   @test -a == 6 + 4*7^1 + 2*7^2 + O(R, 7^3)
   
   @test iszero(-b)

   println("PASS")
end

function test_padics_binary_ops()
   print("Padics.binary_ops()...")  

   R = PadicField(30, 7)

   a = 1 + 2*7 + 4*7^2 + O(R, 7^3)
   b = 7^2 + 3*7^3 + O(R, 7^5)
   c = O(R, 7^3)
   d = R(2)

   @test a + b == 1 + 2*7^1 + 5*7^2 + O(R, 7^3)

   @test a - b == 1 + 2*7^1 + 3*7^2 + O(R, 7^3)

   @test a*b == 1*7^2 + 5*7^3 + 3*7^4 + O(R, 7^5)

   @test b*c == O(R, 7^5)

   @test a*d == 2 + 4*7^1 + 1*7^2 + O(R, 7^3)

   println("PASS")
end

function test_padics_adhoc_binary()
   print("Padics.adhoc_binary()...")  

   R = PadicField(30, 7)

   a = 1 + 2*7 + 4*7^2 + O(R, 7^3)
   b = 7^2 + 3*7^3 + O(R, 7^5)
   c = O(R, 7^3)
   d = R(2)

   @test a + 2 == 3 + 2*7^1 + 4*7^2 + O(R, 7^3)

   @test 3 - b == 3 + 6*7^2 + 3*7^3 + 6*7^4 + O(R, 7^5)

   @test a*ZZ(5) == 5 + 3*7^1 + O(R, 7^3)

   @test ZZ(3)*c == O(R, 7^3)

   @test 2*d == 4

   @test 2 + d == 4

   @test iszero(d - ZZ(2))

   @test a + ZZ(1)//7^2 == ZZ(1)//7^2 + 1 + 2*7^1 + 4*7^2 + O(R, 7^3)

   @test (ZZ(12)//11)*b == 3*7^2 + 3*7^3 + O(R, 7^5)

   @test c*(ZZ(1)//7) == O(R, 7^2)

   println("PASS")
end

function test_padics_comparison()
   print("Padics.comparison()...")  

   R = PadicField(30, 7)

   a = 1 + 2*7 + 4*7^2 + O(R, 7^3)
   b = 3*7^3 + O(R, 7^5)
   c = O(R, 7^3)
   d = R(2)

   @test a == 1 + 2*7 + O(R, 7^2)

   @test b == c

   @test c == R(0)

   @test d == R(2)

   println("PASS")
end

function test_padics_adhoc_comparison()
   print("Padics.adhoc_comparison()...")  

   R = PadicField(30, 7)

   a = 1 + O(R, 7^3)
   b = O(R, 7^5)
   c = R(2)

   @test a == 1

   @test b == ZZ(0)

   @test c == 2

   @test ZZ(2) == c

   @test a == ZZ(344)//1

   println("PASS")
end

function test_padics_powering()
   print("Padics.powering()...")  

   R = PadicField(30, 7)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = O(R, 7^5)
   c = R(2)

   @test a^5 == 1 + 5*7^1 + 6*7^2 + O(R, 7^3)

   @test b^3 == O(R, 7^5)

   @test c^7 == 2 + 4*7^1 + 2*7^2

   println("PASS")
end

function test_padics_exact_division()
   print("Padics.exact_division()...")  

   R = PadicField(30, 7)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = 2 + 3*7 + O(R, 7^5)
   c = 7^2 + 2*7^3 + O(R, 7^4)
   d = 7 + 2*7^2 + O(R, 7^5)

   @test divexact(a, b) == 4 + 1*7^1 + 2*7^2 + O(R, 7^3)

   @test divexact(c, d) == 1*7^1 + O(R, 7^3)

   @test divexact(d, R(7^3)) == ZZ(1)//7^2 + ZZ(2)//7 + O(R, 7^2)

   @test divexact(R(34), R(17)) == 2

   @test divexact(a, 2) == 4 + 1*7^2 + O(R, 7^3)

   @test divexact(b, ZZ(7)) == ZZ(2)//7 + 3 + O(R, 7^4)

   @test divexact(c, ZZ(12)//7^2) == 3*7^4 + 5*7^5 + O(R, 7^6)

   @test divexact(2, d) == ZZ(2)//7 + 3 + 6*7^2 + O(R, 7^3)

   @test divexact(R(3), 3) == 1

   @test divexact(ZZ(5)//7, R(5)) == ZZ(1)//7

   println("PASS")
end

function test_padics_inversion()
   print("Padics.inversion()...")  

   R = PadicField(30, 7)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = 2 + 3*7 + O(R, 7^5)
   c = 7^2 + 2*7^3 + O(R, 7^4)
   d = 7 + 2*7^2 + O(R, 7^5)

   @test inv(a) == 1 + 6*7^1 + 5*7^2 + O(R, 7^3)

   @test inv(b) == 4 + 4*7^1 + 3*7^2 + 1*7^3 + 1*7^4 + O(R, 7^5)

   @test inv(c) == ZZ(1)//7^2 + ZZ(5)//7 + O(R, 7^0)

   @test inv(d) == ZZ(1)//7 + 5 + 3*7^1 + 6*7^2 + O(R, 7^3)

   @test inv(R(1)) == 1

   println("PASS")
end

function test_padics_square_root()
   print("Padics.square_root()...")  

   R = PadicField(30, 7)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = 2 + 3*7 + O(R, 7^5)
   c = 7^2 + 2*7^3 + O(R, 7^4)

   @test sqrt(a) == 1 + 4*7^1 + 3*7^2 + O(R, 7^3)

   @test sqrt(b) == 3 + 5*7^1 + 1*7^2 + 1*7^3 + O(R, 7^5)

   @test sqrt(c) == 1*7^1 + 1*7^2 + O(R, 7^3)

   @test sqrt(R(121)) == 3 + 5*7^1 + 6*7^2 + 6*7^3 + 6*7^4 + 6*7^5 + O(R, 7^6)

   println("PASS")
end

function test_padics_special_functions()
   print("Padics.special_functions()...")  

   R = PadicField(30, 7)

   a = 1 + 7 + 2*7^2 + O(R, 7^3)
   b = 2 + 5*7 + 3*7^2 + O(R, 7^3)
   c = 3*7 + 2*7^2 + O(R, 7^5)

   @test exp(c) == 1 + 3*7^1 + 3*7^2 + 4*7^3 + 4*7^4 + O(R, 7^5)

   @test log(a) == 1*7^1 + 5*7^2 + O(R, 7^3)

   @test exp(R(0)) == 1

   @test log(R(1)) == 0

   @test teichmuller(b) == 2 + 4*7^1 + 6*7^2 + O(R, 7^3)

   println("PASS")
end

function test_padics()
   test_padics_constructors()
   test_padics_manipulation()
   test_padics_unary_ops()
   test_padics_binary_ops()
   test_padics_adhoc_binary()
   test_padics_comparison()
   test_padics_adhoc_comparison()
   test_padics_powering()
   test_padics_exact_division()
   test_padics_inversion()
   test_padics_square_root()
   test_padics_special_functions()

   println("")
end
