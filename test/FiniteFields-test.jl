function test_ffield_constructors()
   print("FiniteField.constructors()...")

   R, x = FiniteField(7, 5, "x")

   @test isa(R, FqNmodFiniteField)

   S, y = FiniteField(ZZ(17), 3, "y")

   @test isa(S, FqFiniteField)

   @test isa(3x^4 + 2x^3 + 4x^2 + x + 1, fq_nmod)

   @test isa(2y^2 + 11y + 16, fq)

   h = R(4)

   @test isa(h, fq_nmod)

   @test isa(S(11), fq)

   @test isa(R(h), fq_nmod)

   println("PASS")
end

function test_ffield_manipulation()
   print("FiniteField.manipulation()...")

   R, x = FiniteField(7, 5, "x")

   @test iszero(zero(R))

   @test isone(one(R))

   @test gen(R) == x # there is no isgen in flint

   @test characteristic(R) == 7

   @test order(R) == ZZ(7)^5

   @test degree(R) == 5

   @test isunit(x + 1)

   println("PASS")
end

function test_ffield_unary_ops()
   print("FiniteField.unary_ops()...")

   R, x = FiniteField(7, 5, "x")

   a = x^4 + 3x^2 + 6x + 1
   
   @test -a == 6*x^4+4*x^2+x+6

   println("PASS")
end

function test_ffield_binary_ops()
   print("FiniteField.binary_ops()...")

   R, x = FiniteField(7, 5, "x")

   a = x^4 + 3x^2 + 6x + 1
   b = 3x^4 + 2x^2 + x + 1

   @test a + b == 4*x^4+5*x^2+2

   @test a - b == 5*x^4+x^2+5*x

   @test a*b == 3*x^3+2

   @test gcd(a, b) == 1

   @test gcd(R(0), R(0)) == 0

   println("PASS")
end

function test_ffield_adhoc_binary()
   print("FiniteField.adhoc_binary()...")

   R, x = FiniteField(7, 5, "x")

   a = x^4 + 3x^2 + 6x + 1

   @test 3a == 3*x^4+2*x^2+4*x+3

   @test a*ZZ(5) == 5*x^4+x^2+2*x+5

   println("PASS")
end

function test_ffield_powering()
   print("FiniteField.powering()...")

   R, x = FiniteField(7, 5, "x")

   a = x^4 + 3x^2 + 6x + 1

   @test a^3 == x^4+6*x^3+5*x^2+5*x+6
   
   @test a^ZZ(-5) == x^4+4*x^3+6*x^2+6*x+2

   println("PASS")
end

function test_ffield_comparison()
   print("FiniteField.comparison()...")
  
   R, x = FiniteField(7, 5, "x")

   a = x^4 + 3x^2 + 6x + 1
   b = 3x^4 + 2x^2 + 2

   @test b != a
   @test R(3) == R(3)

   println("PASS")
end

function test_ffield_exact_division()
   print("FiniteField.exact_division()...")
  
   R, x = FiniteField(7, 5, "x")

   a = x^4 + 3x^2 + 6x + 1
   b = 3x^4 + 2x^2 + 2

   @test divexact(a, b) == 3*x^4+2*x^3+2*x^2+5*x

   @test b//a == 4*x^2+6*x+5

   println("PASS")
end

function test_ffield_inversion()
   print("FiniteField.inversion()...")
  
   R, x = FiniteField(7, 5, "x")

   a = x^4 + 3x^2 + 6x + 1

   b = inv(a)

   @test b == x^4+5*x^3+4*x^2+5*x

   @test b == a^-1

   println("PASS")
end

function test_ffield_special_functions()
   print("FiniteField.special_functions()...")
  
   R, x = FiniteField(7, 5, "x")

   a = x^4 + 3x^2 + 6x + 1

   @test trace(a) == 1

   @test norm(a) == 4

   @test frobenius(a) == x^4+2*x^3+3*x^2+5*x+1

   @test frobenius(a, 3) == 3*x^4+3*x^3+3*x^2+x+4

   @test pth_root(a) == 4*x^4+3*x^3+4*x^2+5*x+2

   println("PASS")
end

function test_ffield()
   test_ffield_constructors()
   test_ffield_manipulation()
   test_ffield_unary_ops()
   test_ffield_binary_ops()
   test_ffield_adhoc_binary()
   test_ffield_powering()
   test_ffield_comparison()
   test_ffield_exact_division()
   test_ffield_inversion()
   test_ffield_special_functions()

   println("")
end
