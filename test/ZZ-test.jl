function test_zz_constructors()
   print("ZZ.constructors...")

   a = ZZ(-123)
   @test isa(a, ZZ)

   b = ZZ(12.0)
   @test isa(b, ZZ)

   c = ZZ("-1234567876545678376545678900000000000000000000000000")
   @test isa(c, ZZ)

   d = ZZ(c)
   @test isa(d, ZZ)

   e = deepcopy(c)
   @test isa(e, ZZ)
   
   println("PASS")
end

function test_zz_manipulation()
   print("ZZ.manipulation...")

   a = one(ZZ)
   b = zero(ZZ)

   @test isa(a, ZZ)
   
   @test isa(b, ZZ)

   @test sign(a) > 0

   @test fits(Int, a)
   
   @test size(a) == 1

   @test canonical_unit(ZZ(-12)) == -1

   @test isunit(ZZ(-1))

   @test iszero(b)

   @test isone(a)
   
   println("PASS")
end

function test_zz_binary_ops()
   print("ZZ.binary_ops...")

   a = ZZ(12)
   b = ZZ(26)

   @test a + b == 38

   @test a - b == -14

   @test a*b == 312

   @test b%a == 2

   @test fdiv(b, a) == 2

   @test cdiv(b, a) == 3

   @test tdiv(b, a) == 2

   @test div(b, a) == 2

   @test mod(b, a) == 2

   @test rem(b, a) == 2

   @test divexact(ZZ(24), ZZ(12)) == 2

   @test gcd(a, b) == 2

   @test lcm(a, b) == 156
 
   @test flog(b, a) == 1

   @test clog(b, a) == 2

   println("PASS")
end

function test_zz_adhoc_binary()
   print("ZZ.adhoc_binary...")

   a = ZZ(-12)

   @test 3 + a == -9

   @test a + 3 == -9

   @test a - 3 == -15

   @test 5 - a == 17

   @test rem(a, 7) == -5

   @test a >> 3 == -2

   @test fdivpow2(a, 2) == -3

   @test a^10 == ZZ("61917364224")
   
   println("PASS")
end

function test_zz_comparison()
   print("ZZ.comparison...")

   a = ZZ(-12)
   b = ZZ(5)

   @test a < b

   @test cmpabs(a, b) > 0

   println("PASS")
end

function test_zz_adhoc_comparison()
   print("ZZ.adhoc_comparison...")

   a = ZZ(-12)
   
   @test a < 7

   @test a == -12

   println("PASS")
end

function test_zz_unary_ops()
   print("ZZ.unary_ops...")

   @test -ZZ(12) == -12

   @test ~ZZ(-5) == 4

   @test abs(ZZ(-12)) == 12

   println("PASS")
end

function test_zz_euclidean_division()
   print("ZZ.euclidean_division...")

   @test fdivrem(ZZ(12), ZZ(5)) == (ZZ(2), ZZ(2))

   @test tdivrem(ZZ(12), ZZ(5)) == (ZZ(2), ZZ(2))

   @test divrem(ZZ(12), ZZ(5)) == (ZZ(2), ZZ(2))

   println("PASS")
end

function test_zz_roots()
   print("ZZ.roots...")

   @test isqrt(ZZ(12)) == 3

   @test isqrtrem(ZZ(12)) == (3, 3)

   @test root(ZZ(1000), 3) == 10

   println("PASS")
end

function test_zz_extended_gcd()
   print("ZZ.extended_gcd...")

   @test xgcd(ZZ(12), ZZ(5)) == (1, -2, 5)

   @test gcdinv(ZZ(5), ZZ(12)) == (1, 5)

   println("PASS")
end

function test_zz_bit_twiddling()
   print("ZZ.bit_twiddling...")

   a = ZZ(12)

   @test popcount(a) == 2

   @test nextpow2(a) == 16

   combit!(a, 2)

   @test a == 8

   println("PASS")
end

function test_zz_bases()
   print("ZZ.bases...")

   a = ZZ(12)

   @test bin(a) == "1100"

   @test base(a, 13) == "c"

   @test nbits(a) == 4

   @test ndigits(a, 3) == 3

   println("PASS")
end

function test_zz_string_io()
   print("ZZ.string_io...")

   a = ZZ(12)

   @test string(a) == "12"

   println("PASS")
end

function test_zz_modular_arithmetic()
   print("ZZ.modular_arithmetic...")

   @test powmod(ZZ(12), ZZ(110), ZZ(13)) == 1

   @test powmod(ZZ(12), 110, ZZ(13)) == 1

   @test invmod(ZZ(12), ZZ(13)) == 12

   @test sqrtmod(ZZ(12), ZZ(13)) == 5

   @test crt(ZZ(5), ZZ(13), ZZ(7), ZZ(37), true) == 44

   println("PASS")
end

function test_zz_number_theoretic()
   print("ZZ.number_theoretic...")

   @test isprime(ZZ(13))

   @test fac(100) == ZZ("93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000")

   @test sigma(ZZ(128), 10) == ZZ("1181745669222511412225")

   @test eulerphi(ZZ(12480)) == 3072

   println("PASS")
end

function test_zz()
   test_zz_constructors()
   test_zz_manipulation()
   test_zz_binary_ops()
   test_zz_adhoc_binary()
   test_zz_comparison()
   test_zz_adhoc_comparison()
   test_zz_unary_ops()
   test_zz_euclidean_division()
   test_zz_roots()
   test_zz_extended_gcd()
   test_zz_bit_twiddling()
   test_zz_bases()
   test_zz_string_io()
   test_zz_modular_arithmetic()
   test_zz_number_theoretic()

   println("")
end
