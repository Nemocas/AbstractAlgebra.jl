function test_gen_res_constructors()
   print("Generic.Res.constructors...")

   B = ZZ

   R = Generic.ResidueRing(B, 16453889)

   @test elem_type(R) == Generic.Res{elem_type(B)}
   @test elem_type(Generic.ResRing{elem_type(B)}) == Generic.Res{elem_type(B)}
   @test parent_type(Generic.Res{elem_type(B)}) == Generic.ResRing{elem_type(B)}

   @test isa(R, Generic.ResRing)

   a = R(123)

   @test isa(a, Generic.Res)

   b = R(a)

   @test isa(b, Generic.Res)

   c = R(ZZ(12))

   @test isa(c, Generic.Res)

   d = R()

   @test isa(d, Generic.Res)

   S, x = PolynomialRing(R, "x")
   T = ResidueRing(S, x^3 + 3x + 1)

   @test isa(T, Generic.ResRing)

   f = T(x^4)

   @test isa(f, Generic.Res)

   g = T(f)

   @test isa(g, Generic.Res)

   S = Generic.ResidueRing(B, 164538890)
   x = R(1)
   y = S(1)
   @test x in [x, y]
   @test x in [y, x]
   @test !(x in [y])

   @test x in keys(Dict(x => 1))
   @test !(y in keys(Dict(x => 1)))

   println("PASS")
end

function test_gen_res_rand()
   print("Generic.Res.rand...")

   R = Generic.ResidueRing(ZZ, 49)
   f = rand(R, 1:9)
   @test f isa Generic.Res
   f = rand(rng, R, 1:9)
   @test f isa Generic.Res

   println("PASS")
end

function test_gen_res_manipulation()
   print("Generic.Res.manipulation...")

   R = Generic.ResidueRing(ZZ, 16453889)

   @test modulus(R) == 16453889

   g = zero(R)

   @test iszero(g)

   @test modulus(g) == 16453889

   S, x = PolynomialRing(R, "x")
   T = ResidueRing(S, x^3 + 3x + 1)

   h = one(T)

   @test isunit(h)

   @test isone(h)

   @test data(h) == 1

   @test canonical_unit(R(11)) == R(11)

   @test canonical_unit(T(x + 1)) == T(x + 1)

   @test deepcopy(h) == h

   println("PASS")
end

function test_gen_res_unary_ops()
   print("Generic.Res.unary_ops...")

   R = Generic.ResidueRing(ZZ, 16453889)

   @test -R(12345) == R(16441544)

   S, x = PolynomialRing(R, "x")
   T = ResidueRing(S, x^3 + 3x + 1)

   @test -T(x^5 + 1) == T(x^2+16453880*x+16453885)

   println("PASS")
end

function test_gen_res_binary_ops()
   print("Generic.Res.binary_ops...")

   R = Generic.ResidueRing(ZZ, 12)

   f = R(4)
   g = R(6)

   @test f + g == R(10)

   @test f - g == R(10)

   @test f*g == R(0)

   Q = Generic.ResidueRing(ZZ, 7)
   S, x = PolynomialRing(Q, "x")
   T = ResidueRing(S, x^3 + 3x + 1)

   n = T(x^5 + 1)
   p = T(x^2 + 2x + 1)

   @test n + p == T(4x + 5)

   @test n - p == T(5x^2 + 3)

   @test n*p == T(3x^2 + 4x + 4)

   println("PASS")
end

function test_gen_res_gcd()
   print("Generic.Res.gcd...")

   R = Generic.ResidueRing(ZZ, 12)

   f = R(4)
   g = R(6)

   @test gcd(f, g) == R(2)

   Q = Generic.ResidueRing(ZZ, 7)
   S, x = PolynomialRing(Q, "x")
   T = ResidueRing(S, x^3 + 3x + 1)

   n = T(x^5 + 1)
   p = T(x^2 + 2x + 1)

   @test gcd(n, p) == 1

   println("PASS")
end

function test_gen_res_adhoc_binary()
   print("Generic.Res.adhoc_binary...")

   R = Generic.ResidueRing(ZZ, 7)

   a = R(3)

   @test a + 3 == R(6)

   @test 3 - a == R(0)

   @test 5a == R(1)

   S, x = PolynomialRing(R, "x")
   T = ResidueRing(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test f + 4 == T(x^5 + 5)

   @test 4 - f == T(x^2+5*x)

   @test f*5 == T(2*x^2+3*x+6)

   println("PASS")
end

function test_gen_res_comparison()
   print("Generic.Res.comparison...")

   R = Generic.ResidueRing(ZZ, 7)

   a = R(3)
   b = a
   c = R(2)

   @test b == a

   @test isequal(b, a)

   @test c != a

   S, x = PolynomialRing(R, "x")
   T = ResidueRing(S, x^3 + 3x + 1)

   f = T(x^5 + 1)
   g = 8f
   h = f + g

   @test f == g
   @test h != g

   @test isequal(f, g)

   println("PASS")
end

function test_gen_res_adhoc_comparison()
   print("Generic.Res.adhoc_comparison...")

   R = Generic.ResidueRing(ZZ, 7)

   a = R(3)

   @test a == 3
   @test 4 != a

   S, x = PolynomialRing(R, "x")
   T = ResidueRing(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test f != 5

   println("PASS")
end

function test_gen_res_powering()
   print("Generic.Res.powering...")

   R = Generic.ResidueRing(ZZ, 7)

   a = R(3)

   @test a^5 == 5

   S, x = PolynomialRing(R, "x")
   T = ResidueRing(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test f^100 == T(x^2 + 2x + 1)

   println("PASS")
end

function test_gen_res_inversion()
   print("Generic.Res.inversion...")

   R = Generic.ResidueRing(ZZ, 49)

   a = R(5)

   @test inv(a) == 10

   R = Generic.ResidueRing(ZZ, 41)
   S, x = PolynomialRing(R, "x")
   T = ResidueRing(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test inv(f) == T(26*x^2+31*x+10)

   println("PASS")
end

function test_gen_res_exact_division()
   print("Generic.Res.exact_division...")

   R = Generic.ResidueRing(ZZ, 49)

   a = R(5)
   b = R(3)

   @test divexact(a, b) == 18

   R = Generic.ResidueRing(ZZ, 41)
   S, x = PolynomialRing(R, "x")
   T = ResidueRing(S, x^3 + 3x + 1)

   f = T(x^5 + 1)
   g = T(x^4 + x + 2)

   @test divexact(f, g) == T(7*x^2+25*x+26)

   for i = 1:100
      n = rand(1:24)
      R = ResidueRing(ZZ, n)

      for iter = 1:100
         a1 = rand(R, 0:n - 1)
         a2 = rand(R, 0:n - 1)
         a2 += Int(a2 == 0) # still works mod 1
         p = a1*a2

         q = divexact(p, a2)

         @test q*a2 == p
      end
   end

   println("PASS")
end

function test_gen_res()
   test_gen_res_constructors()
   test_gen_res_rand()
   test_gen_res_manipulation()
   test_gen_res_unary_ops()
   test_gen_res_binary_ops()
   test_gen_res_gcd()
   test_gen_res_adhoc_binary()
   test_gen_res_comparison()
   test_gen_res_adhoc_comparison()
   test_gen_res_powering()
   test_gen_res_inversion()
   test_gen_res_exact_division()

   println("")
end
