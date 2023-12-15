@testset "Julia.Integers.conformance_tests" begin
   test_Ring_interface_recursive(ZZ)
   test_EuclideanRing_interface(ZZ)
end

@testset "Julia.Integers.manipulation" begin
   R = zz
   S = ZZ

   @test !is_unit(R(3))
   @test !is_unit(S(3))
   @test is_unit(R(-1))
   @test is_unit(S(-1))
   @test is_zero_divisor(R(0))
   @test is_zero_divisor(S(0))
   @test !is_zero_divisor(R(3))
   @test !is_zero_divisor(S(3))

   @test check_parent(S(1), S(1))
end

@testset "Julia.Integers.rand" begin
   test_rand(ZZ, 0:22) do f
      f in 0:22
   end
end

@testset "Julia.Integers.modular_arithmetic" begin
   R = zz
   S = ZZ

   for iter = 1:100
      a = R(1)
      b = S(1)
      modR = 23
      modS = BigInt(23)

      r = rand(R, 0:22)
      s = rand(S, 0:22)

      for n = 0:20
         @test r == 0 || a == powermod(r, n, modR)
         @test s == 0 || b == powermod(s, n, modS)
         @test powermod(r, n, modR) == powermod(r, BigInt(n), modR)
         @test powermod(s, n, modS) == powermod(s, BigInt(n), modS)

         a = mod(a*r, modR)
         b = mod(b*s, modS)
      end
   end

   @test powermod(123, 1, 5) == 3

   @test powermod(123, -1, 5) == 2
   
   n = -rand(2:100)
   
   @test powermod(123, n, 5) == invmod(powermod(123, -n, 5), 5)
end

@testset "Julia.Integers.exact_division" begin
   R = zz
   S = ZZ

   @test_throws ArgumentError divexact(10, 4)
   @test_throws ArgumentError divexact(big(10), big(4))
   @test_throws ArgumentError divexact(big(10), 4)

   for iter = 1:1000
      a1 = rand(R, -100:100)
      a2 = rand(R, -100:100)
      b1 = rand(S, -100:100)
      b2 = rand(S, -100:100)

      @test a2 == 0 || divexact(a1*a2, a2) == a1
      @test b2 == 0 || divexact(b1*b2, b2) == b1
      @test a1 == 0 || divexact(b1*a1, a1) == b1

      if a1 != 0
         flagR, qR = divides(a1*a2, a1)

         @test flagR
         @test qR == a2

         flagT, qT = divides(b1*a1, a1)

         @test flagT
         @test qT == b1
      end
      
      @test is_divisible_by(a1*a2, a1)

      if b1 != 0
         flagS, qS = divides(b1*b2, b1)

         @test flagS
         @test qS == b2
      end

      @test is_divisible_by(b1*b2, b1)

      @test is_divisible_by(b1*a1, a1)
   end
end

@testset "Julia.Integers.inv" begin
   @test AbstractAlgebra.inv(ZZ(1)) == 1
   @test AbstractAlgebra.inv(-ZZ(1)) == -1
   @test AbstractAlgebra.inv(zz(1)) == 1
   @test AbstractAlgebra.inv(-zz(1)) == -1

   @test_throws ArgumentError AbstractAlgebra.inv(ZZ(2))
   @test_throws DivideError AbstractAlgebra.inv(ZZ(0))
end

@testset "Julia.Integers.gcd" begin
   R = zz
   S = ZZ

   n = 23*31

   modR = n
   modS = BigInt(n)

   for iter = 1:1000
      r = rand(R, 0:n - 1)
      s = rand(S, 0:n - 1)

      gR, sR = gcdinv(r, modR)
      gS, sS = gcdinv(s, modS)

      @test gR == gcd(r, modR)
      @test gR != 1 || mod(sR*r, modR) == 1
      @test gS == gcd(s, modS)
      @test gS != 1 || mod(sS*s, modS) == 1
   end
end

@testset "Julia.Integers.crt" begin
   function testit(r, m, check=true)
      n = length(r)
      if n == 2 && rand(Bool)
         a = @inferred crt(r[1], m[1], r[2], m[2]; check=check)
         b, l = @inferred crt_with_lcm(r[1], m[1], r[2], m[2]; check=check)
      else
         a = @inferred crt(r, m; check=check)
         b, l = @inferred crt_with_lcm(r, m; check=check)
      end
      for i in 1:length(r)
         @test is_divisible_by(l, m[i])
         @test is_divisible_by(a - r[i], m[i])
         @test is_divisible_by(b - r[i], m[i])
      end
   end

   testit([ZZ(1)], [ZZ(4)])
   testit([ZZ(1)], [ZZ(0)])

   testit([ZZ(1), ZZ(2)], [ZZ(4), ZZ(5)])
   testit([ZZ(1), ZZ(3)], [ZZ(4), ZZ(6)], false)
   testit([ZZ(1), ZZ(3)], [ZZ(4), ZZ(6)], true)
   testit([ZZ(-1), ZZ(2)], [ZZ(0), ZZ(3)])
   testit([ZZ(-1), ZZ(2)], [ZZ(3), ZZ(0)])
   @test_throws Exception crt([ZZ(1), ZZ(2)], [ZZ(0), ZZ(3)])
   @test_throws Exception crt([ZZ(1), ZZ(2)], [ZZ(3), ZZ(0)])
   @test_throws Exception crt([ZZ(1), ZZ(2)], [ZZ(4), ZZ(6)])
   @test parent(crt([ZZ(1), ZZ(2)], [ZZ(4), ZZ(6)]; check=false)) == ZZ # junk but no throw

   testit([ZZ(1), ZZ(2), ZZ(3)], [ZZ(4), ZZ(5), ZZ(7)])
   testit([ZZ(1), ZZ(2), ZZ(3)], [ZZ(4), ZZ(5), ZZ(6)])
   testit([ZZ(-1), ZZ(2), ZZ(-1)], [ZZ(0), ZZ(3), ZZ(0)])
   @test_throws Exception crt([ZZ(1), ZZ(2), ZZ(2)], [ZZ(4), ZZ(5), ZZ(6)])
   @test_throws Exception crt([ZZ(-1), ZZ(2), ZZ(2)], [ZZ(0), ZZ(3), ZZ(0)])
   @test_throws Exception crt([ZZ(-1), ZZ(-1), ZZ(2)], [ZZ(0), ZZ(0), ZZ(4)])

   testit(Int64[1, 1152921504606847008], Int64[2, 1152921504606847009])
   testit(Int64[10, 100], Int64[1073741827, 536870923])
   testit(Int64[2147483692, 2147483658], Int64[2147483693, 2147483659])
end

@testset "Julia.Integers.square_root" begin
   R = zz
   S = ZZ

   for iter = 1:1000
      r = rand(R, 0:1000)
      s = rand(S, 0:1000)

      f = r^2
      g = s^2

      @test AbstractAlgebra.sqrt(f)^2 == f
      @test AbstractAlgebra.sqrt(g)^2 == g
      @test is_square(f)
      @test is_square(g)
      @test is_square_with_sqrt(f) == (true, r)
      @test is_square_with_sqrt(g) == (true, s)
   end

   @test is_square(-R(1)) == false
   @test is_square(-S(1)) == false

   @test is_square_with_sqrt(-R(1)) == (false, 0)
   @test is_square_with_sqrt(-S(1)) == (false, 0)
   @test is_square_with_sqrt(R(3)) == (false, 0)
   @test is_square_with_sqrt(S(3)) == (false, 0)

   @test_throws ErrorException AbstractAlgebra.sqrt(2)
   @test_throws ErrorException AbstractAlgebra.sqrt(S(2))
end

@testset "Julia.Integers.root" begin
   @test root(BigInt(1000), 3) == 10
   @test root(-BigInt(27), 3) == -3
   @test root(BigInt(27), 3; check=true) == 3
   @test root(BigInt(16), 2; check=true) == 4

   @test_throws DomainError root(-BigInt(1000), 4)
   @test_throws DomainError root(BigInt(1000), -3)
   @test_throws DomainError root(BigInt(-16), 2)

   @test_throws ErrorException root(BigInt(1100), 3)
   @test_throws ErrorException root(-BigInt(40), 3)

   @test iroot(BigInt(1000), 3) == 10
   @test iroot(BigInt(1100), 3) == 10
   @test iroot(-BigInt(40), 3) == -3
   @test iroot(BigInt(17), 2) == 4

   @test_throws DomainError iroot(-BigInt(1000), 4)
   @test_throws DomainError iroot(BigInt(1000), -3)
   @test_throws DomainError iroot(-BigInt(16), 2)

   @test root(1000, 3) == 10
   @test root(-27, 3) == -3
   @test root(27, 3) == 3
   @test root(16, 2) == 4

   @test_throws DomainError root(-1000, 4)
   @test_throws DomainError root(1000, -3)
   @test_throws DomainError root(-16, 2)

   @test_throws ErrorException root(1100, 3)
   @test_throws ErrorException root(-40, 3)

   @test iroot(1000, 3) == 10
   @test iroot(1100, 3) == 10
   @test iroot(-40, 3) == -3
   @test iroot(17, 2) == 4

   @test_throws DomainError iroot(-1000, 4)
   @test_throws DomainError iroot(1000, -3)
   @test_throws DomainError iroot(-16, 2)

   for T in [BigInt, Int]
      for i = 1:1000
         n = rand(1:20)
         a = BigInt(rand(-1000:1000))
         if iseven(n)
            a = abs(a)
         end
         p = a^n
         if T == BigInt || ndigits(p; base=2) < ndigits(typemax(T); base=2)
            p = T(p)

            flag, q = is_power(p, n)

            @test flag && q == a
         end
      end

      @test_throws DomainError is_power(T(5), -1)
      @test_throws DomainError is_power(T(5), 0)
   end
end

@testset "Julia.Integers.exp" begin
   @test AbstractAlgebra.exp(0) == 1
   @test_throws DomainError AbstractAlgebra.exp(1)
   @test_throws DomainError AbstractAlgebra.exp(rand(2:1000))
   @test_throws DomainError AbstractAlgebra.exp(-rand(1:1000))
end
