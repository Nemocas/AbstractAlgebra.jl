function test_elem(R::AbstractAlgebra.Rationals{BigInt})
   n = big(2)^rand(1:100)
   return rand(QQ, -n:n)
end

@testset "Julia.Rationals.conformance_tests" begin
   test_Field_interface_recursive(QQ)
end

@testset "Julia.Rationals.constructors" begin
   R = qq
   S = QQ

   @test R == fraction_field(zz)
   @test S == fraction_field(ZZ)

   a = R(11)
   b = S(11)

   @test isa(R(a), Rational{Int})
   @test isa(S(b), Rational{BigInt})
end

@testset "Julia.Rationals.rand" begin
   test_rand(QQ, 1:9) do f
      @test 1 <= numerator(f) <= 9
      @test 1 <= denominator(f) <= 9
   end
end

@testset "Julia.Rationals.manipulation" begin
   R = qq
   S = QQ

   @test is_unit(R(3))
   @test is_unit(S(3))
   @test is_zero_divisor(R(0))
   @test is_zero_divisor(S(0))
   @test !is_zero_divisor(R(3, 5))
   @test !is_zero_divisor(S(3, 5))
end

@testset "Julia.Rationals.exact_division" begin
   R = qq
   S = QQ

   for iter = 1:1000
      a1 = rand(R, -100:100)
      a2 = rand(R, -100:100)
      b1 = rand(S, -100:100)
      b2 = rand(S, -100:100)
      c1 = rand(zz, -100:100)
      c2 = rand(ZZ, -100:100)

      @test a2 == 0 || divexact(a1, a2)*a2 == a1
      @test b2 == 0 || divexact(b1, b2)*b2 == b1

      if a1 != 0
         flagR, qR = divides(a2, a1)

         @test flagR
         @test qR*a1 == a2
      end

      if b1 != 0
         flagS, qS = divides(b2, b1)

         @test flagS
         @test qS*b1 == b2
      end

      @test c1 == 0 || divexact(a1, c1)*c1 == a1
      @test c2 == 0 || divexact(b1, c2)*c2 == b1

      @test a1 == 0 || divexact(c1, a1)*a1 == c1
      @test b1 == 0 || divexact(c2, b1)*b1 == c2
   end
end

@testset "Julia.Rationals.gcd" begin
   R = qq
   S = QQ

   for iter = 1:1000
      r1 = rand(R, -100:100)
      r2 = rand(R, -100:100)
      r3 = rand(R, -100:100)
      s1 = rand(S, -100:100)
      s2 = rand(S, -100:100)
      s3 = rand(S, -100:100)

      @test gcd(r1, gcd(r2, r3)) == gcd(gcd(r1, r2), r3) == gcd(r1, r2, r3) == gcd([r1, r2, r3])
      @test gcd(s1, gcd(s2, s3)) == gcd(gcd(s1, s2), s3) == gcd(s1, s2, s3) == gcd([s1, s2, s3])
   end
end

@testset "Julia.Rationals.square_root" begin
   R = qq
   S = QQ

   for iter = 1:1000
      r = rand(R, -100:100)
      s = rand(R, -100:100)

      f = r^2
      g = s^2

      @test AbstractAlgebra.sqrt(f)^2 == f
      @test AbstractAlgebra.sqrt(g)^2 == g
      @test is_square(f)
      @test is_square(g)
   
      @test is_square_with_sqrt(f) == (true, abs(r))
      @test is_square_with_sqrt(g) == (true, abs(s))
   end

   @test is_square(-R(1)) == false
   @test is_square(-S(1)) == false

   @test is_square_with_sqrt(-R(1)) == (false, 0)
   @test is_square_with_sqrt(-S(1)) == (false, 0)
   @test is_square_with_sqrt(R(3)) == (false, 0)
   @test is_square_with_sqrt(S(3)) == (false, 0)
   @test is_square_with_sqrt(R(2, 9)) == (false, 0)
   @test is_square_with_sqrt(R(1, 5)) == (false, 0)
   @test is_square_with_sqrt(R(3, 7)) == (false, 0)
   @test is_square_with_sqrt(S(2, 9)) == (false, 0)
   @test is_square_with_sqrt(S(1, 5)) == (false, 0)
   @test is_square_with_sqrt(S(3, 7)) == (false, 0)
end

@testset "Julia.Rationals.root" begin
   @test root(BigInt(1000)//27, 3) == 10//3
   @test root(-BigInt(27)//343, 3) == -3//7
   @test root(BigInt(27)//1331, 3; check=true) == 3//11
   @test root(BigInt(16)//9, 2; check=true) == 4//3

   @test_throws DomainError root(-BigInt(1000)//81, 4)
   @test_throws DomainError root(BigInt(1000)//81, -3)
   @test_throws DomainError root(BigInt(-16)//27, 2)

   @test_throws ErrorException root(BigInt(1100)//27, 3)
   @test_throws ErrorException root(BigInt(27)//1100, 3)
   @test_throws ErrorException root(BigInt(1100)//13, 3)
   @test_throws ErrorException root(-BigInt(40)//27, 3)
   @test_throws ErrorException root(-BigInt(27)//4, 3)

   @test root(1000//27, 3) == 10//3
   @test root(-27//343, 3) == -3//7
   @test root(27//1331, 3) == 3//11
   @test root(16//9, 2) == 4//3

   @test_throws DomainError root(-1000//81, 4)
   @test_throws DomainError root(1000//81, -3)
   @test_throws DomainError root(-16//27, 2)

   @test_throws ErrorException root(1100//27, 3)
   @test_throws ErrorException root(27//1100, 3)
   @test_throws ErrorException root(1100//13, 3)
   @test_throws ErrorException root(-40//27, 3)
   @test_throws ErrorException root(-27//4, 3)

   for T in [BigInt, Int]
      for i = 1:1000
         n = rand(1:20)
         a = BigInt(rand(-1000:1000))
         if iseven(n)
            a = abs(a)
         end
         b = BigInt(rand(1:1000))
         p = a^n
         q = b^n
         if T == BigInt || (ndigits(p; base=2) < ndigits(typemax(T); base=2) &&
                            ndigits(q; base=2) < ndigits(typemax(T); base=2))
            p = T(p)
            q = T(q)

            flag, r = is_power(p//q, n)

            @test flag && r == a//b
         end
      end

      @test_throws DomainError is_power(Rational{T}(5//3), -1)
      @test_throws DomainError is_power(Rational{T}(5//3), 0)
   end
end

@testset "Julia.Rationals.exp" begin
   @test AbstractAlgebra.exp(0//1) == 1
   @test_throws DomainError AbstractAlgebra.exp(1//1)
   @test_throws DomainError AbstractAlgebra.exp(rand(2:1000)//rand(1:1000))
   @test_throws DomainError AbstractAlgebra.exp(-rand(1:1000//rand(1:1000)))
end

@testset "Julia.Rationals.divrem" begin
   R = qq
   S = QQ

   for iter = 1:1000
      r = rand(R, -100:100)
      s = zero(R)
      while iszero(s)
         s = rand(R, -100:100)
      end

      @test AbstractAlgebra.divrem(r,s) == (r/s,0)
   end
end
