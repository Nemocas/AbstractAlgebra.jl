function test_elem(R::AbstractAlgebra.Integers{BigInt})
   n = big(2)^rand(1:100)
   return rand(ZZ, -n:n)
end

@testset "Julia.Integers.conformance_tests" begin
   test_Ring_interface_recursive(ZZ)
   test_EuclideanRing_interface(ZZ)
end

@testset "Julia.Integers.manipulation" begin
   R = zz
   S = ZZ

   @test !isunit(R(3))
   @test !isunit(S(3))
   @test isunit(R(-1))
   @test isunit(S(-1))

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
      
      @test isdivisible_by(a1*a2, a1)

      if b1 != 0
         flagS, qS = divides(b1*b2, b1)

         @test flagS
         @test qS == b2
      end

      @test isdivisible_by(b1*b2, b1)

      @test isdivisible_by(b1*a1, a1)
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
      @test issquare(f)
      @test issquare(g)
      @test issquare_with_sqrt(f) == (true, r)
      @test issquare_with_sqrt(g) == (true, s)
   end

   @test issquare(-R(1)) == false
   @test issquare(-S(1)) == false

   @test issquare_with_sqrt(-R(1)) == (false, 0)
   @test issquare_with_sqrt(-S(1)) == (false, 0)
   @test issquare_with_sqrt(R(3)) == (false, 0)
   @test issquare_with_sqrt(S(3)) == (false, 0)

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

            @test ispower(p, n)

            flag, q = ispower_with_root(p, n)

            @test flag && q == a
         end
      end

      @test_throws DomainError ispower(T(5), -1)
      @test_throws DomainError ispower_with_root(T(5), 0)
   end
end

@testset "Julia.Integers.exp" begin
   @test AbstractAlgebra.exp(0) == 1
   @test_throws DomainError AbstractAlgebra.exp(1)
   @test_throws DomainError AbstractAlgebra.exp(rand(2:1000))
   @test_throws DomainError AbstractAlgebra.exp(-rand(1:1000))
end
