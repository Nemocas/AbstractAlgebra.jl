@testset "Julia.GFElem.conformance_tests" begin
   ConformanceTests.test_Field_interface_recursive(GF(3))
   ConformanceTests.test_Field_interface_recursive(GF(13))
   ConformanceTests.test_Field_interface_recursive(GF(big(13)))
   ConformanceTests.test_Field_interface_recursive(GF(big(10)^20 + 39))
end

@testset "Julia.GFElem.constructors" begin
   @test GF(2, cached = false) !== GF(2, cached = false)
   @test GF(2) === GF(2)
   @test_throws DomainError GF(4)
   @test_throws DomainError GF(4, check=true)
   F = GF(4, check=false)
   @test F isa AbstractAlgebra.GFField{Int64}
end

@testset "Julia.GFElem.printing" begin
   R = GF(13)
   S = GF(BigInt(13))

   @test string(R(3)) == "3"
   @test string(R()) == "0"
   @test string(S(3)) == "3"
   @test string(S()) == "0"

   R, x = polynomial_ring(GF(13), "x")

   @test string(x) == "x"
   @test !occursin("1", string(x^2+2*x))
end

@testset "Julia.GFElem.manipulation" begin
   R = GF(13)
   S = GF(BigInt(13))

   @test characteristic(R) == 13
   @test characteristic(S) == 13

   @test order(R) == 13
   @test order(S) == 13

   @test degree(R) == 1
   @test degree(S) == 1

   @test !is_unit(R())
   @test !is_unit(S())
   @test is_unit(R(3))
   @test is_unit(S(3))
   @test is_zero_divisor(R())
   @test is_zero_divisor(S())
   @test !is_zero_divisor(R(3))
   @test !is_zero_divisor(S(3))

   @test deepcopy(R(3)) == R(3)
   @test deepcopy(S(3)) == S(3)

   R1 = GF(13)
   S1 = GF(BigInt(13))

   @test R === R1
   @test S === S1

   @test data(R(3)) == 3
   @test data(S(3)) == 3
   @test lift(R(3)) == 3
   @test lift(S(3)) == 3
   @test isa(lift(R(3)), BigInt)
   @test isa(lift(S(3)), BigInt)
end

@testset "Julia.GFElem.rand" begin
   test_rand(GF(13))
end

@testset "Julia.GFElem.unary_ops" begin
   for (T, p) in [(Int8, 127),
                  (UInt8, 251),
                  (Int16, 32749),
                  (UInt16, 65521),
                  (Int32, 2147483647),
                  (UInt32, 4294967291),
                  (Int64, 9223372036854775783),
                  (UInt64, 18446744073709551557),
                  (Int128, 170141183460469231731687303715884105727),
                  (UInt128, 340282366920938463463374607431768211297),
                  (BigInt, 13)]
      R = GF(T(p))

      for iter = 1:1000
         a = rand(R)
         @test a == -(-a)
      end
   end
end

@testset "Julia.GFElem.binary_ops" begin
   for (T, p) in [(Int8, 127),
                  (UInt8, 251),
                  (Int16, 32749),
                  (UInt16, 65521),
                  (Int32, 2147483647),
                  (UInt32, 4294967291),
                  (Int64, 9223372036854775783),
                  (UInt64, 18446744073709551557),
                  (Int128, 170141183460469231731687303715884105727),
                  (UInt128, 340282366920938463463374607431768211297),
                  (BigInt, 13)]
      R = GF(T(p))

      for iter = 1:1000
         a1 = rand(R)
         a2 = rand(R)
         a3 = rand(R)

         @test a1 + a2 == a2 + a1
         @test a1 - a2 == -(a2 - a1)
         @test a1 + R() == a1
         @test a1 + (a2 + a3) == (a1 + a2) + a3
         @test a1*(a2 + a3) == a1*a2 + a1*a3
         @test a1*a2 == a2*a1
         @test a1*R(1) == a1
         @test R(1)*a1 == a1
      end
   end

   if Int == Int32
      F = GF(2147483659)
      a = F(2147483659 - 1)
      @test add!(a, a) == F(2147483657)
      @test add!(a, a, a) == F(2147483657)
   else
      F = GF(4611686018427388039)
      a = F(4611686018427388039 - 1)
      @test add!(a, a) == F(4611686018427388037)
      @test add!(a, a, a) == F(4611686018427388037)
   end
end

@testset "Julia.GFElem.adhoc_binary" begin
   for (T, p) in [(Int8, 127),
                  (UInt8, 251),
                  (Int16, 32749),
                  (UInt16, 65521),
                  (Int32, 2147483647),
                  (UInt32, 4294967291),
                  (Int64, 9223372036854775783),
                  (UInt64, 18446744073709551557),
                  (Int128, 170141183460469231731687303715884105727),
                  (UInt128, 340282366920938463463374607431768211297),
                  (BigInt, 13)]
      R = GF(T(p))
      for iter = 1:1000
         a = rand(R)

         c1 = rand(0:100)
         c2 = rand(0:100)
         d1 = rand(BigInt(0):BigInt(100))
         d2 = rand(BigInt(0):BigInt(100))

         @test a + c1 == c1 + a
         @test a + d1 == d1 + a
         @test a - c1 == -(c1 - a)
         @test a - d1 == -(d1 - a)
         @test a*c1 == c1*a
         @test a*d1 == d1*a
         @test a*c1 + a*c2 == a*(c1 + c2)
         @test a*d1 + a*d2 == a*(d1 + d2)
      end
   end
end

@testset "Julia.GFElem.powering" begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:100
      a = R(1)
      b = S(1)

      r = rand(R)
      s = rand(S)

      for n = 0:20
         @test a == r^n
         @test b == s^n

         a *= r
         b *= s
      end
   end

   @test_throws DivideError R(0)^(-1)
   @test_throws DivideError S(0)^(-1)

   for iter = 1:100

      a = R(1)
      b = S(1)

      r = rand(R)
      s = rand(S)

      rinv = r == 0 ? R(0) : inv(r)
      sinv = s == 0 ? S(0) : inv(s)

      for n = 0:20
         @test r == 0 || a == r^(-n)
         @test s == 0 || b == s^(-n)

         a *= rinv
         b *= sinv
      end
   end
end

@testset "Julia.GFElem.comparison" begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:1000
      a = rand(R)
      b = rand(S)

      @test a != a + 1
      @test b != b + 1

      c = rand(0:100)
      d = rand(BigInt(0):BigInt(100))

      @test R(c) == R(c)
      @test S(c) == S(c)
      @test R(d) == R(d)
      @test S(d) == S(d)
   end
end

@testset "Julia.GFElem.adhoc_comparison" begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:1000
      c = rand(0:100)
      d = rand(BigInt(0):BigInt(100))

      @test R(c) == c
      @test c == R(c)
      @test R(d) == d
      @test d == R(d)

      @test S(c) == c
      @test c == S(c)
      @test S(d) == d
      @test d == S(d)
   end
end

@testset "Julia.GFElem.inversion" begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:1000
      a = rand(R)
      b = rand(S)

      @test a == 0 || inv(inv(a)) == a
      @test b == 0 || inv(inv(b)) == b

      @test a == 0 || a*inv(a) == one(R)
      @test b == 0 || b*inv(b) == one(S)
   end
end

@testset "Julia.GFElem.exact_division" begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:1000
      a1 = rand(R)
      a2 = rand(R)
      b1 = rand(S)
      b2 = rand(S)

      @test a2 == 0 || divexact(a1, a2)*a2 == a1
      @test b2 == 0 || divexact(b1, b2)*b2 == b1
   end
end

@testset "Julia.GFElem.square_root" begin
   for p in  [2, 7, 19, 65537, 1073741827, ZZ(2), ZZ(7), ZZ(19), ZZ(65537), ZZ(1073741827)]
      R = GF(p)

      z = rand(R)
      if p != 2
         while is_square(z)
            z = rand(R)
         end
      end

      for i = 1:1000
          a = rand(R)

          f1, s = is_square_with_sqrt(a^2)

          @test f1 && s^2 == a^2

          @test is_square(a^2)

          @test sqrt(a^2)^2 == a^2

          if p != 2 && !iszero(a)
             @test !is_square(z*a^2)

             @test_throws ErrorException sqrt(z*a^2)
          end
       end
   end
end

@testset "Julia.GFElem.iteration" begin
   for n = [2, 3, 5, 13, 31]
      R = GF(n)
      elts = ConformanceTests.test_iterate(R)
      @test elts == R.(0:n-1)
   end
end

@testset "Julia.GFElem.unsafe_operators" begin
   n = rand(ZZ(2)^9:ZZ(2)^10-1)
   while !is_probable_prime(n)
      n = rand(ZZ(2)^9:ZZ(2)^10-1)
   end
   R, x = polynomial_ring(GF(n), "x")
   f = rand(R, 10:10)
   g = rand(R, 10:10)
   @test f*g == g*f
end

