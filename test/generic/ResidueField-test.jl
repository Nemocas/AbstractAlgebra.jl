@testset "Generic.EuclideanRingResidueFieldElem.constructors" begin
   B = ZZ

   R, = Generic.residue_field(B, 16453889)

   S, f1 = quo(Field, B, 16453889)
   U, f2 = quo(Field, B, 16453889; cached=false)

   @test S === R
   @test U !== R

   @test domain(f1) === B
   @test codomain(f1) === S

   v1 = rand(R, -100:100)
   @test f1(preimage(f1, v1)) == v1

   c1 = rand(B, -1000:1000)

   @test f1(c1) == S(c1)
   @test f2(c1) == S(c1)

   @test Generic.residue_field(B, 16453889, cached = true)[1] === Generic.residue_field(B, 16453889, cached = true)[1]
   @test Generic.residue_field(B, 16453889, cached = true)[1] !== Generic.residue_field(B, 16453889, cached = false)[1]

   @test elem_type(R) == Generic.EuclideanRingResidueFieldElem{elem_type(B)}
   @test elem_type(Generic.EuclideanRingResidueField{elem_type(B)}) == Generic.EuclideanRingResidueFieldElem{elem_type(B)}
   @test parent_type(Generic.EuclideanRingResidueFieldElem{elem_type(B)}) == Generic.EuclideanRingResidueField{elem_type(B)}

   @test isa(R, Generic.EuclideanRingResidueField)

   a = R(123)

   @test isa(a, Generic.EuclideanRingResidueFieldElem)

   b = R(a)

   @test isa(b, Generic.EuclideanRingResidueFieldElem)

   c = R(ZZ(12))

   @test isa(c, Generic.EuclideanRingResidueFieldElem)

   d = R()

   @test isa(d, Generic.EuclideanRingResidueFieldElem)

   S, x = polynomial_ring(R, "x")
   T, = residue_field(S, x^3 + 3x + 1)

   @test isa(T, Generic.EuclideanRingResidueField)

   f = T(x^4)

   @test isa(f, Generic.EuclideanRingResidueFieldElem)

   g = T(f)

   @test isa(g, Generic.EuclideanRingResidueFieldElem)

   S, = Generic.residue_ring(B, 2)
   x = R(1)
   y = S(1)
   @test x in [x, y]
   @test x in [y, x]
   @test !(x in [y])

   @test x in keys(Dict(x => 1))
   @test !(y in keys(Dict(x => 1)))
end

@testset "Generic.EuclideanRingResidueFieldElem.printing" begin
   R, = Generic.residue_field(ZZ, 16453889)

   @test string(zero(R)) == "0"
   @test string(one(R)) == "1"

   S, x = polynomial_ring(R, "x")

   @test string(zero(S)) == "0"
   @test string(one(S)) == "1"
   @test string(x) == "x"
   @test string(5*x^5+3*x^3+2*x^2+x+1) == "5*x^5 + 3*x^3 + 2*x^2 + x + 1"
end

@testset "Generic.EuclideanRingResidueFieldElem.rand" begin
   R, = Generic.residue_field(ZZ, 16453889)

   test_rand(R, 1:9) do f
      @test 1 <= f.data <= 9
   end
end

@testset "Generic.EuclideanRingResidueFieldElem.manipulation" begin
   R, = Generic.residue_field(ZZ, 16453889)

   @test modulus(R) == 16453889

   g = zero(R)

   @test iszero(g)

   @test modulus(g) == 16453889

   S, x = polynomial_ring(R, "x")
   T, = residue_field(S, x^3 + 3x + 1)

   h = one(T)

   @test is_unit(h)
   @test !is_zero_divisor(h)

   @test isone(h)

   @test lift(h) == 1

   @test data(h) == 1

   @test canonical_unit(R(11)) == R(11)

   @test canonical_unit(T(x + 1)) == T(x + 1)

   @test deepcopy(h) == h

   S, = residue_field(zz, 23)

   @test lift(S(1)) == 1
   @test isa(lift(S(1)), BigInt)
end

@testset "Generic.EuclideanRingResidueFieldElem.unary_ops" begin
   R, = Generic.residue_field(ZZ, 16453889)

   @test -R(12345) == R(16441544)

   S, x = polynomial_ring(R, "x")
   T, = residue_field(S, x^3 + 3x + 1)

   @test -T(x^5 + 1) == T(x^2+16453880*x+16453885)
end

@testset "Generic.EuclideanRingResidueFieldElem.binary_ops" begin
   R, = Generic.residue_field(ZZ, 13)

   f = R(4)
   g = R(6)

   @test f + g == R(10)

   @test f - g == R(11)

   @test f*g == R(11)

   Q, = Generic.residue_field(ZZ, 7)
   S, x = polynomial_ring(Q, "x")
   T, = residue_field(S, x^3 + 3x + 1)

   n = T(x^5 + 1)
   p = T(x^2 + 2x + 1)

   @test n + p == T(4x + 5)

   @test n - p == T(5x^2 + 3)

   @test n*p == T(3x^2 + 4x + 4)
end

@testset "Generic.EuclideanRingResidueFieldElem.gcd" begin
   R, = Generic.residue_field(ZZ, 13)

   f = R(4)
   g = R(6)

   @test gcd(f, g) == R(1)

   Q, = Generic.residue_field(ZZ, 7)
   S, x = polynomial_ring(Q, "x")
   T, = residue_field(S, x^3 + 3x + 1)

   n = T(x^5 + 1)
   p = T(x^2 + 2x + 1)

   @test gcd(n, p) == 1
end

@testset "Generic.EuclideanRingResidueFieldElem.adhoc_binary" begin
   R, = Generic.residue_field(ZZ, 7)

   a = R(3)

   @test a + 3 == R(6)

   @test 3 - a == R(0)

   @test 5a == R(1)

   S, x = polynomial_ring(R, "x")
   T, = residue_field(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test f + 4 == T(x^5 + 5)

   @test 4 - f == T(x^2+5*x)

   @test f*5 == T(2*x^2+3*x+6)
end

@testset "Generic.EuclideanRingResidueFieldElem.comparison" begin
   R, = Generic.residue_field(ZZ, 7)

   a = R(3)
   b = a
   c = R(2)

   @test b == a

   @test isequal(b, a)

   @test c != a

   S, x = polynomial_ring(R, "x")
   T, = residue_field(S, x^3 + 3x + 1)

   f = T(x^5 + 1)
   g = 8f
   h = f + g

   @test f == g
   @test h != g

   @test isequal(f, g)
end

@testset "Generic.EuclideanRingResidueFieldElem.adhoc_comparison" begin
   R, = Generic.residue_field(ZZ, 7)

   a = R(3)

   @test a == 3
   @test 4 != a

   S, x = polynomial_ring(R, "x")
   T, = residue_field(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test f != 5
end

@testset "Generic.EuclideanRingResidueFieldElem.powering" begin
   R, = Generic.residue_field(ZZ, 7)

   a = R(3)

   @test a^5 == 5

   S, x = polynomial_ring(R, "x")
   T, = residue_field(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test f^100 == T(x^2 + 2x + 1)
end

@testset "Generic.EuclideanRingResidueFieldElem.inversion" begin
   R, = Generic.residue_field(ZZ, 47)

   a = R(5)

   @test inv(a) == 19

   R, = Generic.residue_field(ZZ, 41)
   S, x = polynomial_ring(R, "x")
   T, = residue_field(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test inv(f) == T(26*x^2+31*x+10)
end

@testset "Generic.EuclideanRingResidueFieldElem.exact_division" begin
   R, = Generic.residue_field(ZZ, 47)

   a = R(5)
   b = R(3)

   @test divexact(a, b) == 33

   R, = Generic.residue_field(ZZ, 41)
   S, x = polynomial_ring(R, "x")
   T, = residue_field(S, x^3 + 3x + 1)

   f = T(x^5 + 1)
   g = T(x^4 + x + 2)

   @test divexact(f, g) == T(7*x^2+25*x+26)
end

@testset "Generic.EuclideanRingResidueFieldElem.square_root" begin
   for p in [3, 47, 733, 13913, 168937, 3980299, 57586577]
       R, = Generic.residue_field(ZZ, p)

       for i = 1:10
          a = rand(R, 0:p - 1)^2
          b = sqrt(a)
          @test b^2 == a

          if !iszero(a)
             z = rand(R, 1:p - 1)
             while is_square(z)
                z = rand(R, 1:p - 1)
             end

             @test !is_square(z*a)
             @test_throws ErrorException sqrt(z*a)
         end
      end
   end

   for p in [ZZ(3), ZZ(53), ZZ(727), ZZ(8893), ZZ(191339), ZZ(2369093), ZZ(52694921)]
       R, = Generic.residue_field(ZZ, p)

       for i = 1:10
          a = rand(R, 0:Int(p - 1))^2
          b = sqrt(a)
          @test b^2 == a

          if !iszero(a)
             z = rand(R, 1:Int(p - 1))
             while is_square(z)
                z = rand(R, 1:Int(p - 1))
             end

             @test !is_square(z*a)
             @test_throws ErrorException sqrt(z*a)
         end
      end
   end
end
