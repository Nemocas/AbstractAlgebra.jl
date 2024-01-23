function test_elem(R::AbstractAlgebra.Generic.EuclideanRingResidueRing{BigInt})
   return rand(R, 0:characteristic(R))
end

function test_elem(R::AbstractAlgebra.Generic.EuclideanRingResidueRing{AbstractAlgebra.Generic.Poly{T}}) where T
   return rand(R, 0:100, -100:100)
end

@testset "Generic.EuclideanRingResidueRingElem.conformance_tests" begin
   test_Ring_interface(residue_ring(ZZ, 1)[1])   # is_gen fails on polys
   test_Ring_interface_recursive(residue_ring(ZZ, -4)[1])

   #
   R, = residue_ring(ZZ, 16453889)
   test_Ring_interface_recursive(R)

   #
   S, x = polynomial_ring(R, "x")
   T, = residue_ring(S, x^3 + 3x + 1)
   test_Ring_interface_recursive(T)

   #
   S, x = polynomial_ring(ZZ, "x")
   T, = residue_ring(S, x^2 + 1)
   #test_Ring_interface_recursive(T)   # TODO: currently fails because `inv(one(T))` fails

   @test !occursin("\n", sprint(show, T))
end

@testset "Generic.EuclideanRingResidueRingElem.constructors" begin
   B = ZZ

   R, f = Generic.residue_ring(B, 16453889)
   @test domain(f) === B
   @test codomain(f) === R

   S, f1 = quo(B, 16453889)
   U, f2 = quo(B, 16453889; cached=false)

   @test S === R
   @test U !== R

   v1 = rand(R, -100:100)
   @test f1(preimage(f1, v1)) == v1

   c1 = rand(B, -1000:1000)

   @test f1(c1) == S(c1)
   @test f2(c1) == S(c1)

   @test Generic.residue_ring(B, 16453889, cached = true)[1] === Generic.residue_ring(B, 16453889, cached = true)[1]
   @test Generic.residue_ring(B, 16453889, cached = false)[1] !== Generic.residue_ring(B, 16453889, cached = true)[1]

   @test_throws DomainError Generic.residue_ring(B, 0)

   @test elem_type(R) == Generic.EuclideanRingResidueRingElem{elem_type(B)}
   @test elem_type(Generic.EuclideanRingResidueRing{elem_type(B)}) == Generic.EuclideanRingResidueRingElem{elem_type(B)}
   @test parent_type(Generic.EuclideanRingResidueRingElem{elem_type(B)}) == Generic.EuclideanRingResidueRing{elem_type(B)}

   @test isa(R, Generic.EuclideanRingResidueRing)

   a = R(123)

   @test isa(a, Generic.EuclideanRingResidueRingElem)

   b = R(a)

   @test isa(b, Generic.EuclideanRingResidueRingElem)

   c = R(ZZ(12))

   @test isa(c, Generic.EuclideanRingResidueRingElem)

   d = R()

   @test isa(d, Generic.EuclideanRingResidueRingElem)

   S, x = polynomial_ring(R, "x")
   T, = residue_ring(S, x^3 + 3x + 1)

   @test isa(T, Generic.EuclideanRingResidueRing)

   f = T(x^4)

   @test isa(f, Generic.EuclideanRingResidueRingElem)

   g = T(f)

   @test isa(g, Generic.EuclideanRingResidueRingElem)

   # Poly modulus, invertible lc 
   S, x = polynomial_ring(ZZ, "x")
   T, = residue_ring(S, x^2 + 1)

   @test isa(T, Generic.EuclideanRingResidueRing)

   f = T(x^4)

   @test isa(f, Generic.EuclideanRingResidueRingElem)

   g = T(f)

   @test isa(g, Generic.EuclideanRingResidueRingElem)

   h = T()

   @test isa(h, Generic.EuclideanRingResidueRingElem)

   k = T(1)

   @test isa(k, Generic.EuclideanRingResidueRingElem)

   S, = Generic.residue_ring(B, 164538890)
   x = R(1)
   y = S(1)
   @test x in [x, y]
   @test x in [y, x]
   @test !(x in [y])

   @test x in keys(Dict(x => 1))
   @test !(y in keys(Dict(x => 1)))
end

@testset "Generic.EuclideanRingResidueRingElem.rand" begin
   R, = Generic.residue_ring(ZZ, 49)

   test_rand(R, 1:9) do f
      @test 1 <= f.data <= 9
   end

   # make with 3 arguments
   P, x = polynomial_ring(RealField, "x")
   R, = Generic.residue_ring(P, x^3)

   test_rand(R, -1:9, -3:3)
end

@testset "Generic.EuclideanRingResidueRingElem.manipulation" begin
   R, = Generic.residue_ring(ZZ, 16453889)

   @test modulus(R) == 16453889

   g = zero(R)

   @test iszero(g)

   @test modulus(g) == 16453889

   S, x = polynomial_ring(R, "x")
   T, = residue_ring(S, x^3 + 3x + 1)

   h = one(T)

   @test is_unit(h)

   @test is_zero_divisor(T(x + 9223338))
   @test !is_zero_divisor(T(x + 1))

   @test isone(h)

   @test data(h) == 1
   @test lift(h) == 1

   @test canonical_unit(R(11)) == R(11)

   @test canonical_unit(T(x + 1)) == T(x + 1)

   @test deepcopy(h) == h

   @test characteristic(R) == 16453889

   # Poly modulus, invertible lc
   S, x = polynomial_ring(ZZ, "x")
   T, = residue_ring(S, x^2 + 1)

   @test isone(one(T))
   @test iszero(zero(T))
   @test is_unit(T(x))
   @test canonical_unit(T(x)) == T(x)
   @test modulus(T) == x^2 + 1

   S, = residue_ring(zz, 23)

   @test lift(S(1)) == 1
   @test isa(lift(S(1)), BigInt)
end

@testset "Generic.EuclideanRingResidueRingElem.unary_ops" begin
   R, = Generic.residue_ring(ZZ, 16453889)

   @test -R(12345) == R(16441544)

   S, x = polynomial_ring(R, "x")
   T, = residue_ring(S, x^3 + 3x + 1)

   @test -T(x^5 + 1) == T(x^2+16453880*x+16453885)

   # Poly modulus, invertible lc
   S, x = polynomial_ring(ZZ, "x")
   T, = residue_ring(S, x^2 + 1)

   @test -T(x + 1) == T(-x - 1)
end

@testset "Generic.EuclideanRingResidueRingElem.binary_ops" begin
   R, = Generic.residue_ring(ZZ, 12)

   f = R(4)
   g = R(6)

   @test f + g == R(10)

   @test f - g == R(10)

   @test f*g == R(0)

   Q, = Generic.residue_ring(ZZ, 7)
   S, x = polynomial_ring(Q, "x")
   T, = residue_ring(S, x^3 + 3x + 1)

   n = T(x^5 + 1)
   p = T(x^2 + 2x + 1)

   @test n + p == T(4x + 5)

   @test n - p == T(5x^2 + 3)

   @test n*p == T(3x^2 + 4x + 4)

   # Poly modulus, invertible lc
   S, x = polynomial_ring(ZZ, "x")
   T, = residue_ring(S, x^2 + 1)

   n = T(x^5 + 1)
   p = T(x^2 + 2x + 1)

   @test n + p == T(3x + 1)

   @test n - p == T(-x + 1)

   @test n*p == T(2x - 2)
end

@testset "Generic.EuclideanRingResidueRingElem.gcd" begin
   R, = Generic.residue_ring(ZZ, 12)

   f = R(4)
   g = R(6)

   @test gcd(f, g) == R(2)

   Q, = Generic.residue_ring(ZZ, 7)
   S, x = polynomial_ring(Q, "x")
   T, = residue_ring(S, x^3 + 3x + 1)

   n = T(x^5 + 1)
   p = T(x^2 + 2x + 1)

   @test gcd(n, p) == 1

   # Poly modulus, invertible lc
   S, x = polynomial_ring(ZZ, "x")
   T, = residue_ring(S, x^2 + 2x + 1)
   
   @test gcd(T(x + 1), T(x + 1)) == T(x + 1)

   T, = residue_ring(S, x^2 + 1)

   @test gcd(T(x + 1), T(x + 1)) == T(1)
end

@testset "Generic.EuclideanRingResidueRingElem.adhoc_binary" begin
   R, = Generic.residue_ring(ZZ, 7)

   a = R(3)

   @test a + 3 == R(6)

   @test 3 - a == R(0)

   @test 5a == R(1)

   S, x = polynomial_ring(R, "x")
   T, = residue_ring(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test f + 4 == T(x^5 + 5)

   @test 4 - f == T(x^2+5*x)

   @test f*5 == T(2*x^2+3*x+6)

   # Poly modulus, invertible lc
   S, x = polynomial_ring(ZZ, "x")
   T, = residue_ring(S, x^2 + 2x + 1)

   f = T(x + 1)

   @test f + 4 == T(x + 5)

   @test 4 - f == T(-x + 3)

   @test f*5 == T(5x + 5)
end

@testset "Generic.EuclideanRingResidueRingElem.comparison" begin
   R, = Generic.residue_ring(ZZ, 7)

   a = R(3)
   b = a
   c = R(2)

   @test b == a

   @test isequal(b, a)

   @test c != a

   S, x = polynomial_ring(R, "x")
   T, = residue_ring(S, x^3 + 3x + 1)

   f = T(x^5 + 1)
   g = 8f
   h = f + g

   @test f == g
   @test h != g

   @test isequal(f, g)

   # Poly modulus, invertible lc
   S, x = polynomial_ring(ZZ, "x")
   T, = residue_ring(S, x^2 + 1)

   @test T(x + 1) == T(x + 1)
   @test isequal(T(x + 2), T(x + 2))
end

@testset "Generic.EuclideanRingResidueRingElem.adhoc_comparison" begin
   R, = Generic.residue_ring(ZZ, 7)

   a = R(3)

   @test a == 3
   @test 4 != a

   S, x = polynomial_ring(R, "x")
   T, = residue_ring(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test f != 5

   # Poly modulus, invertible lc
   S, x = polynomial_ring(ZZ, "x")
   T, = residue_ring(S, x^2 + 1)

   @test T(2) == 2
   @test T(x) != 2
end

@testset "Generic.EuclideanRingResidueRingElem.powering" begin
   R, = Generic.residue_ring(ZZ, 7)

   a = R(3)

   @test a^5 == 5

   S, x = polynomial_ring(R, "x")
   T, = residue_ring(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test f^100 == T(x^2 + 2x + 1)

   # Poly modulus, invertible lc
   S, x = polynomial_ring(ZZ, "x")
   T, = residue_ring(S, x^2 + 1)

   @test T(x + 1)^0 == T(1)
   @test T(x + 1)^1 == T(x + 1)
   @test T(x + 1)^2 == T(2x)
   @test T(x + 1)^3 == T(2x - 2)


   R, = residue_ring(ZZ, ZZ(4))
   @test_throws NotInvertibleError R(4)^-1
   @test_throws NotInvertibleError R(2)^-1
   try
      R(2)^-1
   catch e
      @test e isa NotInvertibleError
      @test e.data == 2
      @test modulus(e.mod) == 4
   end

   R, = residue_ring(ZZ, 4)
   @test_throws NotInvertibleError R(4)^-1
   @test_throws NotInvertibleError R(2)^-1
   try
      R(2)^-1
   catch e
      @test e isa NotInvertibleError
      @test e.data == 2
      @test modulus(e.mod) == 4
   end

   R, = residue_ring(ZZ, ZZ(5))
   @test_throws NotInvertibleError R(5)^-1

   R, = residue_ring(ZZ, 5)
   @test_throws NotInvertibleError R(5)^-1
end

@testset "Generic.EuclideanRingResidueRingElem.inversion" begin
   R, = Generic.residue_ring(ZZ, 49)

   a = R(5)

   @test inv(a) == 10

   R, = Generic.residue_ring(ZZ, 41)
   S, x = polynomial_ring(R, "x")
   T, = residue_ring(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test inv(f) == T(26*x^2+31*x+10)
end

@testset "Generic.EuclideanRingResidueRingElem.exact_division" begin
   R, = Generic.residue_ring(ZZ, 49)

   a = R(5)
   b = R(3)

   @test divexact(a, b) == 18

   R, = Generic.residue_ring(ZZ, 41)
   S, x = polynomial_ring(R, "x")
   T, = residue_ring(S, x^3 + 3x + 1)

   f = T(x^5 + 1)
   g = T(x^4 + x + 2)

   @test divexact(f, g) == T(7*x^2+25*x+26)

   for i = 1:100
      n = rand(1:24)
      R, = residue_ring(ZZ, n)

      for iter = 1:100
         a1 = rand(R, 0:n - 1)
         a2 = rand(R, 0:n - 1)
         a2 += Int(a2 == 0) # still works mod 1
         p = a1*a2

         q = divexact(p, a2)

         @test q*a2 == p
      end
   end
end
