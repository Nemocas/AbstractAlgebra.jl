@testset "Generic.TotFrac.constructors" begin
   S, = residue_ring(ZZ, 12)
   T = total_ring_of_fractions(S)

   @test total_ring_of_fractions(S, cached = true) === total_ring_of_fractions(S, cached = true)
   @test total_ring_of_fractions(S, cached = false) !== total_ring_of_fractions(S, cached = true)

   @test elem_type(T) == Generic.TotFrac{elem_type(S)}
   @test elem_type(Generic.TotFracRing{elem_type(S)}) == Generic.TotFrac{elem_type(S)}
   @test parent_type(Generic.TotFrac{elem_type(S)}) == Generic.TotFracRing{elem_type(S)}

   @test isa(T, Generic.TotFracRing)

   @test isa(T(3), Generic.TotFrac)

   @test isa(T(BigInt(7)), Generic.TotFrac)

   @test isa(T(3, 7), Generic.TotFrac)

   @test isa(T(S(3), 7), Generic.TotFrac)

   @test isa(T(3, S(7)), Generic.TotFrac)

   @test isa(T(), Generic.TotFrac)

   @test isa(T(S(5)), Generic.TotFrac)

   @test isa(T(T(5)), Generic.TotFrac)

   @test isa(T(3, 5)//7, Generic.TotFrac)

   @test isa(5//T(7, 5), Generic.TotFrac)

   @test isa(T(3, 5)//S(7), Generic.TotFrac)

   @test isa(S(5)//T(7, 5), Generic.TotFrac)
end

@testset "Generic.TotFrac.printing" begin
   S, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
   K = total_ring_of_fractions(S)

   @test string(K(x+y, z)) == "(x + y)//z"
   @test string(K(x, y)//z) == "x//(y*z)"
end

@testset "Generic.TotFrac.rand" begin
   S, = residue_ring(ZZ, 12)
   K = total_ring_of_fractions(S)

   test_rand(K, 0:11)
end

@testset "Generic.TotFrac.manipulation" begin
   R, = residue_ring(ZZ, 12)
   S = total_ring_of_fractions(R)

   @test denominator(S(5, 7)) == 7

   @test numerator(S(5, 7)) == 5

   @test iszero(zero(S))

   @test isone(one(S))

   @test canonical_unit(S(3)) == 1

   @test is_unit(S(5, 7))
   @test !is_unit(S(3, 7))

   a = S(5, 7)
   @test !is_zero_divisor(a)
   (f, h) = is_zero_divisor_with_annihilator(a)
   @test !f
   @test parent(h) == S

   a = S(3, 7)
   @test is_zero_divisor(a)
   (f, h) = is_zero_divisor_with_annihilator(a)
   @test f
   @test parent(h) == S
   @test !is_zero(h)
   @test is_zero(h*a)

   @test deepcopy(S(5, 7)) == S(5, 7)

   @test characteristic(S) == 12
end

@testset "Generic.TotFrac.unary_ops" begin
   R, = residue_ring(ZZ, 12)
   S = total_ring_of_fractions(R)

   @test -(S(5, 7)) == 1
end

@testset "Generic.TotFrac.binary_ops" begin
   R, = residue_ring(ZZ, 12)
   S = total_ring_of_fractions(R)

   for iter = 1:100
      a = rand(S, 0:11)
      b = rand(S, 0:11)
      c = rand(S, 0:11)

      @test a + b - b == a
      @test c*(a + b) == c*a + c*b
      @test c*(a - b) == c*(a - b)
      @test a - b == -(b - a)
   end
end

@testset "Generic.TotFrac.adhoc_binary" begin
   R, = residue_ring(ZZ, 12)
   S = total_ring_of_fractions(R)

   a = S(3, 7)
   b = S(1, 5)

   @test a + 2 == S(5, 7)

   @test a - 2 == S(1, 7)

   @test 2 + a == S(5, 7)

   @test 2 - a == S(11, 7)

   @test a + R(2) == S(5, 7)

   @test a - R(2) == S(1, 7)

   @test R(2) + a == S(5, 7)

   @test R(2) - a == S(11, 7)

   @test b*2 == S(2, 5)

   @test 2*b == S(2, 5)

   @test 7*a == 3

   @test a*7 == 3

   @test b*R(2) == S(2, 5)

   @test R(2)*b == S(2, 5)

   @test R(7)*a == 3

   @test a*R(7) == 3
end

@testset "Generic.TotFrac.comparison" begin
   R, = residue_ring(ZZ, 12)
   S = total_ring_of_fractions(R)

   a = S(3, 7)

   @test a == S(3, 7)
   @test S(10, 5) == S(2, 7)   

   @test isequal(a, S(3, 7))
end

@testset "Generic.TotFrac.adhoc_comparison" begin
   R, = residue_ring(ZZ, 12)
   S = total_ring_of_fractions(R)

   @test S(3, 1) == 3

   @test 3 == S(3, 1)

   @test S(3, 1) == R(3)

   @test R(3) == S(3, 1)
end

@testset "Generic.TotFrac.powering" begin
   R, = residue_ring(ZZ, 12)
   S = total_ring_of_fractions(R)

   a = S(5, 7)

   @test a^-12 == 1
end

@testset "Generic.TotFrac.inversion" begin
   R, = residue_ring(ZZ, 12)
   S = total_ring_of_fractions(R)

   a = S(5, 7)

   @test inv(a)*a == 1

   @test_throws(NotInvertibleError, inv(S(3, 7)))
end

@testset "Generic.TotFrac.promotion" begin
   R, = residue_ring(ZZ, 12)
   S = total_ring_of_fractions(R)

   T = elem_type(S)
   @test AbstractAlgebra.promote_rule(T, T) == T
end

@testset "Generic.TotFrac.unsafe_operators" begin
   R, = residue_ring(ZZ, 12)
   S = total_ring_of_fractions(R)

   a = S(5, 7)
   a = zero!(a)

   @test iszero(a)

   a = S(5, 7)
   b = S(1, 5)
   c = S(3, 5)

   c = mul!(c, a, b)

   @test c == 7

   c = mul!(c, b, a)

   @test c == 7

   a = S(5, 7)
   b = S(3, 1)
   c = S(3, 5)

   c = mul!(c, a, b)

   @test c == S(3, 7)

   a = S(5, 7)
   b = S(1, 7)

   a = addeq!(a, b)

   @test a == S(6, 7)

   a = S(3, 1)

   a = addeq!(a, a)

   @test a == S(6, 1)

   a = S(3, 1)
   b = S(5, 7)

   a = addeq!(a, b)

   @test a == S(2, 7)

   a = S(5, 7)
   b = S(3, 1)

   a = addeq!(a, b)

   @test a == S(2, 7)

   a = S(3, 5)
   b = S(1, 7)

   a = addeq!(a, b)

   @test a == S(2, 11)

   c = S(2, 7)
   a = S(5, 7)
   b = S(1, 7)

   c = add!(c, a, b)

   @test c == S(6, 7)

   a = S(3, 1)
   b = S(5, 7)

   c = add!(c, a, b)

   @test c == S(2, 7)

   a = S(5, 7)
   b = S(3, 1)

   c = add!(c, a, b)

   @test c == S(2, 7)

   a = S(3, 5)
   b = S(1, 7)

   c = add!(c, a, b)

   @test c == S(2, 11)
end

