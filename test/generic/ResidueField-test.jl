@testset "Generic.ResF.constructors..." begin
   B = ZZ

   R = Generic.ResidueField(B, 16453889)

   @test elem_type(R) == Generic.ResF{elem_type(B)}
   @test elem_type(Generic.ResField{elem_type(B)}) == Generic.ResF{elem_type(B)}
   @test parent_type(Generic.ResF{elem_type(B)}) == Generic.ResField{elem_type(B)}

   @test isa(R, Generic.ResField)

   a = R(123)

   @test isa(a, Generic.ResF)

   b = R(a)

   @test isa(b, Generic.ResF)

   c = R(ZZ(12))

   @test isa(c, Generic.ResF)

   d = R()

   @test isa(d, Generic.ResF)

   S, x = PolynomialRing(R, "x")
   T = ResidueField(S, x^3 + 3x + 1)

   @test isa(T, Generic.ResField)

   f = T(x^4)

   @test isa(f, Generic.ResF)

   g = T(f)

   @test isa(g, Generic.ResF)

   S = Generic.ResidueRing(B, 2)
   x = R(1)
   y = S(1)
   @test x in [x, y]
   @test x in [y, x]
   @test !(x in [y])

   @test x in keys(Dict(x => 1))
   @test !(y in keys(Dict(x => 1)))
end

@testset "Generic.ResF.rand..." begin
   R = Generic.ResidueField(ZZ, 16453889)
   f = rand(R, 1:9)
   @test f isa Generic.ResF
   f = rand(rng, R, 1:9)
   @test f isa Generic.ResF
end

@testset "Generic.ResF.manipulation..." begin
   R = Generic.ResidueField(ZZ, 16453889)

   @test modulus(R) == 16453889

   g = zero(R)

   @test iszero(g)

   @test modulus(g) == 16453889

   S, x = PolynomialRing(R, "x")
   T = ResidueField(S, x^3 + 3x + 1)

   h = one(T)

   @test isunit(h)

   @test isone(h)

   @test data(h) == 1

   @test canonical_unit(R(11)) == R(11)

   @test canonical_unit(T(x + 1)) == T(x + 1)

   @test deepcopy(h) == h
end

@testset "Generic.ResF.unary_ops..." begin
   R = Generic.ResidueField(ZZ, 16453889)

   @test -R(12345) == R(16441544)

   S, x = PolynomialRing(R, "x")
   T = ResidueField(S, x^3 + 3x + 1)

   @test -T(x^5 + 1) == T(x^2+16453880*x+16453885)
end

@testset "Generic.ResF.binary_ops..." begin
   R = Generic.ResidueField(ZZ, 13)

   f = R(4)
   g = R(6)

   @test f + g == R(10)

   @test f - g == R(11)

   @test f*g == R(11)

   Q = Generic.ResidueField(ZZ, 7)
   S, x = PolynomialRing(Q, "x")
   T = ResidueField(S, x^3 + 3x + 1)

   n = T(x^5 + 1)
   p = T(x^2 + 2x + 1)

   @test n + p == T(4x + 5)

   @test n - p == T(5x^2 + 3)

   @test n*p == T(3x^2 + 4x + 4)
end

@testset "Generic.ResF.gcd..." begin
   R = Generic.ResidueField(ZZ, 13)

   f = R(4)
   g = R(6)

   @test gcd(f, g) == R(1)

   Q = Generic.ResidueField(ZZ, 7)
   S, x = PolynomialRing(Q, "x")
   T = ResidueField(S, x^3 + 3x + 1)

   n = T(x^5 + 1)
   p = T(x^2 + 2x + 1)

   @test gcd(n, p) == 1
end

@testset "Generic.ResF.adhoc_binary..." begin
   R = Generic.ResidueField(ZZ, 7)

   a = R(3)

   @test a + 3 == R(6)

   @test 3 - a == R(0)

   @test 5a == R(1)

   S, x = PolynomialRing(R, "x")
   T = ResidueField(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test f + 4 == T(x^5 + 5)

   @test 4 - f == T(x^2+5*x)

   @test f*5 == T(2*x^2+3*x+6)
end

@testset "Generic.ResF.comparison..." begin
   R = Generic.ResidueField(ZZ, 7)

   a = R(3)
   b = a
   c = R(2)

   @test b == a

   @test isequal(b, a)

   @test c != a

   S, x = PolynomialRing(R, "x")
   T = ResidueField(S, x^3 + 3x + 1)

   f = T(x^5 + 1)
   g = 8f
   h = f + g

   @test f == g
   @test h != g

   @test isequal(f, g)
end

@testset "Generic.ResF.adhoc_comparison..." begin
   R = Generic.ResidueField(ZZ, 7)

   a = R(3)

   @test a == 3
   @test 4 != a

   S, x = PolynomialRing(R, "x")
   T = ResidueField(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test f != 5
end

@testset "Generic.ResF.powering..." begin
   R = Generic.ResidueField(ZZ, 7)

   a = R(3)

   @test a^5 == 5

   S, x = PolynomialRing(R, "x")
   T = ResidueField(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test f^100 == T(x^2 + 2x + 1)
end

@testset "Generic.ResF.inversion..." begin
   R = Generic.ResidueField(ZZ, 47)

   a = R(5)

   @test inv(a) == 19

   R = Generic.ResidueField(ZZ, 41)
   S, x = PolynomialRing(R, "x")
   T = ResidueField(S, x^3 + 3x + 1)

   f = T(x^5 + 1)

   @test inv(f) == T(26*x^2+31*x+10)
end

@testset "Generic.ResF.exact_division..." begin
   R = Generic.ResidueField(ZZ, 47)

   a = R(5)
   b = R(3)

   @test divexact(a, b) == 33

   R = Generic.ResidueField(ZZ, 41)
   S, x = PolynomialRing(R, "x")
   T = ResidueField(S, x^3 + 3x + 1)

   f = T(x^5 + 1)
   g = T(x^4 + x + 2)

   @test divexact(f, g) == T(7*x^2+25*x+26)
end
