@testset "Generic.LaurentMPoly.conformance" begin
    L, (x, y) = laurent_polynomial_ring(ZZ, ["x", "y"])
    ConformanceTests.test_Ring_interface(L)
    ConformanceTests.test_Ring_interface_recursive(L)

    L, (x, y) = laurent_polynomial_ring(residue_ring(ZZ, ZZ(6))[1], ["x", "y"])
    ConformanceTests.test_Ring_interface(L)
end

@testset "Generic.LaurentMPoly.constructors" begin
    L, (x, y) = laurent_polynomial_ring(GF(5), 2, "x", cached = true)
    @test L != laurent_polynomial_ring(GF(5), 2, 'x', cached = false)[1]
    @test L == laurent_polynomial_ring(GF(5), 2, :x, cached = true)[1]

    @test base_ring(L) == GF(5)
    @test coefficient_ring(L) == GF(5)
    @test coefficient_ring_type(L) === typeof(GF(5))
    @test base_ring(x) == GF(5)
    @test coefficient_ring(x) == GF(5)
    @test coefficient_ring_type(x) === typeof(GF(5))

    L, (x, y) = laurent_polynomial_ring(GF(5), ["x", "y"])
    @test L == laurent_polynomial_ring(GF(5), ['x', 'y'])[1]
    @test L != laurent_polynomial_ring(GF(5), [:x, :y], cached = false)[1]

    # only works because of the caching
    R, (X, Y) = polynomial_ring(coefficient_ring(L), symbols(L))
    @test one(L) == L(one(R))
    @test x == L(X)
    @test y == L(Y)
    @test X + x == 2*x
end

@testset "Generic.LaurentMPoly.characteristic" for R in (GF(5), ZZ, residue_ring(ZZ, 6)[1])
   L, (x, y) = laurent_polynomial_ring(R, 2, "x", cached = true)
   @test characteristic(L) == characteristic(R)
end

@testset "Generic.LaurentMPoly.printing" begin
   R, (x,) = laurent_polynomial_ring(residue_ring(ZZ, 6)[1], ["x"])
   @test !occursin("\n", sprint(show, R))
end

@testset "Generic.LaurentMPoly.derivative" begin
    L, (x, y) = laurent_polynomial_ring(ZZ, ["x", "y"])

    @test derivative(x, x) == 1
    @test derivative(y, x) == 0
    @test derivative(y, y) == 1
    @test derivative(x^2*y + x + x^3*y^2 - y, x) == 2*x*y + 1 + 3*x^2*y^2
    @test derivative(x^-2*y + x + x^3*y^2 - y, x) == -2*x^-3*y + 1 + 3*x^2*y^2
end

@testset "Generic.LaurentMPoly.euclidean" begin
    L, (x, y) = laurent_polynomial_ring(ZZ, ["x", "y"])
    @test isone(gcd(x, y))
    @test isone(gcd(inv(x), inv(y)))
    @test_throws Exception divrem(x, y)
end

@testset "Generic.LaurentMPoly.mpoly" begin
    L, (x, y) = laurent_polynomial_ring(ZZ, ["x", "y"])

    @test is_gen(x)
    @test is_gen(y)
    @test !is_gen(one(L))
    @test !is_gen(inv(x))
    @test !is_gen(inv(y))
    @test !is_gen(x*y)

    @test iszero(L(elem_type(coefficient_ring(L))[], Vector{Int}[]))

    a = x^-9 + y^9
    @test divides(a^2, a) == (true, a)

    @test evaluate(x^-1 + y, [QQ(2), QQ(3)]) == 1//2 + 3
    @test evaluate(inv(x)*x, [QQ(2), QQ(3)]) == 1

    a = 2*x^-2*y + 3*x*y^-3
    le = leading_exponent_vector(a)
    @test le == [-2, 1] || le == [1, -3]
    @test le == leading_exponent_vector(leading_monomial(a))
    @test le == leading_exponent_vector(leading_term(a))
    @test leading_term(a) == leading_coefficient(a)*leading_monomial(a)
    @test leading_coefficient(a) == first(collect(coefficients(a)))
    @test leading_monomial(a) == first(collect(monomials(a)))
    @test leading_term(a) == first(collect(terms(a)))
    @test a == sum(terms(a))
    @test a == sum(coefficients(a) .* monomials(a))
    @test a == L(collect(coefficients(a)), collect(exponent_vectors(a)))
    @test iszero((@inferred constant_coefficient(a)))
    @test coeff(a, [-2, 1]) == 2
    @test coeff(a, [1, -3]) == 3
    @test coeff(a, [1, 1]) == 0

    b = MPolyBuildCtx(L)
    for (c, e) in zip(coefficients(a), exponent_vectors(a))
        push_term!(b, c, e + [1, 2])
    end
    @test a*x*y^2 == finish(b)
    @test iszero(finish(b))

    @test map_coefficients(x->x^2, a) == 4*x^-2*y + 9*x*y^-3

    Q, (X, Y) = laurent_polynomial_ring(QQ, ["x", "y"])
    @test change_base_ring(QQ, a) == 2*X^-2*Y + 3*X*Y^-3

    b = MPolyBuildCtx(L)
    push_term!(b, 1, [0, 0])
    push_term!(b, 2, [-2, -1])
    p = @inferred finish(b)
    @test p == 1 + 2 * x^-2 * y^-1
    @test constant_coefficient(p) == 1

    p = inv(inv(x))
    @test constant_coefficient(p) == 0

    @test is_monomial(x^2*y^-2)
    @test !is_monomial(2*x^2*y^-2)
    @test !is_monomial(x+y)
    @test is_term(x^2*y^-2)
    @test is_term(2*x^2*y^-2)
    @test !is_term(x+y)

    @test inflate(y+y^-1, [1, 2]) == y^2+y^-2
    @test inflate(y+y^-1, [0, 3], [1, 2]) == y^5+y
    @test inflate(y+y^-1, [2], [3], [2]) == y^5+y
    @test inflate(x*y^2+y, [2, 2]) == x^2*y^4+y^2
    @test inflate(x*y^2+y, [1, 2], [2, 2]) == x^3*y^6+x*y^4
end


# -------------------------------------------------------

# Coeff rings for the tests below
ZeroRing,_ = residue_ring(ZZ,1);
ZZmod720,_ = residue_ring(ZZ, 720);

# [2024-12-12  laurent_polynomial_ring currently gives error when coeff ring is zero ring]
# ## LaurentMPoly over ZeroRing
# @testset "Nilpotent/unit for ZeroRing[x,y, x^(-1),y^(-1)]" begin
#   P,(x,y) = laurent_polynomial_ring(ZeroRing, ["x","y"]);
#   @test is_nilpotent(P(0))
#   @test is_nilpotent(P(1))
#   @test is_nilpotent(x)
#   @test is_nilpotent(-x)
#   @test is_nilpotent(x+y)
#   @test is_nilpotent(x-y)
#   @test is_nilpotent(x*y)

#   @test is_unit(P(0))
#   @test is_unit(P(1))
#   @test is_unit(x)
#   @test is_unit(-x)
#   @test is_unit(x+y)
#   @test is_unit(x-y)
#   @test is_unit(x*y)
# end

## LaurentMPoly over ZZ
@testset "Nilpotent/unit for ZZ[x,y, x^(-1),y^(-1)]" begin
  P,(x,y) = laurent_polynomial_ring(ZZ, ["x","y"]);
  @test is_nilpotent(P(0))
  @test !is_nilpotent(P(1))
  @test !is_nilpotent(x)
  @test !is_nilpotent(-x)
  @test !is_nilpotent(x+y)
  @test !is_nilpotent(x-y)
  @test !is_nilpotent(x*y)

  @test !is_unit(P(0))
  @test is_unit(P(1))
  @test is_unit(P(-1))
  @test !is_unit(P(-2))
  @test !is_unit(P(-2))
  @test is_unit(x)
  @test is_unit(-x)
  @test is_unit(1/x)
  @test is_unit(-1/x)
  @test !is_unit(2/x)
  @test !is_unit(-2/x)
  @test !is_unit(x+1)
  @test !is_unit(x-1)
  @test !is_unit(x+y)
  @test !is_unit(x-y)
  @test is_unit(x*y)
end

## LaurentMPoly over QQ
@testset "Nilpotent/unit for QQ[x,y, x^(-1),y^(-1)]" begin
  P,(x,y) = laurent_polynomial_ring(QQ, ["x","y"]);
  @test is_nilpotent(P(0))
  @test !is_nilpotent(P(1))
  @test !is_nilpotent(x)
  @test !is_nilpotent(-x)
  @test !is_nilpotent(x+y)
  @test !is_nilpotent(x-y)
  @test !is_nilpotent(x*y)

  @test !is_unit(P(0))
  @test is_unit(P(1))
  @test is_unit(P(-1))
  @test is_unit(P(2))
  @test is_unit(P(-2))
  @test is_unit(x)
  @test is_unit(-x)
  @test is_unit(2*x)
  @test is_unit(-2*x)
  @test is_unit(1/x)
  @test is_unit(-1/x)
  @test is_unit(2/x)
  @test is_unit(-2/x)
  @test !is_unit(x+1)
  @test !is_unit(x-1)
  @test !is_unit(x+y)
  @test !is_unit(x-y)
  @test is_unit(x*y)
end

## LaurentMPoly over ZZ/720
@testset "Nilpotent/unit for ZZ/(720)[x,y, x^(-1), y^(-1)]" begin
  P,(x,y) = laurent_polynomial_ring(ZZmod720, ["x","y"]);
  @test is_nilpotent(P(0))
  @test !is_nilpotent(P(1))
  @test is_nilpotent(P(30))
  @test !is_nilpotent(x)
  @test !is_nilpotent(-x)
  @test is_nilpotent(30*x)
  @test is_nilpotent(30/x)
  @test is_nilpotent(30*x+120*y)
  @test is_nilpotent(30*x-120*y)
  @test !is_nilpotent(x*y)
  @test is_nilpotent(30*x*y)
  @test is_nilpotent(30*x/y)

  @test !is_unit(P(0))
  @test is_unit(P(1))
  @test is_unit(P(-1))
  @test !is_unit(P(2))
  @test !is_unit(P(-2))
  @test is_unit(P(7))
  @test is_unit(P(-7))
  @test is_unit(x)
  @test is_unit(-x)
  @test !is_unit(35*x)
  @test !is_unit(35/x)
  @test !is_unit(30*x)
  @test !is_unit(30/x)
  @test !is_unit(x+1)
  @test !is_unit(x-1)
  @test is_unit(x+30)
  @test is_unit(x-30)
  @test is_unit(1+30*x)
  @test is_unit(1-30*x)
  @test is_unit(7+60*x)
  @test is_unit(7-60*x)
  @test is_unit(600+7*x+30*x^2)
  @test is_unit(600-7*x+30*x^2)
  @test is_unit(x+30/y)
  @test is_unit(x-30/y)
  @test is_unit(1+30*x/y)
  @test is_unit(1-30*x*y)
  @test is_unit(7+60*x+210/y)
  @test is_unit(7-60*x+210/y)
  @test is_unit(600+7*x/y+30*x^2)
  @test is_unit(600-7*x*y+30*x^2)
  @test !is_unit(30*x+120*y)
  @test !is_unit(30*x-120*y)
  @test is_unit(x*y)
  @test !is_unit(30*x*y)
  @test !is_unit(30*x/y)
end
