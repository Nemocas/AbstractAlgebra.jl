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

    @test is_domain_type(L)
    @test is_domain_type(elem_type(L))
    @test is_exact_type(elem_type(L))
    @test !is_exact_type(elem_type(laurent_polynomial_ring(RealField, ["a"])[1]))
    @test !is_univariate(L)
    @test is_univariate(laurent_polynomial_ring(GF(5), ["x"])[1])

    @test coefficient_ring(L) == GF(5)
    @test coefficient_ring_type(L) === typeof(GF(5))
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

@testset "Generic.LaurentMPoly.more_mpoly" begin
    L, (x, y, z) = laurent_polynomial_ring(ZZ, ["x", "y", "z"])
    f = 2*x^-2*y + 3*x*y^-3 - 4

    @test degree(f, 1) == 1
    @test degree(f, y) == 1
    @test degree(f, 3) == 0
    @test degree(x^-2*y^-1 + x^-5, 1) == -2
    @test degrees(f) == [1, 1, 0]
    @test degrees(x^-2 + y^-3) == [0, 0, 0]
    @test total_degree(f) == 0
    @test total_degree(x^-2*y^-1 + x^-5) == -3
    @test_throws ArgumentError degree(zero(L), 1)
    @test_throws ArgumentError degrees(zero(L))
    @test_throws ArgumentError total_degree(zero(L))

    @test var_indices(f) == [1, 2]
    @test vars(f) == [x, y]
    @test vars(3*x*inv(x)) == []
    @test is_univariate(x^-3 + 2x)
    @test is_univariate(L(7))
    @test !is_univariate(f)

    @test tail(f) == f - leading_term(f)
    @test iszero(tail(x^-1))
    @test iszero(tail(zero(L)))
    @test trailing_coefficient(f) == last(collect(coefficients(f)))
    @test_throws ArgumentError trailing_coefficient(zero(L))

    @test content(f) == 1
    @test content(6*x^-1 - 9*y) == 3
    @test is_homogeneous(x^-1*y^2 + z)
    @test !is_homogeneous(f)

    Q, (u, v, w) = laurent_polynomial_ring(QQ, ["u", "v", "w"])
    g = u^-2*v + 3*u*v^-1*w
    @test evaluate(g, [1], [QQ(2)]) == 1//4*v + 6*v^-1*w
    @test evaluate(g, [2, 3], [QQ(1//2), QQ(3)]) == 1//2*u^-2 + 18*u
    @test evaluate(g, [v], [u]) == u^-1 + 3*w
    @test evaluate(g, [v, w], [w, v]) == u^-2*w + 3*u*w^-1*v
    @test evaluate(g, Int[], elem_type(QQ)[]) == g
    @test_throws ArgumentError evaluate(g, [1, 1], [QQ(1), QQ(2)])
    @test_throws ArgumentError evaluate(g, [1], [QQ(1), QQ(2)])
    @test_throws ArgumentError evaluate(g, [4], [QQ(1)])
    @test_throws ErrorException evaluate(x^-1, [1], [ZZ(2)])

    @test change_coefficient_ring(QQ, f) == change_base_ring(QQ, f)
    @test change_coefficient_ring(QQ, f; parent = Q) == 2*u^-2*v + 3*u*v^-3 - 4

    h = x^-3*y^2 + x*y^6 + x^5*y^-2
    @test deflation(h) == ([-3, -2, 0], [4, 4, 0])
    shift, defl = deflation(h)
    @test inflate(deflate(h, shift, defl), shift, defl) == h
    @test deflate(x^-2 + y^4, [2, 2, 1]) == x^-1 + y^2
    @test deflation(zero(L)) == ([0, 0, 0], [0, 0, 0])

    s = x^-1*y + 2 - z^3
    @test is_square(s^2)
    @test is_square(x^-2*y^4)
    @test is_square(zero(L))
    @test !is_square(x^-1)
    @test !is_square(-x^2)
    @test !is_square(s^2*x)
    @test sqrt(s^2) == s || sqrt(s^2) == -s
    @test sqrt(4*x^-2*y^4) == 2*x^-1*y^2 || sqrt(4*x^-2*y^4) == -2*x^-1*y^2
    @test_throws ErrorException sqrt(x)
    ok, r = is_square_with_sqrt(s^2*x^-4)
    @test ok && r^2 == s^2*x^-4

    @test isless(x^-1, x) == isless(one(L), x^2)
    @test isless(y^-1, x^-1) == isless(x, y)
    @test !isless(x^-1*y, x^-1*y)
    @test_throws ErrorException isless(x + 1, y)
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
  @test_broken !is_unit(x+1)
  @test_broken !is_unit(x-1)
  @test_broken is_unit(x+30)
  @test_broken is_unit(x-30)
  @test_broken is_unit(1+30*x)
  @test_broken is_unit(1-30*x)
  @test_broken is_unit(7+60*x)
  @test_broken is_unit(7-60*x)
  @test_broken is_unit(600+7*x+30*x^2)
  @test_broken is_unit(600-7*x+30*x^2)
  @test_broken is_unit(x+30/y)
  @test_broken is_unit(x-30/y)
  @test_broken is_unit(1+30*x/y)
  @test_broken is_unit(1-30*x*y)
  @test_broken is_unit(7+60*x+210/y)
  @test_broken is_unit(7-60*x+210/y)
  @test_broken is_unit(600+7*x/y+30*x^2)
  @test_broken is_unit(600-7*x*y+30*x^2)
  @test_broken !is_unit(30*x+120*y)
  @test_broken !is_unit(30*x-120*y)
  @test is_unit(x*y)
  @test !is_unit(30*x*y)
  @test !is_unit(30*x/y)
end

@testset "#2378" begin
  L1, (a, b) = laurent_polynomial_ring(QQ, [:a, :b]);
  R2, (k, l) = laurent_polynomial_ring(L1, [:k, :l]);
  @test a * k == R2(a) * k
end
