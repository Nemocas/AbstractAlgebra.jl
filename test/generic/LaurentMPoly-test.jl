function test_elem(R::AbstractAlgebra.LaurentMPolyRing{BigInt})
    n = rand(1:10)
    # R: length between 1 and 9
    # R: exponents between -n and n
    # ZZ: coeffs between -99 and 99
    rand(R, 1:9, -n:n, -99:99)
end

@testset "Generic.LaurentMPoly.conformance" begin
    L, (x, y) = LaurentPolynomialRing(ZZ, ["x", "y"])
    test_Ring_interface(L)
    test_Ring_interface_recursive(L)
end

@testset "Generic.LaurentMPoly.constructors" begin
    L, (x, y) = LaurentPolynomialRing(GF(5), 2, "x", cached = true)
    @test L != LaurentPolynomialRing(GF(5), 2, 'x', cached = false)[1]
    @test L == LaurentPolynomialRing(GF(5), 2, :x, cached = true)[1]

    L, (x, y) = LaurentPolynomialRing(GF(5), ["x", "y"])
    @test L == LaurentPolynomialRing(GF(5), ['x', 'y'])[1]
    @test L != LaurentPolynomialRing(GF(5), [:x, :y], cached = false)[1]

    #only works because of the caching
    R, (X, Y) = PolynomialRing(coefficient_ring(L), symbols(L))
    @test one(L) == L(one(R))
    @test x == L(X)
    @test y == L(Y)
end

@testset "Generic.LaurentMPoly.derivative" begin
    L, (x, y) = LaurentPolynomialRing(ZZ, ["x", "y"])

    @test derivative(x, x) == 1
    @test derivative(y, x) == 0
    @test derivative(y, y) == 1
    @test derivative(x^2*y + x + x^3*y^2 - y, x) == 2*x*y + 1 + 3*x^2*y^2
    @test derivative(x^-2*y + x + x^3*y^2 - y, x) == -2*x^-3*y + 1 + 3*x^2*y^2
end

@testset "Generic.LaurentMPoly.euclidean" begin
    L, (x, y) = LaurentPolynomialRing(ZZ, ["x", "y"])
    @test isone(gcd(x, y))
    @test_throws Exception divrem(x, y)
end

@testset "Generic.LaurentMPoly.mpoly" begin
    L, (x, y) = LaurentPolynomialRing(ZZ, ["x", "y"])

    @test isgen(x)
    @test isgen(y)
    @test !isgen(one(L))
    @test !isgen(inv(x))
    @test !isgen(inv(y))
    @test !isgen(x*y)

    @test iszero(L(elem_type(coefficient_ring(L))[], Vector{Int}[]))

    a = x^-9 + y^9
    @test divides(a^2, a) == (true, a)

    @test evaluate(x^-1 + y, [QQ(2), QQ(3)]) == 1//2 + 3

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

    b = MPolyBuildCtx(L)
    for (c, e) in zip(coefficients(a), exponent_vectors(a))
        push_term!(b, c, e + [1, 2])
    end
    @test a*x*y^2 == finish(b)
    @test iszero(finish(b))

    @test map_coefficients(x->x^2, a) == 4*x^-2*y + 9*x*y^-3

    Q, (X, Y) = LaurentPolynomialRing(QQ, ["x", "y"])
    @test change_base_ring(QQ, a) == 2*X^-2*Y + 3*X*Y^-3
end
