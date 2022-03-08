@testset "Generic.FactoredFrac.ZZ" begin
    test_reps = 20

    FF = FactoredFractionField(ZZ)

    function randomFF(limit::Int)
        t = one(FF)
        for i in 1:abs(rand(Int)%limit)
            s = FF(rand(Int)%(20*limit))
            e = rand(Int)%limit
            t *= iszero(s) ? s^abs(e) : s^e
        end
        return t
    end

    for i in 1:test_reps
        local a, b
        a = randomFF(10)
        b = a - FF(1)
        if isone(a)
            @test iszero(b)
        else
            @test a != FF(1)
        end
        @test a == numerator(a)//denominator(a)
        @test AbstractAlgebra.Generic._bases_are_nice(a)
        @test AbstractAlgebra.Generic._bases_are_nice(b)
    end

    for i in 1:test_reps
        local a, b, c
        a = randomFF(10)
        b = randomFF(10)
        c = a*b
        if iszero(b)
            d = c
            @test iszero(c)
        else
            d = divexact(c, b)
            @test d == a
            if !iszero(a)
                e = divexact(divexact(c, a), b)
                @test isone(e)
            end
        end
        @test AbstractAlgebra.Generic._bases_are_nice(a)
        @test AbstractAlgebra.Generic._bases_are_nice(b)
        @test AbstractAlgebra.Generic._bases_are_nice(c)
        @test AbstractAlgebra.Generic._bases_are_nice(d)
    end

    for i in 1:test_reps
        local a, b, c, d
        a = randomFF(10)
        b = randomFF(10)
        c = a + b
        d = c - b
        @test d == a
        @test AbstractAlgebra.Generic._bases_are_nice(a)
        @test AbstractAlgebra.Generic._bases_are_nice(b)
        @test AbstractAlgebra.Generic._bases_are_nice(c)
        @test AbstractAlgebra.Generic._bases_are_nice(d)
    end
end

@testset "Generic.FactoredFrac.ZZx" begin
    Zx, x = PolynomialRing(ZZ, "x")
    F = FactoredFractionField(Zx)
    x = F(x)
    a = divexact(x + 1, (x + 2)^2) + x
    @test (a - x)*(x + 1)^-1 == (x + 2)^-2
end

@testset "Generic.FactoredFrac.ZZxyz" begin
    Zxyz, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])
    F = FactoredFractionField(Zxyz)
    (x, y, z) = (F(x), F(y), F(z))
    @test (inv(x+y) + inv(x+y+z))*(x+y)*(x+y+z) == 2*x+2*y+z
    @test det(matrix(F, [x y z; x^2 y^2 z^2; x^3 y^3 z^3])) == -x*y*z*(y - z)*(x - y)*(x - z)
end

@testset "Generic.FactoredFrac.constructors" begin
    S, x = PolynomialRing(ZZ, "x")
    T = FactoredFractionField(S)

    @test FactoredFractionField(S, cached = true) === FactoredFractionField(S, cached = true)
    @test FactoredFractionField(S, cached = false) !== FactoredFractionField(S, cached = true)

    @test elem_type(T) == Generic.FactoredFrac{elem_type(S)}
    @test elem_type(Generic.FactoredFracField{elem_type(S)}) == Generic.FactoredFrac{elem_type(S)}
    @test parent_type(Generic.FactoredFrac{elem_type(S)}) == Generic.FactoredFracField{elem_type(S)}

    @test isa(T, Generic.FactoredFracField)
    @test isa(T(3), Generic.FactoredFrac)
    @test isa(T(BigInt(7)), Generic.FactoredFrac)
    @test isa(T(x + 2), Generic.FactoredFrac)
    @test isa(T(3, 7), Generic.FactoredFrac)
    @test isa(T(x + 2, x + 1), Generic.FactoredFrac)
    @test isa(T(x + 2, 4), Generic.FactoredFrac)
    @test isa(T(3, x + 1), Generic.FactoredFrac)
    @test isa(T(T(x + 2)), Generic.FactoredFrac)
    @test isa(T(), Generic.FactoredFrac)
    @test isa(T(x + 3)//(x^2 + 2), Generic.FactoredFrac)
    @test isa(T(x + 3)//12, Generic.FactoredFrac)
    @test isa(12//(x + 2), Generic.Frac)
    @test isa(12//(T(x) + 2), Generic.FactoredFrac)
    @test isa((x + 1)//T(x + 2, x + 1), Generic.FactoredFrac)
    @test isa(T(x + 2, x + 1)//(x + 1), Generic.FactoredFrac)
    @test isa(T(x + 2, x + 1)//T(x, x + 2), Generic.FactoredFrac)

    TT = FactoredFractionField(PolynomialRing(QQ, "x")[1])
    a = TT(1)
    b = T(2)

    @test a in [a, b]
    @test a in [b, a]
    @test !(a in [b])
    @test a in keys(Dict(a => 1))
    @test !(b in keys(Dict(a => 1)))
end

@testset "Generic.FactoredFrac.printing" begin
    S, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
    F = FactoredFractionField(S)
    (x, y, z) = (F(x), F(y), F(z))

    s = string(x//y//z)
    @test s == "x/(y*z)" || s == "x/(z*y)"

    s = string((x+y)*(z+1)//z)
    @test s == "(x + y)*(z + 1)/z" || s == "(z + 1)*(x + y)/z"
end

@testset "Generic.FactoredFrac.derivative" begin
    S, (x, y) = PolynomialRing(ZZ, ["x", "y"])
    F = FactoredFractionField(S)
    (X, Y) = (F(x), F(y))

    @test derivative((1+X)*(2+X)^2//(3+X)^3//(4+X)^4, x) ==
                              -(2*(2+X)*(11*X+10*X^2+2*X^3))//((3+X)^4*(4+X)^5)
end

