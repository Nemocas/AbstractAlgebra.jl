function test_elem(FF::Generic.FactoredFracField{BigInt})
    limit = 10
    t = one(FF)
    for i in 1:abs(rand(Int)%limit)
        s = FF(rand(Int)%(20*limit))
        e = rand(Int)%limit
        t *= iszero(s) ? s^abs(e) : s^e
    end
    return t
end

@testset "Generic.FactoredFracFieldElem.ZZ.conformance" begin
    FF = FactoredFractionField(ZZ)
    test_Field_interface(FF)
end

@testset "Generic.FactoredFracFieldElem.ZZ.adhoc" begin
    FF = FactoredFractionField(ZZ)

    @test FF(2//3) == FF(2, 3)
    @test FF(2//-1) == FF(2, -1)
    @test FF(2, 3) + ZZ(1) == 5//3
    @test FF(2, 3) + ZZ(1) == 5//3
    @test FF(2, 3) - ZZ(1) == -1//3
    @test ZZ(1) + FF(2, 3) == 5//3
    @test ZZ(1) - FF(2, 3) == 1//3
    @test FF(2, 3) * ZZ(1) == 2//3
    @test ZZ(1) * FF(2, 3) == 2//3
    @test divexact(FF(2, 3), ZZ(1)) == 2//3
    @test divexact(ZZ(1), FF(2, 3)) == 3//2
    @test !iszero(gcd(FF(2//3), FF(4//5)))

    Zx, x = polynomial_ring(ZZ, "x")
    F = FactoredFractionField(Zx)

    @test F(2//3) == F(2, 3)
    @test F(2//-1) == F(2, -1)
    @test F(2, 3) + Zx(1) == 5//3
    @test F(2, 3) - Zx(1) == -1//3
    @test Zx(1) + F(2, 3) == 5//3
    @test Zx(1) - F(2, 3) == 1//3
    @test F(2, 3) * Zx(1) == 2//3
    @test Zx(1) * F(2, 3) == 2//3
    @test F(2, 3) * ZZ(1) == 2//3
    @test ZZ(1) * F(2, 3) == 2//3
    @test divexact(F(2, 3), Zx(1)) == 2//3
    @test divexact(Zx(1), F(2, 3)) == 3//2
    @test divexact(F(2, 3), ZZ(1)) == 2//3
    @test divexact(ZZ(1), F(2, 3)) == 3//2

    a = FF(-1//3)
    @test (a + 1//2) - 1//2 == a
    @test 1//2 - (1//2 + a) == -a
    @test divexact((a * 1//2), 1//2) == a
    @test divexact(1//2, 1//2 * a) == inv(a)
end

@testset "Generic.FactoredFracFieldElem.ZZ.evaluate" begin
    Zx, x = polynomial_ring(ZZ, "x")
    F = FactoredFractionField(Zx)
    x = F(x)
    @test evaluate(x//(x+1), 2//3) == 2//5

    Zxy, (x, y) = polynomial_ring(ZZ, ["x", "y"])
    F = FactoredFractionField(Zxy)
    (x, y) = (F(x), F(y))
    @test evaluate(x//(x+y)^2, [1//3, 1//2]) == 12//25
end

@testset "Generic.FactoredFracFieldElem.ZZ.valuation" begin
    Zxy, (x, y) = polynomial_ring(ZZ, ["x", "y"])
    F = FactoredFractionField(Zxy)
    (X, Y) = (F(x), F(y))
    @test remove(X//(X+Y)^2, x)   == (1, 1//(X+Y)^2)
    @test remove(X//(X+Y)^2, x+y) == (-2, X)
    @test remove(X//(X+Y)^2, x-y) == (0, X//(X+Y)^2)
end

@testset "Generic.FactoredFracFieldElem.ZZ._bases_are_nice" begin
    FF = FactoredFractionField(ZZ)
    test_reps = 20

    for i in 1:test_reps
        a = test_elem(FF)
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
        a = test_elem(FF)
        b = test_elem(FF)
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
        a = test_elem(FF)
        b = test_elem(FF)
        c = a + b
        d = c - b
        @test d == a
        @test AbstractAlgebra.Generic._bases_are_nice(a)
        @test AbstractAlgebra.Generic._bases_are_nice(b)
        @test AbstractAlgebra.Generic._bases_are_nice(c)
        @test AbstractAlgebra.Generic._bases_are_nice(d)
    end
end

@testset "Generic.FactoredFracFieldElem.ZZ.normalise" begin
    FF = FactoredFractionField(ZZ)

    a = one(FF)
    @test length(a) == 0

    push_term!(a, 4, 3)
    @test length(a) == 1
    @test a == 64

    push_term!(a, 2, -6)
    @test length(a) == 2
    @test isone(a)
    @test a == 1

    @test collect(a) == [(ZZ(4), 3), (ZZ(2), -6)]
    @test !AbstractAlgebra.Generic._bases_are_nice(a)

    a = normalise(a)
    @test length(a) == 0
    @test isone(a)
    @test is_unit(a)
    @test a == 1

    a = one(FF)
    push_term!(a, 4, 3)
    push_term!(a, 6, -3)
    a = normalise(a)
    @test a == 4^3//6^3
    @test length(a) == 2

    push_term!(a, 0, 3)
    a = normalise(a)
    @test iszero(a)
    @test length(a) == 0
end

@testset "Generic.FactoredFracFieldElem.ZZx" begin
    Zx, x = polynomial_ring(ZZ, "x")
    F = FactoredFractionField(Zx)
    x = F(x)
    a = divexact(x + 1, (x + 2)^2) + x
    @test (a - x)*(x + 1)^-1 == (x + 2)^-2
    @test F(2//3) == F(2, 3)
    @test F(2//gen(Zx)) == F(2, gen(Zx))
end

@testset "Generic.FactoredFracFieldElem.ZZxyz" begin
    Zxyz, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"])
    F = FactoredFractionField(Zxyz)
    (x, y, z) = (F(x), F(y), F(z))
    @test (inv(x+y) + inv(x+y+z))*(x+y)*(x+y+z) == 2*x+2*y+z
    @test det(matrix(F, [x y z; x^2 y^2 z^2; x^3 y^3 z^3])) == -x*y*z*(y - z)*(x - y)*(x - z)
end

@testset "Generic.FactoredFracFieldElem.constructors" begin
    S, x = polynomial_ring(ZZ, "x")
    T = FactoredFractionField(S)

    @test FactoredFractionField(S, cached = true) === FactoredFractionField(S, cached = true)
    @test FactoredFractionField(S, cached = false) !== FactoredFractionField(S, cached = true)

    @test elem_type(T) == Generic.FactoredFracFieldElem{elem_type(S)}
    @test elem_type(Generic.FactoredFracField{elem_type(S)}) == Generic.FactoredFracFieldElem{elem_type(S)}
    @test parent_type(Generic.FactoredFracFieldElem{elem_type(S)}) == Generic.FactoredFracField{elem_type(S)}

    @test isa(T, Generic.FactoredFracField)
    @test isa(T(3), Generic.FactoredFracFieldElem)
    @test isa(T(BigInt(7)), Generic.FactoredFracFieldElem)
    @test isa(T(x + 2), Generic.FactoredFracFieldElem)
    @test isa(T(3, 7), Generic.FactoredFracFieldElem)
    @test isa(T(x + 2, x + 1), Generic.FactoredFracFieldElem)
    @test T(x + 2, x + 1)*T(x + 1) == T(x + 2)
    @test isa(T(x + 2, 4), Generic.FactoredFracFieldElem)
    @test isa(T(3, x + 1), Generic.FactoredFracFieldElem)
    @test isa(T(T(x + 2)), Generic.FactoredFracFieldElem)
    @test isa(T(), Generic.FactoredFracFieldElem)
    @test isa(T(x + 3)//(x^2 + 2), Generic.FactoredFracFieldElem)
    @test isa(T(x + 3)//12, Generic.FactoredFracFieldElem)
    @test isa(12//(x + 2), Generic.FracFieldElem)
    @test isa(12//(T(x) + 2), Generic.FactoredFracFieldElem)
    @test isa((x + 1)//T(x + 2, x + 1), Generic.FactoredFracFieldElem)
    @test isa(T(x + 2, x + 1)//(x + 1), Generic.FactoredFracFieldElem)
    @test isa(T(x + 2, x + 1)//T(x, x + 2), Generic.FactoredFracFieldElem)

    @test characteristic(T) == 0

    TT = FactoredFractionField(polynomial_ring(QQ, "x")[1])
    a = TT(1)
    b = T(2)

    @test a in [a, b]
    @test a in [b, a]
    @test !(a in [b])
    @test a in keys(Dict(a => 1))
    @test !(b in keys(Dict(a => 1)))
end

@testset "Generic.FactoredFracFieldElem.printing" begin
    S, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
    F = FactoredFractionField(S)
    (x, y, z) = (F(x), F(y), F(z))

    s = string(x//y//z)
    @test s == "x/(y*z)" || s == "x/(z*y)"

    s = string((x+y)*(z+1)//z)
    @test s == "(x + y)*(z + 1)/z" || s == "(z + 1)*(x + y)/z"

    @test length(string(F)) > 3
end

@testset "Generic.FactoredFracFieldElem.derivative" begin
    S, (x, y) = polynomial_ring(ZZ, ["x", "y"])
    F = FactoredFractionField(S)
    (X, Y) = (F(x), F(y))

    @test derivative((1+X)*(2+X)^2//(3+X)^3//(4+X)^4, x) ==
                              -(2*(2+X)*(11*X+10*X^2+2*X^3))//((3+X)^4*(4+X)^5)
end

