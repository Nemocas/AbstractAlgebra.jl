using Test

import AbstractAlgebra: Generic.puiseux_polynomial_ring_elem
import AbstractAlgebra: Generic.rescale
import AbstractAlgebra: Generic.normalize!

@testset "PuiseuxPolynomials.jl" begin
    @testset "Construction" begin
        K, (t1,t2,t3) = laurent_polynomial_ring(QQ, ["t1","t2","t3"])
        K_p, (tp1,tp2,tp3) = puiseux_polynomial_ring(QQ, ["t1","t2","t3"])
        @test K_p.baseRing == K

        h = 1+t1 + 2*t2+3*t1^4+t1*t2^4+t3^2
        g = puiseux_polynomial_ring_elem(K_p, h)
        @test g.scale == 1

        g = PuiseuxMPolyRingElem(K_p,h,3)
        @test g.scale==3

        h = t1^(2)
        g = puiseux_polynomial_ring_elem(K_p,h)
        @test g.scale == 1
        @test normalize!(g) == false

        h = t1^2*(1 + t1)
        g = puiseux_polynomial_ring_elem(K_p,h)
        @test g.scale == 1

        h = t1*t2^(2)*t3^(3)*(1+t1+t2+t3)
        g = puiseux_polynomial_ring_elem(K_p,h)
        @test g.scale == 1

        h = t1*t2^(2)*t3^(3)*(1+t1+t2+t3)
        g = puiseux_polynomial_ring_elem(K_p,h,skip_normalization=true)
        @test g.scale == 1
        @test normalize!(g) == false

        h = t1^4*t2^2 + t3^6
        g = PuiseuxMPolyRingElem(K_p,h,2)
        @test normalize!(g) == true

        h = t1^(-4) + t2^2
        g = PuiseuxMPolyRingElem(K_p,h,2)
        @test normalize!(g) == true

        h = (1+t1+t2+t3)
        g = puiseux_polynomial_ring_elem(K_p,h,3,skip_normalization=true)
        @test g.scale == 3

        K, _ = laurent_polynomial_ring(QQ, ["t1","t2","t3"])
        Kt, _ = puiseux_polynomial_ring(QQ,["t1","t2","t3"])
        @test Kt.baseRing == K
    end

    @testset "Getters" begin
        K, (t1,t2,t3) = laurent_polynomial_ring(QQ, ["t1","t2","t3"])
        Kp, (tp1,tp2,tp3) = puiseux_polynomial_ring(QQ,["t1","t2","t3"])

	      @test K == base_ring(Kp)
        @test QQ == coefficient_ring(Kp)
        @test ngens(Kp) == 3
        @test gens(Kp) == [tp1,tp2,tp3]
        @test !is_univariate(Kp)

        g = tp1^(1//2)+tp3^(1//3)
        @test elem_type(Kp) == typeof(g)
        @test parent_type(g) == typeof(Kp)
        @test coefficient_ring_type(Kp) == typeof(QQ)
        @test parent(g) == Kp
        @test poly(g) == t1^3 + t3^2
        @test scale(g) == 6

        g = tp1^(2//3)*tp1*tp2^(1//2)*tp3 + tp3^(3//7)*tp1*tp2^(1//2)*tp3 + tp2^(1//2)*tp1*tp2^(1//2)*tp3
        @test parent(g) == Kp
        @test poly(g) == t1^70*t2^21*t3^42 + t1^42*t2^42*t3^42 + t1^42*t2^21*t3^60
        @test scale(g) == 2*3*7

        @test @inferred is_constant(zero(Kp))
        @test @inferred is_constant(Kp(5))
        @test is_constant(3 * tp1^(1//2) * tp1^(-1//2))
        @test !is_constant(tp1^(1//2))
        @test !is_constant(tp1^(-1//3))
        @test !is_constant(1 + tp1)

        @test @inferred(constant_coefficient(zero(Kp))) == 0
        @test @inferred(constant_coefficient(tp1^(-1//2) + 5 + tp2^(1//3))) == 5
        @test constant_coefficient(tp1^(1//2)) == 0

        K, (t,) = puiseux_polynomial_ring(QQ,["t"])
        @test is_univariate(K)
        @test valuation(K(0)) == PosInf()
        @test valuation(t^(-1)) == -1
        @test valuation(t^(-2//3)+t^(-1//2)) == -2//3
        @test_throws ArgumentError valuation(g)
    end

    @testset "Arithmetic" begin
        F, (up,vp,wp) = laurent_polynomial_ring(QQ,["u","v","w"])
        K, (u,v,w) = puiseux_polynomial_ring(QQ,["u","v","w"])
        g = v^(1//3)+u^(1//2)
        h = v^(1//3) + w^(1//3)

        @test monomials(g) == [u^(1//2),v^(1//3)]
        @test monomials(h) == [v^(1//3),w^(1//3)]
        @test collect(coefficients(g)) == [1,1]
        @test collect(coefficients(h)) == [1,1]
        @test collect(exponent_vectors(g)) == [[QQ(1//2),QQ(0),QQ(0)],[QQ(0),QQ(1//3),QQ(0)]]
        @test 0*g == 0
        @test 1*g == g
        @test g*1 == g
        @test 4*g == 4*u^(1//2) + 4*v^(1//3)
        @test g+0 == g
        @test 0+g == g
        @test h+g == u^(1//2) + 2*v^(1//3) + w^(1//3)
        @test h-g == w^(1//3)-u^(1//2)
        @test h*g == u^(1//2)*v^(1//3) + u^(1//2)*w^(1//3) + v^(2//3) + v^(1//3)*w^(1//3)
        @test (g)^3 == u^(3//2) + 3*u*v^(1//3) + 3*u^(1//2)*v^(2//3) + v
        @test (g)^1 == g
        @test (g)^0 == 1
        @test (g)^3 == g^QQ(3)
        @test (g)^1 == g^QQ(1)
        @test (g)^0 == g^QQ(0)

        @test is_unit(u^(1//2))
        @test is_unit(K(2))
        @test is_unit(2*u^(-1//2)*v)
        @test !is_unit(zero(K))
        @test !is_unit(g)

        @test divexact(g, 2) == (1//2)*u^(1//2) + (1//2)*v^(1//3)
        @test divexact(g, QQ(2)) == (1//2)*u^(1//2) + (1//2)*v^(1//3)
        @test divexact(2*g, QQ(2)) == g
        @test_throws ArgumentError divexact(g, QQ(0))
        # TODO: add some more divexact tests

        g = u^(1//2)*v^(2//3) + w^(1//4)
        h = u^(2//3)
        @test g*h == u^(7//6)*v^(2//3)+w^(1//4)*u^(2//3)
        @test (u^(1//2)+u^(-1//2))^2 == u+2+u^(-1)

        g = puiseux_polynomial_ring_elem(K,up*vp*wp*((up^4)*(vp^4)*(wp^4)+1),4,skip_normalization=true)
        @test monomials(g) == [u^(5//4)*v^(5//4)*w^(5//4),u^(1//4)*v^(1//4)*w^(1//4)]
        @test collect(exponent_vectors(g)) == [[5//4,5//4,5//4],[1//4,1//4,1//4]]
        @test poly(g) == up*vp*wp*(up^4*vp^4*wp^4+1)
        @test scale(g) == 4
        @test normalize!(g) == false

        g = puiseux_polynomial_ring_elem(K,up^2*vp^4*wp^8*(1+up^4 - vp^6 + wp^10),4)
        @test scale(g) == 2
        g_c = rescale(g,10)
        @test scale(g_c) == 10
        @test poly(g_c) == up^5*vp^10*wp^20*(1+up^10 - vp^15 + wp^25)
        @test collect(exponent_vectors(g_c)) == [[3//2, 1, 2], [1//2, 5//2, 2], [1//2, 1, 9//2], [1//2, 1, 2]]
        @test normalize!(K(0)) == false
    end

    @testset "Ring properties and delegations" begin
        K, (u,v) = puiseux_polynomial_ring(QQ,["u","v"])
        g = 2*u^(-1//2)*v + 3*v^(1//3)

        @test !is_noetherian(K)
        @test is_noetherian(puiseux_polynomial_ring(QQ, String[])[1])

        @test leading_coefficient(g) in (2, 3)
        @test leading_coefficient(g) == first(coefficients(g))
        @test_throws ArgumentError leading_coefficient(zero(K))

        @test !is_zero_divisor(g)
        @test is_zero_divisor(zero(K))

        cu = canonical_unit(g)
        @test is_unit(cu)
        @test canonical_unit(divexact(g, cu)) == 1
        @test canonical_unit(-g) == -cu

        h = @inferred map_coefficients(c -> 2*c, g)
        @test parent(h) === K
        @test h == 2*g
        @test iszero(map_coefficients(c -> 0*c, g))
        h = map_coefficients(c -> c == 2 ? zero(c) : c, 2*u^(1//2) + v)
        @test h == v
        @test scale(h) == 1

        L, (x, y) = puiseux_polynomial_ring(RealField, ["u","v"])
        h = change_base_ring(RealField, g)
        @test parent(h) === L
        @test h == 2*x^(-1//2)*y + 3*y^(1//3)
        @test map_coefficients(c -> RealField(c), g; parent = L) == h

        ok, q = divides(g*(u^(1//5) - v), u^(1//5) - v)
        @test ok && q == g
        ok, q = divides(u^(1//2), u^(1//3))
        @test ok && q == u^(1//6)
        @test !divides(1 + u^(1//2), 1 + u)[1]
        @test divides(zero(K), g) == (true, zero(K))
        @test divides(zero(K), zero(K)) == (true, zero(K))
        @test !divides(g, zero(K))[1]
    end

    @testset "Conversions" begin
        F, (up,vp,wp) = laurent_polynomial_ring(QQ,["u","v","w"])
        K, (u,v,w) = puiseux_polynomial_ring(QQ,["u","v","w"])
        @test typeof(3) == Int64
        @test K(3) == puiseux_polynomial_ring_elem(K,F(3))
        @test typeof(4//5) == Rational{Int64}
        @test K(4//5) == puiseux_polynomial_ring_elem(K,F(4//5))
        @test K(0) == zero(K)
        @test K(1) == one(K)
    end

    @testset "Ordering" begin
        R, (t,) = puiseux_polynomial_ring(QQ, [:t])
        @test isless(t^2-t,0)
        @test !isless(0,t^2-t)
        K = fraction_field(R)
        @test isless(K(-t),0)
        @test !isless(0,K(-t))
        @test isless(K(t)/K(-1),0)
        @test !isless(0,K(t)/K(-1))
    end

    @testset "Printing" begin
        R, _ = puiseux_polynomial_ring(QQ, ["x", "y"])
        @test sprint(show, R) == "Puiseux polynomial ring in 2 variables over rationals"
        @test sprint(show, R; context = :supercompact => true) == "Puiseux polynomial ring"
    end

    @testset "Conformance tests" begin
        K_p,(tp1,tp2,tp3) = puiseux_polynomial_ring(QQ, ["t1","t2","t3"])
        ConformanceTests.test_Ring_interface(K_p) # basic tests

        K_p,(tp1,) = puiseux_polynomial_ring(QQ, ["t1"])
        ConformanceTests.test_Ring_interface_recursive(K_p) # also tests constructions like mpoly over your ring; if you have this, you don't need the line above
    end
end
