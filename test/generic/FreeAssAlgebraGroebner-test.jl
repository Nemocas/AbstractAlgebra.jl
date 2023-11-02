include("AhoCorasick-test.jl")
using AbstractAlgebra.Generic: AhoCorasickAutomaton
import AbstractAlgebra.Generic: normal_form_weak
@testset "Generic.FreeAssAlgebra.groebner" begin
 
    R, (x, y, u, v, t, s) = free_associative_algebra(GF(2), ["x", "y", "u", "v", "t", "s"])
    g = AbstractAlgebra.groebner_basis([u*(x*y)^3 + u*(x*y)^2 + u + v,
                                       (y*x)^3*t + (y*x)^2*t + t + s])
    @test length(g) >= 5
 
    # Example 6.1 Kreuzer & Xiu
    R, (a, b) = free_associative_algebra(QQ, ["a", "b"])
 
    g = AbstractAlgebra.groebner_basis([a^2 - 1, b^3 - 1, (a*b*a*b^2)^2 - 1])
    AbstractAlgebra.interreduce!(g)
    @test length(g) == 5
 
    g = AbstractAlgebra.groebner_basis([a^2 - 1, b^3 - 1, (a*b*a*b*a*b^2)^2 - 1])
    AbstractAlgebra.interreduce!(g)
    @test length(g) == 15

    g2 = AbstractAlgebra.groebner_basis([a^2 - 1, b^3 - 1, (a*b*a*b*a*b^2)^2 - 1], typemax(Int), true)
    @test all(u ->iszero(normal_form(u, g2)), g) # make sure removing redundant obstructions in the computation does not change the groebner basis
    @test all(u ->iszero(normal_form(u, g)), g2) 

end

@testset "Generic.free_associative_algebra.groebner.normal_form" begin
   R, (x, y, u, v, t, s) = free_associative_algebra(QQ, ["x", "y", "u", "v", "t", "s"]) # x > y > ... > s
   ideal_generators = [x*y, u*y*t, s*t - t*s, x*y*y + x*x*y - one(R)]
   aut = AhoCorasickAutomaton([g_i.exps[1] for g_i in ideal_generators])
   @test normal_form(x*y, ideal_generators, aut) == zero(R)
   @test normal_form(u*y*t, ideal_generators, aut) == zero(R)
   @test normal_form(s*t - t*s, ideal_generators, aut) == zero(R)
   @test x*y*y + x*x*y - one(R) in ideal_generators
   @test normal_form(x*y*t + y, ideal_generators, aut) == y
   @test normal_form(one(R), ideal_generators, aut) == one(R)
   @test normal_form(v*y, ideal_generators, aut) == v*y
   @test normal_form(x*y, ideal_generators) == zero(R)
   @test normal_form(u*y*t, ideal_generators) == zero(R)
   @test normal_form(s*t - t*s, ideal_generators) == zero(R)
   @test normal_form(x*y*t + y, ideal_generators) == y
   @test normal_form(one(R), ideal_generators) == one(R)
   @test normal_form(v*y, ideal_generators) == v*y
   @test normal_form(x*y*v*v + t*s*x*y*v + y*s*t - y*t*s + v*x*y*y*s + v*x*x*y*s, ideal_generators, aut) == zero(R)
   @test normal_form_weak(x*y*v*v + t*s*x*y*v + y*s*t - y*t*s + v*x*y*y*s + v*x*x*y*s, ideal_generators) == zero(R)
   @test normal_form_weak(x*y + u*y, ideal_generators) <= x*y + u*y

#   @test gb_divides_leftmost((x*y*u*v*t).exps[1], aut) == (true, [], [3, 4, 5], 1)
#   @test gb_divides_leftmost((x*u*y*t*(s*t - t*s)).exps[1], aut) == (true, [1], [5, 6], 2)
#   @test gb_divides_leftmost((x*s*t).exps[1], aut) == (false, [], [], -1)
end

@testset "Generic.free_associative_algebra.groebner.overlaps_and_obstructions" begin
    w1 = [1, 1, 2, 1, 3]
    w2 = [2, 1, 3, 4, 3, 4]
    w3 = [1, 3, 4]
    w4 = [1, 1, 2, 1]
    w5 = [5, 1, 3, 4, 4, 2, 1]
    @test AbstractAlgebra.Generic.check_left_overlap(w2, w1, 3)
    @test !AbstractAlgebra.Generic.check_left_overlap(w1, w2, 3)
    @test AbstractAlgebra.Generic.check_center_overlap(w3, w2, 2)
    @test AbstractAlgebra.Generic.check_center_overlap(w4, w1, 1)
    @test !AbstractAlgebra.Generic.check_left_overlap(w4, w1, 1)
    R, (x, y, z) = free_associative_algebra(QQ, ["x", "y", "z"])
    poly1 = x*y*x*x*z
    poly2 = y*x*x*z*y*y
    poly3 = x*y*x*y
    poly4 = y*x*y*z
    poly5 = x*y
    poly6 = x*y*x*y*z*x*y
    lw1 = AbstractAlgebra.Generic._leading_word(poly1)
    lw2 = AbstractAlgebra.Generic._leading_word(poly2)
    lw3 = AbstractAlgebra.Generic._leading_word(poly3)
    lw4 = AbstractAlgebra.Generic._leading_word(poly4)
    lw5 = AbstractAlgebra.Generic._leading_word(poly5)
    lw6 = AbstractAlgebra.Generic._leading_word(poly6)
    ot = AbstractAlgebra.Generic.ObstructionTriple{Rational{BigInt}}( poly1, poly2, [], [2, 2], [1], [], 1, 2)

    @test AbstractAlgebra.Generic.has_overlap(ot)
    @test length(AbstractAlgebra.Generic.left_obstructions(lw2, lw1)) == 1
    @test isempty(AbstractAlgebra.Generic.left_obstructions(lw1, lw2))
    @test length(AbstractAlgebra.Generic.right_obstructions(lw3, lw4)) == 2
    @test isempty(AbstractAlgebra.Generic.left_obstructions(lw3, lw4))
    @test length(AbstractAlgebra.Generic.center_obstructions(lw5, lw6)) == 3
    @test length(AbstractAlgebra.Generic.get_obstructions([poly1, poly2])) == 2
end
