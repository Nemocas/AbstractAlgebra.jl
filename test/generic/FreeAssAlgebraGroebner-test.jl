include("AhoCorasick-test.jl")
using AbstractAlgebra.Generic: AhoCorasickAutomaton
import AbstractAlgebra.Generic: normal_form_weak
@testset "Generic.FreeAssAlgebra.groebner" begin
 
    R, (x, y, u, v, t, s) = FreeAssociativeAlgebra(GF(2), ["x", "y", "u", "v", "t", "s"])
    g = AbstractAlgebra.groebner_basis([u*(x*y)^3 + u*(x*y)^2 + u + v,
                                       (y*x)^3*t + (y*x)^2*t + t + s])
    @test length(g) >= 5
 
    # Example 6.1 Kreuzer & Xiu
    R, (a, b) = FreeAssociativeAlgebra(QQ, ["a", "b"])
 
    g = AbstractAlgebra.groebner_basis([a^2 - 1, b^3 - 1, (a*b*a*b^2)^2 - 1])
    AbstractAlgebra.interreduce!(g)
    @test length(g) == 5
 
    g = AbstractAlgebra.groebner_basis([a^2 - 1, b^3 - 1, (a*b*a*b*a*b^2)^2 - 1])
    AbstractAlgebra.interreduce!(g)
    @test length(g) == 15
end

@testset "Generic.FreeAssociativeAlgebra.groebner.normal_form" begin
   R, (x, y, u, v, t, s) = FreeAssociativeAlgebra(QQ, ["x", "y", "u", "v", "t", "s"]) # x > y > ... > s
   ideal_generators = [x*y, u*y*t, s*t - t*s, x*y*y + x*x*y - one(R)]
   aut = AhoCorasickAutomaton([g_i.exps[1] for g_i in ideal_generators])
   @test normal_form(x*y, ideal_generators, aut) == zero(R)
   @test normal_form(u*y*t, ideal_generators, aut) == zero(R)
   @test normal_form(s*t - t*s, ideal_generators, aut) == zero(R)
   @test x*y*y + x*x*y - one(R) in ideal_generators
   @test normal_form(x*y*t + y, ideal_generators, aut) == y
   @test normal_form(one(R), ideal_generators, aut) == one(R)
   @test normal_form(v*y, ideal_generators, aut) == v*y
   @test normal_form(x*y*v*v + t*s*x*y*v + y*s*t - y*t*s + v*x*y*y*s + v*x*x*y*s, ideal_generators, aut) == zero(R)
   @test normal_form_weak(x*y*v*v + t*s*x*y*v + y*s*t - y*t*s + v*x*y*y*s + v*x*x*y*s, ideal_generators) == zero(R)
   @test normal_form_weak(x*y + u*y, ideal_generators) <= x*y + u*y

#   @test gb_divides_leftmost((x*y*u*v*t).exps[1], aut) == (true, [], [3, 4, 5], 1)
#   @test gb_divides_leftmost((x*u*y*t*(s*t - t*s)).exps[1], aut) == (true, [1], [5, 6], 2)
#   @test gb_divides_leftmost((x*s*t).exps[1], aut) == (false, [], [], -1)
end

