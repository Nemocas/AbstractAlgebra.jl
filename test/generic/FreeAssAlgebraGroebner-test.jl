#@testset "Generic.FreeAssAlgebra.groebner" begin
#
#   R, (x, y, u, v, t, s) = FreeAssociativeAlgebra(GF(2), ["x", "y", "u", "v", "t", "s"])
#   g = AbstractAlgebra.groebner_basis([u*(x*y)^3 + u*(x*y)^2 + u + v,
#                                      (y*x)^3*t + (y*x)^2*t + t + s])
#   @test length(g) >= 5
#
#   # Example 6.1 Kreuzer & Xiu
#   R, (a, b) = FreeAssociativeAlgebra(QQ, ["a", "b"])
#
#   g = AbstractAlgebra.groebner_basis([a^2 - 1, b^3 - 1, (a*b*a*b^2)^2 - 1])
#   AbstractAlgebra.interreduce!(g)
#   @test length(g) == 5
#
#   g = AbstractAlgebra.groebner_basis([a^2 - 1, b^3 - 1, (a*b*a*b*a*b^2)^2 - 1])
#   AbstractAlgebra.interreduce!(g)
#   @test length(g) == 15
#end

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
   @test gb_divides_leftmost_aho_corasick((x*y*u*v*t).exps[1], aut) == (true, [], [3, 4, 5], 1)
   @test gb_divides_leftmost_aho_corasick((x*u*y*t*(s*t - t*s)).exps[1], aut) == (true, [1], [5, 6], 2)
   @test gb_divides_leftmost_aho_corasick((x*s*t).exps[1], aut) == (false, [], [], -1)
end

include("/home/julien/uni/mathe/Promotion/ncmpoly/graph_to_ass_alg.jl")
using DataStructures

@testset "Generic.FreeAssociativeAlgebra.groebner_basis" begin
#    M = [0 1 0 0 1; 1 0 1 1 1; 0 1 0 1 1; 0 1 1 0 1; 1 1 1 1 0]
    M = [0 0 0 0 0 1 1; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1; 0 0 0 0 0 0 1; 0 0 0 0 0 0 1; 1 1 0 0 0 0 0; 1 0 1 1 1 0 0]
    A, u, r = free_alg_for_graph(M)
    g = copy(r)
    checked_obstructions = 0
    nonzero_reductions = 0
    reduction_bound = 5000
    # compute the aho corasick automaton
    # to make normal form computation more efficient
    aut = AhoCorasickAutomaton([g_i.exps[1] for g_i in g])
 
    suffix_match_vectors = Vector{Vector{Int}}(undef, length(g))
       for i in 1:length(g)
           suffix_match_vectors[i] = calc_suffix_match_vector(g[i].exps[1])
       end
    # step 1
    obstruction_queue = get_obstructions(g) 
    while !isempty(obstruction_queue)
 #       #println(length(obstruction_queue))
       obstruction = dequeue!(obstruction_queue)
       # step3 
       S = s_polynomial(obstruction)
       Sp = normal_form(S, g, aut) # or normal_form_weak
       Sp2 = normal_form(S, g, suffix_match_vectors)
       @test Sp == Sp2
       if iszero(Sp)
           #println("is zero")
          continue
       end
       #println(Sp)
       nonzero_reductions += 1
       # step4
       push!(g, Sp)
       aut = AhoCorasickAutomaton([g_i.exps[1] for g_i in g])
       push!(suffix_match_vectors, calc_suffix_match_vector(Sp2.exps[1]))
 
       if nonzero_reductions >= reduction_bound
               return g
       end
       add_obstructions!(obstruction_queue, g)
    end
    return g
end
