using AbstractAlgebra.Generic: AhoCorasickAutomaton, search, AhoCorasickMatch, aho_corasick_automaton
@testset "Generic.AhoCorasick" begin
    keywords = [[1, 2, 3, 4], [1, 5, 4], [4, 1, 2], [1, 2]]
    aut = aho_corasick_automaton(keywords)
    @test search(aut, [10, 4, 1, 2, 3, 4]) == AhoCorasickMatch(6, 1, [1, 2, 3, 4])
    @test hash(search(aut, [10, 4, 1, 2, 3, 4])) == hash(AhoCorasickMatch(6, 1, [1, 2, 3, 4]))
    @test isnothing(search(aut, []))
    @test search(aut, [1, 5, 4, 1, 1, 1, 4, 4]) == AhoCorasickMatch(3, 2, [1, 5, 4])
    @test search(aut, [1, 2, 3, 1, 4, 1, 2, 1, 4, 1, 2]) == AhoCorasickMatch(7, 3, [4, 1, 2])
    @test search(aut, [2, 1, 2, 3, 1]) == AhoCorasickMatch(3, 4, [1, 2])
    @test isnothing(search(aut, [1, 3, 1, 5, 1, 4, 8]))
    @test isnothing(search(aut, [8, 8, 7, 10, 456]))
    @test search(aut, [4, 1, 5, 4]) == AhoCorasickMatch(4, 2, [1, 5, 4])
    @test isnothing(search(aut, [4, 1, 5, 10]))
    @test !isnothing(AhoCorasickAutomaton(Vector{Int}[]))
end
