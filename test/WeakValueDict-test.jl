WeakValueDict = AbstractAlgebra.WeakValueDict

@testset "WeakValueDict" begin
    A = [1]
    B = [2]
    C = [3]

    # construction
    wkd = WeakValueDict()
    wkd[2] = A
    wkd[3] = B
    wkd[4] = C
    dd = convert(Dict{Any,Any},wkd)
    @test WeakValueDict(dd) == wkd
    @test convert(WeakValueDict{Any, Any}, dd) == wkd
    @test isa(WeakValueDict(dd), WeakValueDict{Any,Any})
    @test WeakValueDict(2=>A, 3=>B, 4=>C) == wkd
    @test isa(WeakValueDict(2=>A, 3=>B, 4=>C), WeakValueDict{Int, Vector{Int}})
    @test WeakValueDict((i+1)=>a for (i,a) in enumerate([A,B,C]) ) == wkd
    @test WeakValueDict([(2,A), (3,B), (4,C)]) == wkd
    @test WeakValueDict{Int, typeof(A)}(Pair(2,A), Pair(3,B), Pair(4,C)) == wkd
    @test WeakValueDict(Pair(2,A), Pair(3,B), Pair(4,C)) == wkd
    D = [[4.0]]
    @test WeakValueDict(Pair(2,A), Pair(3,B), Pair(4.0,D)) isa WeakValueDict{Any, Any}
    @test isa(WeakValueDict(Pair(2,A), Pair(3.0,B), Pair(4,C)), WeakValueDict{Any, Vector{Int}})
    @test copy(wkd) == wkd

    @test length(wkd) == 3
    @test !isempty(wkd)
    res = pop!(wkd, 3)
    @test res == B
    @test length(wkd) == 2
    res = pop!(wkd, 3, C)
    @test res == C
    @test 3 ∉ keys(wkd)
    @test B ∉ values(wkd)
    @test length(wkd) == 2
    @test !isempty(wkd)
    wkd = filter!( p -> p.first != 2, wkd)
    @test 3 ∉ keys(wkd)
    @test B ∉ values(wkd)
    @test length(wkd) == 1
    @test WeakValueDict(Pair(4, C)) == wkd
    @test !isempty(wkd)

    wkd = empty!(wkd)
    @test wkd == empty(wkd)
    @test typeof(wkd) == typeof(empty(wkd))
    @test length(wkd) == 0
    @test isempty(wkd)
    @test isa(wkd, WeakValueDict)

    @test_throws ArgumentError WeakValueDict([1, 2, 3])

    wkd = WeakValueDict(1=>A)
    @test delete!(wkd, 1) == empty(wkd)
    @test delete!(wkd, 1) === wkd

    GC.@preserve A B C D nothing
end

