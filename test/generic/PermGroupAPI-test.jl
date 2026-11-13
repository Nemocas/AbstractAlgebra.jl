@testset "Groups API PermGroup" begin
    @testset "Sym($n)" for n in [1,2,5,10]
        G = SymmetricGroup(n)
        ConformanceTests.test_Group_interface(G)
        ConformanceTests.test_GroupElem_interface(rand(G, 2)...)
    end
end
