@testset "GroupsCore API PermGroup" begin

    include(joinpath(dirname(dirname(pathof(AbstractAlgebra))), "test", "Groups-conformance-tests.jl"))
    @testset "Sym($n)" for n in [1,2,5]
        G = SymmetricGroup(n)
        test_Group_interface(G)
        test_GroupElem_interface(rand(G, 2)...)
    end
end
