import GroupsCore

@testset "GroupsCore API PermGroup" begin

    include(joinpath(pathof(GroupsCore), "..", "..", "test", "conformance_test.jl"))
    let G = SymmetricGroup(5)
        test_Group_interface(G)
        test_GroupElement_interface(rand(G, 2)...)
    end
end
