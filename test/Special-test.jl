import AbstractAlgebra: get_special, get_special!, set_special, hasspecial

module Tmp
    import AbstractAlgebra: @declare_other
    mutable struct Foo
        x::Int
        @declare_other
        Foo(x::Int) = new(x)
    end
end

@testset "declare_other and get_special" begin
    x = Tmp.Foo(1)

    @test !hasspecial(x)[1]
    @test get_special(x, :bar) === nothing
    @test get_special(x, :qux) === nothing
    @test get_special(x, :frob) === nothing

    set_special(x, :bar => 42, :qux => "test")
    @test get_special(x, :bar) == 42
    @test get_special(x, :qux) == "test"

    set_special(x, :bar => 43)
    @test get_special(x, :bar) == 43
    @test get_special(x, :qux) == "test"

    # test get_special! with default value for pre-existing entry
    @test get_special!(x, :bar, 44) == 43
    @test get_special(x, :bar) == 43

    # test get_special! with default value for new entry
    @test get_special(x, :frob) == nothing
    @test get_special!(x, :frob, 45) == 45
    @test get_special(x, :frob) == 45

    # test get_special! with callback for pre-existing entry
    @test get_special!(x, :bar, 44) == 43
    @test get_special(x, :bar) == 43

    # test get_special! with callback value for new entry
    @test get_special(x, :blam) == nothing
    @test 45 == get_special!(x, :blam) do
                    return 45
                end
    @test get_special(x, :blam) == 45
end
