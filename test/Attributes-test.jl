module Tmp
    using AbstractAlgebra

    # test @attributes applied to a struct definition, using internal storage
    @attributes mutable struct Foo
        x::Int
        Foo(x::Int) = new(x)
    end

    # test @attributes applied to a struct typename, using external storage
    mutable struct Bar
        x::Int
    end
    @attributes Bar

    mutable struct Quux
        x::Int
    end

    struct Singleton
    end

    struct NotSupported
        x::Int
    end

    # To have the :curly in the expression
    mutable struct FooBar{T}
        x::Int
    end
    @attributes Tmp.FooBar{Bar}
end

# test @attributes applied to a struct typename in another module
@attributes Tmp.Quux
@attributes Tmp.FooBar{Tmp.Quux}

@testset "@attributes input validation" begin

    @test_throws Exception @attributes Int
    @test_throws Exception @attributes Tmp.NotSupported

end

@testset "attributes for singletons" begin

    T = Tmp.Singleton

    # test querying attributes when no attribute storage exists
    x = T()
    @test get_attribute(x, :bar) === nothing

    # test setting one attribute
    x = T()
    set_attribute!(x, :bar, 17)
    @test get_attribute(x, :bar) == 17
    @test get_attribute(x, :qux) === nothing

    # verify that `@attributes` does not reset attribute storage for a type which
    # already has attribute storage enabled
    @attributes T
    @test get_attribute(x, :bar) == 17

    # test setting two attributes
    x = T()
    set_attribute!(x, :bar => 42, :qux => "test")
    @test get_attribute(x, :bar) == 42
    @test get_attribute(x, :qux) == "test"

    # test modifying one attribute
    set_attribute!(x, :bar => 43)
    @test get_attribute(x, :bar) == 43
    @test get_attribute(x, :qux) == "test"

    # test get_attribute with default value for new entry
    x = T()
    @test get_attribute(x, :bar2) == nothing
    @test get_attribute(x, :bar2, 42) == 42
    @test get_attribute(x, :bar2) == nothing

    # test get_attribute with default value for pre-existing entry
    x = T()
    set_attribute!(x, :bar3 => 0)
    @test get_attribute(x, :bar3) == 0
    @test get_attribute(x, :bar3, 42) == 0
    @test get_attribute(x, :bar3) == 0

    # test get_attribute! with default value for new entry
    x = T()
    @test get_attribute(x, :bar4) == nothing
    @test get_attribute!(x, :bar4, 42) == 42
    @test get_attribute(x, :bar4) == 42

    # test get_attribute! with default value for pre-existing entry
    x = T()
    set_attribute!(x, :bar5 => 0)
    @test get_attribute(x, :bar5) == 0
    @test get_attribute!(x, :bar5, 42) == 0
    @test get_attribute(x, :bar5) == 0

    # test get_attribute! with callback for new entry
    x = T()
    @test get_attribute(x, :bar6) == nothing
    @test get_attribute!(() -> 42, x, :bar6) == 42
    @test get_attribute(x, :bar6) == 42

    # test get_attribute! with callback for pre-existing entry
    x = T()
    set_attribute!(x, :bar7 => 0)
    @test get_attribute(x, :bar7) == 0
    @test get_attribute!(() -> 42, x, :bar7) == 0
    @test get_attribute(x, :bar7) == 0

end

@testset "attributes for $T" for T in (Tmp.Foo, Tmp.Bar, Tmp.Quux, Tmp.FooBar{Tmp.Bar}, Tmp.FooBar{Tmp.Quux})

    # test querying attributes when no attribute storage exists
    x = T(1)
    @test get_attribute(x, :bar) === nothing

    # test setting one attribute
    x = T(1)
    set_attribute!(x, :bar, 17)
    @test get_attribute(x, :bar) == 17
    @test get_attribute(x, :qux) === nothing

    # verify that `@attributes` does not reset attribute storage for a type which
    # already has attribute storage enabled
    @attributes T
    @test get_attribute(x, :bar) == 17

    # test setting two attributes
    x = T(1)
    set_attribute!(x, :bar => 42, :qux => "test")
    @test get_attribute(x, :bar) == 42
    @test get_attribute(x, :qux) == "test"

    # test modifying one attribute
    set_attribute!(x, :bar => 43)
    @test get_attribute(x, :bar) == 43
    @test get_attribute(x, :qux) == "test"

    # test get_attribute with default value for new entry
    x = T(1)
    @test get_attribute(x, :bar) == nothing
    @test get_attribute(x, :bar, 42) == 42
    @test get_attribute(x, :bar) == nothing

    # test get_attribute with default value for pre-existing entry
    x = T(1)
    set_attribute!(x, :bar => 0)
    @test get_attribute(x, :bar) == 0
    @test get_attribute(x, :bar, 42) == 0
    @test get_attribute(x, :bar) == 0

    # test get_attribute! with default value for new entry
    x = T(1)
    @test get_attribute(x, :bar) == nothing
    @test get_attribute!(x, :bar, 42) == 42
    @test get_attribute(x, :bar) == 42

    # test get_attribute! with default value for pre-existing entry
    x = T(1)
    set_attribute!(x, :bar => 0)
    @test get_attribute(x, :bar) == 0
    @test get_attribute!(x, :bar, 42) == 0
    @test get_attribute(x, :bar) == 0

    # test get_attribute! with callback for new entry
    x = T(1)
    @test get_attribute(x, :bar) == nothing
    @test get_attribute!(() -> 42, x, :bar) == 42
    @test get_attribute(x, :bar) == 42

    # test get_attribute! with callback for pre-existing entry
    x = T(1)
    set_attribute!(x, :bar => 0)
    @test get_attribute(x, :bar) == 0
    @test get_attribute!(() -> 42, x, :bar) == 0
    @test get_attribute(x, :bar) == 0
end
