module Tmp
    using AbstractAlgebra

    # test @attributes applied to a struct definition, using internal storage
    @attributes mutable struct Foo
        x::Int
        Foo() = new(0)
    end

    # test @attributes applied to a struct typename, using external storage
    mutable struct Bar
        x::Int
        Bar() = new(0)
    end
    @attributes Bar

    mutable struct Quux
        x::Int
        Quux() = new(0)
    end

    struct Singleton
    end

    struct NotSupported
        x::Int
        NotSupported() = new(0)
    end

    # To have the :curly in the expression
    mutable struct FooBar{T}
        x::Int
        FooBar{T}() where T = new(0)
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

@testset "attributes for $T" for T in (Tmp.Foo, Tmp.Bar, Tmp.Quux, Tmp.FooBar{Tmp.Bar}, Tmp.FooBar{Tmp.Quux})

    x = T()

    # test querying attributes when no attribute storage exists
    @test get_attribute(x, :bar) === nothing

    # test setting one attribute
    set_attribute!(x, :bar, 17)
    @test get_attribute(x, :bar) == 17
    @test get_attribute(x, :qux) === nothing

    # verify that `@attributes` does not reset attribute storage for a type which
    # already has attribute storage enabled
    @attributes T
    @test get_attribute(x, :bar) == 17

    # test setting two attributes
    set_attribute!(x, :bar => 42, :qux => "test")
    @test get_attribute(x, :bar) == 42
    @test get_attribute(x, :qux) == "test"

    # test modifying one attribute
    set_attribute!(x, :bar => 43)
    @test get_attribute(x, :bar) == 43
    @test get_attribute(x, :qux) == "test"

    # test get_attribute with default value for new entry
    @test get_attribute(x, :bar2) == nothing
    @test get_attribute(x, :bar2, 42) == 42
    @test get_attribute(x, :bar2) == nothing

    # test get_attribute with default value for pre-existing entry
    set_attribute!(x, :bar3 => 0)
    @test get_attribute(x, :bar3) == 0
    @test get_attribute(x, :bar3, 42) == 0
    @test get_attribute(x, :bar3) == 0

    # test get_attribute with callback for new entry
    @test get_attribute(x, :bar8) == nothing
    @test get_attribute(() -> 42, x, :bar8) == 42
    @test get_attribute(x, :bar8) == nothing

    # test get_attribute with callback for pre-existing entry
    set_attribute!(x, :bar9 => 0)
    @test get_attribute(x, :bar9) == 0
    @test get_attribute(() -> 42, x, :bar9) == 0
    @test get_attribute(x, :bar9) == 0

    # test get_attribute! with default value for new entry
    @test get_attribute(x, :bar4) == nothing
    @test get_attribute!(x, :bar4, 42) == 42
    @test get_attribute(x, :bar4) == 42

    # test get_attribute! with default value for pre-existing entry
    set_attribute!(x, :bar5 => 0)
    @test get_attribute(x, :bar5) == 0
    @test get_attribute!(x, :bar5, 42) == 0
    @test get_attribute(x, :bar5) == 0

    # test get_attribute! with callback for new entry
    @test get_attribute(x, :bar6) == nothing
    @test get_attribute!(() -> 42, x, :bar6) == 42
    @test get_attribute(x, :bar6) == 42

    # test get_attribute! with callback for pre-existing entry
    set_attribute!(x, :bar7 => 0)
    @test get_attribute(x, :bar7) == 0
    @test get_attribute!(() -> 42, x, :bar7) == 0
    @test get_attribute(x, :bar7) == 0

    # test get_attribute! with callback a type
    @test get_attribute(x, :bar8) == nothing
    @test get_attribute!(Vector{Int}, x, :bar8) isa Vector{Int}
    @test get_attribute(x, :bar8) == []

end


# attribute caching


uncached_attr(obj::T) where T = [obj,T]

"""
    cached_attr(obj::T) where T

A cached attribute.
"""
@attr cached_attr(obj::T) where T = [obj,T]

@testset "attribute caching for $T" for T in (Tmp.Foo, Tmp.Bar, Tmp.Quux, Tmp.FooBar{Tmp.Bar}, Tmp.FooBar{Tmp.Quux})

    x = T()

    # check uncached case: multiple calls return equal but not identical results
    y = uncached_attr(x)
    @test y == [x,T]
    @test uncached_attr(x) == y
    @test uncached_attr(x) !== y

    # check cached case: multiple calls return identical results
    y = cached_attr(x)
    @test y == [x,T]
    @test cached_attr(x) == y
    @test cached_attr(x) === y

    # verify docstring is correctly attached
    @test string(@doc cached_attr) ==
        """
        ```
        cached_attr(obj::T) where T
        ```

        A cached attribute.
        """

    # test function location is tracked accurately (this requires that the
    # definition of uncached_attr is before that of cached_attr)
    @test functionloc(uncached_attr)[1] == functionloc(cached_attr)[1]
    @test functionloc(uncached_attr)[2] < functionloc(cached_attr)[2]
end
