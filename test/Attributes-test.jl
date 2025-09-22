using REPL # needed due to https://github.com/JuliaLang/julia/issues/53349

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

    # applying @attributes to a singleton type definition is supported but does nothing
    @attributes struct AnotherSingleton
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


    @attributes mutable struct Container{T}
        x::T
        Container(x::T) where T = new{T}(x)
    end
end

# test @attributes applied to a struct typename in another module
@attributes Tmp.Quux
@attributes Tmp.FooBar{Tmp.Quux}

@testset "is_attribute_storing_type" begin
  @test is_attribute_storing_type(Tmp.Foo)
  @test is_attribute_storing_type(Tmp.Bar)
  @test is_attribute_storing_type(Tmp.Quux)

  @test !is_attribute_storing_type(Tmp.FooBar)
  @test !is_attribute_storing_type(Tmp.FooBar{Tmp.Foo})
  @test is_attribute_storing_type(Tmp.FooBar{Tmp.Bar})
  @test is_attribute_storing_type(Tmp.FooBar{Tmp.Quux})

  @test is_attribute_storing_type(Tmp.Container)
  @test is_attribute_storing_type(Tmp.Container{Tmp.Bar})
  @test is_attribute_storing_type(Tmp.Container{Tmp.Quux})
end

# applying @attributes to a singleton typename is supported but does nothing
@attributes Tmp.Singleton

@testset "@attributes input validation" begin

    @test_throws Exception @attributes Int
    @test_throws Exception @attributes Tmp.NotSupported

end

@testset "attributes for $T" for T in (Tmp.Foo, Tmp.Bar, Tmp.Quux, Tmp.Singleton, Tmp.AnotherSingleton, Tmp.FooBar{Tmp.Bar}, Tmp.FooBar{Tmp.Quux})

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

    # test get_attribute with symbol entries
    set_attribute!(x, :bar10 => :somesymbol)
    @test get_attribute(x, :bar10) == :somesymbol
    @test get_attribute(x, :bar11) == nothing
    @test get_attribute(x, :bar11, :defaultsymbol) == :defaultsymbol

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

    # test get_attribute with symbol entries
    set_attribute!(x, :bar12 => :somesymbol)
    @test get_attribute(x, :bar12) == :somesymbol
    @test get_attribute!(x, :bar12, :defaultsymbol) == :somesymbol
    @test get_attribute(x, :bar12) == :somesymbol
    @test get_attribute!(x, :bar13, :defaultsymbol) == :defaultsymbol
    @test get_attribute(x, :bar13) == :defaultsymbol
end


# attribute caching


uncached_attr(obj::T) where T = (obj,T,[])

"""
    cached_attr(obj::T) where T

A cached attribute.
"""
@attr Tuple{T,DataType,Vector{Any}} cached_attr(obj::T) where T = (obj,T,[])

# cached attribute with bad return type specification
@attr Any cached_attr2(obj::T) where T = (obj,T,[])

# cached attribute with return type specification depending on the input type
my_derived_type(::Type{Tmp.Container{T}}) where T = T
@attr my_derived_type(T) cached_attr3(obj::T) where T <: Tmp.Container = obj.x

@attr Tuple{T,DataType,Bool} cached_attr_with_kwarg1(obj::T; some_kwarg::Bool) where T = (obj,T,some_kwarg)
@attr Tuple{T,DataType,Bool} cached_attr_with_kwarg2(obj::T; some_kwarg::Bool=true) where T = (obj,T,some_kwarg)

@testset "attribute caching for $T" for T in (Tmp.Foo, Tmp.Bar, Tmp.Quux, Tmp.FooBar{Tmp.Bar}, Tmp.FooBar{Tmp.Quux})

    x = T()

    # check uncached case: multiple calls return equal but not identical results
    y = @inferred uncached_attr(x)
    @test y == (x,T,[])
    @test uncached_attr(x) == y
    @test uncached_attr(x) !== y

    # check cached case: multiple calls return identical results
    y = @inferred cached_attr(x)
    @test y == (x,T,[])
    @test cached_attr(x) == y
    @test cached_attr(x) === y

    # check cached with bad type specification is not inferring return types correctly
    @test_throws ErrorException @inferred cached_attr2(x)

    # check when return type is derived via a function from input type
    z = Tmp.Container(x)
    y = @inferred cached_attr3(z)
    @test y === x

    # check case of ignored keyword arguments
    x = T()
    y = @inferred cached_attr_with_kwarg1(x; some_kwarg=true)
    @test y == (x,T,true)
    @test cached_attr_with_kwarg1(x; some_kwarg=true) === y
    @test cached_attr_with_kwarg1(x; some_kwarg=false) === y

    x = T()
    y = @inferred cached_attr_with_kwarg1(x; some_kwarg=false)
    @test y == (x,T,false)
    @test cached_attr_with_kwarg1(x; some_kwarg=true) === y
    @test cached_attr_with_kwarg1(x; some_kwarg=false) === y

    x = T()
    y = @inferred cached_attr_with_kwarg2(x; some_kwarg=true)
    @test y == (x,T,true)
    @test cached_attr_with_kwarg2(x; some_kwarg=true) === y
    @test cached_attr_with_kwarg2(x; some_kwarg=false) === y
    @test cached_attr_with_kwarg2(x) === y

    x = T()
    y = @inferred cached_attr_with_kwarg2(x; some_kwarg=false)
    @test y == (x,T,false)
    @test cached_attr_with_kwarg2(x; some_kwarg=true) === y
    @test cached_attr_with_kwarg2(x; some_kwarg=false) === y
    @test cached_attr_with_kwarg2(x) === y

    x = T()
    y = @inferred cached_attr_with_kwarg2(x)
    @test y == (x,T,true)
    @test cached_attr_with_kwarg2(x; some_kwarg=true) === y
    @test cached_attr_with_kwarg2(x; some_kwarg=false) === y
    @test cached_attr_with_kwarg2(x) === y

    # verify docstring is correctly attached
    if VERSION >= v"1.12.0-DEV.1223"
        @test string(@doc cached_attr) ==
            """
            ```julia
            cached_attr(obj::T) where T
            ```

            A cached attribute.
            """
    else
        @test string(@doc cached_attr) ==
            """
            ```
            cached_attr(obj::T) where T
            ```

            A cached attribute.
            """
    end

    # test function location is tracked accurately (this requires that the
    # definition of uncached_attr is before that of cached_attr)
    @test functionloc(uncached_attr)[1] == functionloc(cached_attr)[1]
    @test functionloc(uncached_attr)[2] < functionloc(cached_attr)[2]
end

@testset "@attr error handling" begin
    # wrong number of arguments
    @test_throws ArgumentError @macroexpand @attr Any foo() = 1
    @test_throws ArgumentError @macroexpand @attr Any foo(x::Int, y::Int) = 1
    @test_throws ArgumentError @macroexpand @attr Int foo() = 1
    @test_throws ArgumentError @macroexpand @attr Int foo(x::Int, y::Int) = 1
    @test_throws ArgumentError @macroexpand @attr Any foo(; some_kwarg::Bool) = 1
    @test_throws ArgumentError @macroexpand @attr Any foo(; some_kwarg::Bool=true) = 1
    @test_throws ArgumentError @macroexpand @attr Any foo(x::Int, y::Int; some_kwarg::Bool) = 1
    @test_throws ArgumentError @macroexpand @attr Any foo(x::Int, y::Int; some_kwarg::Bool=true) = 1
    @test_throws MethodError @macroexpand @attr Int foo(x::Int) = 1 Any
    @test_throws MethodError @macroexpand @attr Int Int Int

    # wrong kind of arguments
    #@test_throws ArgumentError @macroexpand @attr Int Int
    #@test_throws ArgumentError @macroexpand @attr foo(x::Int) = 1 Int
end
