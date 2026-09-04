module Bla
    using AbstractAlgebra
    const sfsdf = 1

    is_foo(x::Int) = x
    is_bar(x::Int) = x

    export is_foo

    @alias isfoo is_foo
    @alias isbar is_bar

    # verify idempotence
    @alias isbar is_bar

    # verify symmetric idempotence
    @alias is_bar isbar
end

@testset "Aliases" begin

    @test Bla.isfoo === Bla.is_foo
    @test :is_foo in names(Bla)
    @test :isfoo in names(Bla)

    @test Bla.isbar === Bla.is_bar
    @test !(:is_bar in names(Bla))
    @test !(:isbar in names(Bla))

    # verify docstring is correctly attached
    if VERSION >= v"1.12.0-DEV.1223"
      @test string(@doc Bla.isfoo) ==
          """
          ```julia
          isfoo
          ```

          Alias for `is_foo`.
          """
    else
      @test string(@doc Bla.isfoo) ==
          """
          ```
          isfoo
          ```

          Alias for `is_foo`.
          """
    end
end
