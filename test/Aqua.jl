using Aqua

@testset "Aqua.jl" begin
    Aqua.test_all(
        AbstractAlgebra;
        # 370 ambiguities remain, all of them between a `rand(S::SomeRing, v...)`
        # and one of the `rand(X, args...)` methods of Random/RandomExtensions,
        # which leave `X` untyped. The `v...` there is a sampler spec for the
        # base ring and cannot be typed, so resolving these needs a decision on
        # what e.g. `rand(S, 2, 3)` should mean. See docs/src/rand.md.
        ambiguities=(exclude=[rand],),   # TODO: fix the `rand` ambiguities
        piracies=false                   # TODO: fix piracy
    )
end
