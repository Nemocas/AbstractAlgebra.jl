using Random: MersenneTwister, randsubseq

const rng = MersenneTwister()

include("Groups-test.jl")
include("NCRings-test.jl")
include("Rings-test.jl")
include("Fields-test.jl")
include("Modules-test.jl")
include("Maps-test.jl")
include("Benchmark-test.jl")
