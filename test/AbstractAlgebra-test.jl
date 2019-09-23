using Random: MersenneTwister

const rng = MersenneTwister()

include("Maps-test.jl")
include("Groups-test.jl")
include("NCRings-test.jl")
include("Rings-test.jl")
include("Fields-test.jl")
include("Modules-test.jl")
include("Benchmark-test.jl")

function test_all()
   test_maps()
   test_groups()
   test_ncrings()
   test_rings()
   test_fields()
   test_modules()
   test_benchmarks()
end
