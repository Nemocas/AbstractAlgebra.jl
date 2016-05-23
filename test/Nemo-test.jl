include("Groups-test.jl")
include("Rings-test.jl")
include("Fields-test.jl")
include("Benchmark-test.jl")

function test_all()
   test_groups()
   test_rings()
   test_fields()
   test_benchmarks()
end
