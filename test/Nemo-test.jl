include("Groups-test.jl")
include("Rings-test.jl")
include("Fields-test.jl")
include("Benchmark-test.jl")

function test_all()
   test_fmpq_abs_series()
   test_gen_abs_series()
   test_fmpz_abs_series()

   test_groups()
   test_rings()
   test_fields()
   test_benchmarks()
end
