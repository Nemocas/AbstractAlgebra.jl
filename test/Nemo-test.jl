include("Groups-test.jl")
include("Rings-test.jl")
include("Fields-test.jl")
include("Benchmark-test.jl")

ZZ = Nemo.ZZ
QQ = Nemo.QQ
FiniteField = Nemo.FiniteField
PadicField = Nemo.PadicField

function test_all()
   test_groups()
   test_rings()
   test_fields()
   test_benchmarks()
end
