include("generic/Map-test.jl")
include("generic/MapCache-test.jl")
include("generic/MapWithInverse-test.jl")

function test_maps()
   test_gen_map()
   test_gen_map_cache()
   test_gen_map_with_inverse()

end
