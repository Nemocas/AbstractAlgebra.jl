include("generic/Map-test.jl")
include("generic/MapCache-test.jl")

function test_maps()
   test_gen_map()
   test_gen_map_cache()
end
