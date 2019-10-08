using Random: MersenneTwister

const rng = MersenneTwister()

include("Maps-test.jl")
include("Groups-test.jl")
include("NCRings-test.jl")
include("Rings-test.jl")
include("Fields-test.jl")
include("Modules-test.jl")
include("Benchmark-test.jl")

@testset "parent_type/elem_type..." begin
   @test_throws MethodError parent_type('c')
   @test_throws MethodError parent_type(Char)
   @test parent_type(big(1)) == AbstractAlgebra.Integers{BigInt}

   @test_throws MethodError elem_type('c')
   @test_throws MethodError elem_type(Char)
end
