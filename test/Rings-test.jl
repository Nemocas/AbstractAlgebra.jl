include("julia/Integers-test.jl")

include("generic/Poly-test.jl")
include("generic/Residue-test.jl")
include("generic/ResidueField-test.jl")
include("generic/RelSeries-test.jl")
include("generic/AbsSeries-test.jl")
include("generic/LaurentSeries-test.jl")
include("generic/PuiseuxSeries-test.jl")
include("generic/Matrix-test.jl")
include("generic/MPoly-test.jl")

@testset "Generic.Rings.broadcast..." begin
   F = GF(3)
   @test F(2) .* [F(1), F(2)] == [F(2), F(1)]
end
