module Test

using Base.Test, Rings, Fields

export test_all

include("Fraction-test.jl")
include("FiniteFields-test.jl")

function test_all()
   test_fraction()
   test_ffield()
end

end # module
