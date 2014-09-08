module Test

using Base.Test, Rings, Fields

export test_all

include("Fraction-test.jl")
include("FiniteFields-test.jl")
include("Padics2-test.jl")

function test_all()
   test_fraction()
   test_ffield()
   test_padics2()
end

end # module
