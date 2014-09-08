module Test

using Base.Test, Rings, Fields

export test_all

include("Fraction-test.jl")

function test_all()
   test_fraction()
end

end # module
