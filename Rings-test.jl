module Test

using Base.Test, Rings

export test_all

include("Poly-test.jl")

function test_all()
   test_poly()
end

end # module
