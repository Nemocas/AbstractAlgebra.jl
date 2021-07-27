include("julia/Integers-test.jl")

include("generic/Poly-test.jl")
include("generic/SparsePoly-test.jl")
include("generic/LaurentPoly-test.jl")
include("generic/Residue-test.jl")
include("generic/ResidueField-test.jl")
include("generic/RelSeries-test.jl")
include("generic/AbsSeries-test.jl")
include("generic/LaurentSeries-test.jl")
include("generic/PuiseuxSeries-test.jl")
include("generic/AbsMSeries-test.jl")
include("generic/Matrix-test.jl")
include("generic/MPoly-test.jl")
include("algorithms/MPolyFactor-test.jl")
include("algorithms/DensePoly-test.jl")

@testset "Generic.Rings.broadcast" begin
   F = GF(3)
   @test F(2) .* [F(1), F(2)] == [F(2), F(1)]
end

@testset "Generic.Rings.elem/parent_type" begin
   @test_throws MethodError parent_type('c')
   @test_throws MethodError parent_type(Char)
   @test parent_type(big(1)) == AbstractAlgebra.Integers{BigInt}

   @test_throws MethodError elem_type('c')
   @test_throws MethodError elem_type(Char)
end

include("julia/Integers-test.jl")

# very generic testing: just define test_elem(R) and call test_Ring_interface(R)

include("Rings-conformance-tests.jl")

function test_elem(R::AbstractAlgebra.Floats{Float64})
   return rand(Float64)*rand(-100:100)
end
test_Ring_interface(RDF)

function test_elem(R::AbstractAlgebra.Floats{BigFloat})
   return rand(BigFloat)*rand(-100:100)
end
test_Ring_interface(RealField)

function test_elem(R::AbstractAlgebra.Integers{BigInt})
   n = big(2)^rand(1:100)
   return ZZ(rand(-n:n))
end
test_Ring_interface_recursive(ZZ)

function test_elem(R::AbstractAlgebra.Rationals{BigInt})
   n = big(2)^rand(1:100)
   return QQ(rand(-n:n)//rand(1:n))
end
test_Field_interface_recursive(QQ)

function test_elem(R::AbstractAlgebra.GFField)
   return R(rand(0:characteristic(R)))
end
test_Field_interface_recursive(GF(3))
test_Field_interface_recursive(GF(13))
test_Field_interface_recursive(GF(big(13)))
test_Field_interface_recursive(GF(big(10)^20 + 39))
