include("julia/Integers-test.jl")
include("julia/Matrix-test.jl")
include("broadcasting-test.jl")

# artificially low cutoffs for testing purposes
function AbstractAlgebra.hgcd_prefers_basecase(a::T, b::T) where T <: PolyRingElem
   return degree(b) < 5
end

function AbstractAlgebra.mat22_mul_prefers_classical(
   a11::T, a12::T, a21::T, a22::T,
   b11::T, b12::T, b21::T, b22::T
) where T <: PolyRingElem

   return degree(a11) + degree(a22) < 4 || degree(b11) + degree(b22) < 4
end

include("generic/Poly-test.jl")
include("generic/SparsePoly-test.jl")
include("generic/LaurentPoly-test.jl")
include("generic/Residue-test.jl")
include("generic/ResidueField-test.jl")
include("generic/TotalFraction-test.jl")
include("generic/RelSeries-test.jl")
include("generic/AbsSeries-test.jl")
include("generic/LaurentSeries-test.jl")
include("generic/PuiseuxSeries-test.jl")
include("generic/AbsMSeries-test.jl")
include("generic/Matrix-test.jl")
include("generic/MPoly-test.jl")
include("generic/LaurentMPoly-test.jl")
include("generic/UnivPoly-test.jl")
include("algorithms/MPolyEvaluate-test.jl")
include("algorithms/MPolyFactor-test.jl")
include("algorithms/MPolyNested-test.jl")
include("algorithms/DensePoly-test.jl")
include("algorithms/coprime_base-test.jl")
include("generic/PolyRingHom-test.jl")

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

@testset "ppio" begin
   Qx, x = QQ["x"]
   p, q = @inferred ppio(zero(Qx), x)
   @test isone(p) && q == zero(Qx)

   @testset "ppio($a,$b) for $T" for a in -25:25, b in -25:25, T in [Int, BigInt, Float64]
     c, d = ppio(T(a),T(b))
     @test a == c * d
     @test a == 0 || c == gcd(T(a), T(b)^5)
   end
end

@testset "sqrt" begin
  Base.sqrt(a::EuclideanRingResidueRingElem{BigInt}) = parent(a)(BigInt(sqrt(a.data)))
  R = residue_ring(ZZ, 101)[1]

  @test sqrt(R(81), check=true) == R(9)
  @test_throws ErrorException sqrt(R(82), check=true)
  @test is_square_with_sqrt(R(16)) == (true, R(4))
  @test is_square_with_sqrt(R(15))[1] == false
end

@testset "properties" begin
  @test is_perfect(QQ)
  @test is_perfect(GF(2))
  @test is_finite(GF(2))
  @test !is_finite(QQ)
end

@testset "is_associated" begin
   F = GF(13)
   for x in F, y in F
      @test is_associated(x, y) == (iszero(x) == iszero(y))
   end

   Fx, x = polynomial_ring(F)
   f = x^3 + 4x + 3
   @test is_associated(f, f)
   @test is_associated(f, 6*f)
   @test is_associated(7*f, f)
   @test !is_associated(f, 0*f)
   @test !is_associated(0*f, f)
   @test !is_associated(f, 2x^2-x+2)
end
