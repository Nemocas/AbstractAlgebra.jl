include("julia/Integers-test.jl")
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
include("UniversalRing-test.jl")
include("generic/Localization-test.jl")
include("algorithms/MPolyEvaluate-test.jl")
include("algorithms/MPolyFactor-test.jl")
include("algorithms/MPolyNested-test.jl")
include("algorithms/DensePoly-test.jl")
include("algorithms/coprime_base-test.jl")
include("generic/PolyRingHom-test.jl")
include("PuiseuxPolynomials-test.jl")

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

@testset "properties" begin
  @test is_perfect(QQ)
  @test is_perfect(GF(2))
  @test is_finite(GF(2))
  @test !is_finite(QQ)
  R1, = laurent_polynomial_ring(QQ[:x][1], :y)
  @test krull_dim(R1) == 2
  R2, = power_series_ring(ZZ, 3, :x)
  @test krull_dim(R2) == 2
  S, = universal_polynomial_ring(ZZ, [:x,:y])
  @test_throws NotImplementedError is_noetherian(S)
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

# `promote` returns copies although the operands already have a common type
mutable struct CopyPromoteRingElem <: RingElem end

Base.promote(x::CopyPromoteRingElem, y::CopyPromoteRingElem) =
   CopyPromoteRingElem(), CopyPromoteRingElem()

@testset "divexact fallback" begin
   @test_throws NotImplementedError divexact(CopyPromoteRingElem(), CopyPromoteRingElem())
end

@testset "equality across parents" begin
   # pairs of elements with the same representation in two incompatible
   # parents, so that their hashes may collide
   twice(f) = (f(), f())
   Zx, = polynomial_ring(ZZ, :x)
   U = universal_polynomial_ring(QQ; cached = false)
   V = universal_polynomial_ring(QQ; cached = false)

   pairs = [
      "Poly" => twice(() -> gen(polynomial_ring(QQ, :x; cached = false)[1])),
      "NCPoly" => twice(() -> gen(polynomial_ring(matrix_ring(QQ, 2), :x; cached = false)[1])),
      "MPoly" => twice(() -> gen(polynomial_ring(QQ, [:x, :y]; cached = false)[1], 1)),
      "LaurentPoly" => twice(() -> gen(laurent_polynomial_ring(QQ, :x; cached = false)[1])),
      "LaurentMPoly" => twice(() -> gen(laurent_polynomial_ring(QQ, [:x, :y]; cached = false)[1], 1)),
      "SparsePoly" => twice(() -> Generic.SparsePolynomialRing(QQ, :x; cached = false)[2]),
      "AbsSeries" => twice(() -> gen(power_series_ring(QQ, 10, :t; model = :capped_absolute, cached = false)[1])),
      "RelSeries" => twice(() -> gen(power_series_ring(QQ, 10, :t; cached = false)[1])),
      "LaurentSeriesRing" => twice(() -> gen(laurent_series_ring(ZZ, 10, :t; cached = false)[1])),
      "LaurentSeriesField" => twice(() -> gen(laurent_series_field(QQ, 10, :t; cached = false)[1])),
      "PuiseuxSeries" => twice(() -> gen(puiseux_series_ring(ZZ, 10, :t; cached = false)[1])),
      "AbsMSeries" => twice(() -> gen(power_series_ring(QQ, [5, 5], [:x, :y]; cached = false)[1], 1)),
      "Fraction" => twice(() -> fraction_field(Zx; cached = false)(gen(Zx))),
      "FactoredFraction" => twice(() -> Generic.FactoredFractionField(ZZ; cached = false)(3)),
      "TotFrac" => twice(() -> total_ring_of_fractions(Zx; cached = false)(gen(Zx))),
      "Residue" => (residue_ring(ZZ, 7)[1](3), residue_ring(ZZ, 11)[1](3)),
      "ResidueField" => (residue_field(ZZ, 7)[1](3), residue_field(ZZ, 11)[1](3)),
      "RationalFunctionField" => twice(() -> gen(rational_function_field(QQ, :x; cached = false)[1])),
      "FreeAssociativeAlgebra" => twice(() -> gen(free_associative_algebra(QQ, [:x, :y]; cached = false)[1], 1)),
      "MatRing" => (one(matrix_ring(ZZ, 2)), one(matrix_ring(ZZ, 3))),
      "GF" => (GF(7)(3), GF(11)(3)),
      "UnivPoly" => (gen(U, :x), gen(V, :x)),
   ]

   @testset "$name" for (name, (a, b)) in pairs
      @test a == a
      @test isequal(a, a)
      @test !isequal(a, b)
      @test !isequal(b, a)
      @test length(Dict(a => 1, b => 2)) == 2
   end
end
