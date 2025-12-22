@testset "Localization" begin

  Qx, x = QQ["x"]

   @testset "Conformance" begin
      L = localization(Qx, x^2 + 1)
      ConformanceTests.test_Ring_interface(L)
   end

  @testset "Constructor" begin

    @test parent_type(LocalizedEuclideanRingElem{elem_type(Qx)}) == LocalizedEuclideanRing{elem_type(Qx)}

    L = localization(Qx, x^2 + 1)
    @test elem_type(L) == LocalizedEuclideanRingElem{elem_type(Qx)}
    @test base_ring(L) == Qx
    @test base_ring(L()) == Qx

    L = localization(Qx, 4 * x^2 + x + 1)
    @test elem_type(L) == LocalizedEuclideanRingElem{elem_type(Qx)}
    @test base_ring(L) == Qx
    @test base_ring(L()) == Qx

    L = localization(Qx, [3 * x^4 + x + 18, x^9 + 3])
    @test elem_type(L) == LocalizedEuclideanRingElem{elem_type(Qx)}
    @test base_ring(L) == Qx
    @test base_ring(L()) == Qx
  end

  @testset "Canonicalisation" begin

    L = localization(Qx, 5x^2 + 3x + 2)
    @test canonical_unit(L(5x^2 + 3x + 2)) == L(5)
    @test canonical_unit(L(28)) == L(28)
    @test canonical_unit(L(5 * x^4 + 3 * x^3 + 12 * x^2 + 6 * x + 4)) == L(5 * (x^2 + 2))

    L = localization(Qx, [5x^2 + 3x + 2, x^16 + x + 3])
    @test canonical_unit(L(x^16 + x + 3)) == one(L)
    @test canonical_unit(L(89)) == L(89)
    @test canonical_unit(L((5 * x^4 + 3 * x^3 + 12 * x^2 + 6 * x + 4) // (x^3 + 9))) == inv(L((1 // 5 * x^3 + 9 // 5) // (x^2 + 2)))
  end

  @testset "Parent object call overloading" begin

    L1 = localization(Qx, 5x^2 + 3x + 2)
    L2 = localization(Qx, 7x^3 + 2)
    @test_throws ErrorException L1(L2(5))
    @test parent(L1(4x // 7)) == L1

    L = localization(Qx, [x^2 + 2, 3x^3 + 5])
    @test parent(L(3x^2 + 6)) == L
    @test_throws ErrorException L(15 // (3x^2 + 6))
    @test_throws ErrorException L(17 // (9 * x^5 + 18 * x^3 + 15 * x^2 + 30))
    @test_throws ErrorException L(19 // (9 * x^6 + 54 * x^5 + 18 * x^4 + 123 * x^3 + 90 * x^2 + 30 * x + 180))
  end

  @testset "GCDX" begin

    L = localization(Qx, x^6 + 108)
    a = L(x^10 + 2 * x^8 + 14 * x^7 + 49 * x^5 + 42 * x^3 + 294 * x^2 + 588)
    b = L(x^11 + 14 * x^8 + 21 * x^6 + 294 * x^3 + 6)
    (g, u, v) = gcdx(a, b)
    @test g == u * a + v * b

    L = localization(Qx, [x^6 + 108, x^8 + 12x^3 + 3])
    a = L(x^10 + 2 * x^8 + 14 * x^7 + 49 * x^5 + 42 * x^3 + 294 * x^2 + 588)
    b = L(x^11 + 14 * x^8 + 21 * x^6 + 294 * x^3 + 6)
    (g, u, v) = gcdx(a, b)
    @test g == u * a + v * b
  end

  @testset "GCD and LCM" begin

    L = localization(Qx, x^6 + 108)
    a = L(x^10 + 2 * x^8 + 14 * x^7 + 49 * x^5 + 42 * x^3 + 294 * x^2 + 588)
    b = L(x^11 + 14 * x^8 + 21 * x^6 + 294 * x^3 + 6)
    g = gcd(a, b)
    @test g == L(1)
    a = L(x^16 + 2 * x^14 + 14 * x^13 + 49 * x^11 + 108 * x^10 + 42 * x^9 + 510 * x^8 + 1512 * x^7 + 588 * x^6 + 5292 * x^5 + 4536 * x^3 + 31752 * x^2 + 63504)
    b = L(x^17 + 14 * x^14 + 21 * x^12 + 108 * x^11 + 294 * x^9 + 1512 * x^8 + 2274 * x^6 + 31752 * x^3 + 648)
    g = gcd(a, b)
    @test g == L(x^6 + 108)

    L = localization(Qx, [x^6 + 108, x^8 + 12x^3 + 3])
    a = L(x^10 + 2 * x^8 + 14 * x^7 + 49 * x^5 + 42 * x^3 + 294 * x^2 + 588)
    b = L(x^13 + 2 * x^11 + 14 * x^10 + 49 * x^8 + 42 * x^6 + 294 * x^5 + 588 * x^3 + 6 * x^2 + 12)
    g = gcd(a, b)
    @test g == L(1)
  end

  @testset "Exact division" begin
    L = localization(Qx, x^6 + 108)
    @test_throws ErrorException divexact(L(9), L(2x^6 + 216))
    a = L(x^10 + 2 * x^8 + 14 * x^7 + 49 * x^5 + 42 * x^3 + 294 * x^2 + 588)
    b = L(x^11 + 14 * x^8 + 21 * x^6 + 294 * x^3 + 6)
    g = divides(a, b)
    @test a == b * g[2]
    @test g[1] == true
    a = L(x^16 + 2 * x^14 + 14 * x^13 + 49 * x^11 + 108 * x^10 + 42 * x^9 + 510 * x^8 + 1512 * x^7 + 588 * x^6 + 5292 * x^5 + 4536 * x^3 + 31752 * x^2 + 63504)
    b = L(x^17 + 14 * x^14 + 21 * x^12 + 108 * x^11 + 294 * x^9 + 1512 * x^8 + 2274 * x^6 + 31752 * x^3 + 648)
    @test a == b * g[2]
    @test g[1] == true
  end

  @testset "Inversion" begin
    L = localization(Qx, x^6 + 108)
    @test_throws ErrorException inv(L(2x^6 + 216))
    @test inv(L(x^2 + 108)) == L(1 // (x^2 + 108))

    L = localization(Qx, [x^6 + 3, x^3 + 5])
    @test_throws ErrorException inv(L(2x^6 + 6))
    @test inv(L(x^2 + 108)) == L(1 // (x^2 + 108))
  end

  @testset "Binary operators" begin
    L = localization(Qx, x^6 + 108)
    @test L(18x // 2) + L(2x // 1) == L(11x)
    @test L(18x // 3) - L(1x) == L(5x)
    @test L(32x) * L(4x) == L(128x^2)

    L = localization(Qx, [x^6 + 3, x^3 + 5])
    @test L(18x // 2) + L(2x // 1) == L(11x)
    @test L(18x // 3) - L(1x) == L(5x)
    @test L(32x) * L(4x) == L(128x^2)
  end

  @testset "Basic manipulation" begin
    L = localization(Qx, x^6 + 108)
    @test iszero(L(x^2 + x)) == false
    @test isone(L(x^2 // x^2)) == true
    @test is_unit(L(x)) == true

    L = localization(Qx, [x^6 + 3, x^3 + 5])
    @test iszero(L(x - x)) == true
    @test isone(L(x^2 + 3x)) == false
    @test is_unit(L((x^3 + 5) // (x^3 + 5)))
  end
end
