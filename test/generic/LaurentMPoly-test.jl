function test_elem(R::AbstractAlgebra.LaurentMPolyRing{BigInt})
    n = rand(1:10)
    # R: length between 1 and 9
    # R: exponents between -n and n
    # ZZ: coeffs between -99 and 99
    rand(R, 1:9, -n:n, -99:99)
end

@testset "conformance" begin
    L, (x, y) = LaurentPolynomialRing(ZZ, ["x", "y"])
    test_Ring_interface(L)
    test_Ring_interface_recursive(L)
end
