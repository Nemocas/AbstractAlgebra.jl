@testset "broadcasting" begin
  A = ZZ[1 2; 3 4]
  B = QQ[0 0; 0 0]
  f = x -> QQ(x)^3
  B .= (f.(3 .+ A .* A .* ZZ(3)))
  @test B == map_entries(x -> f(3 + x * x * 3), A)
  k = GF(5)
  kx, x = k["x"]
  A = kx[x x^2; 1 x]
  B = ((x, y) -> x(3) * y(1)).(k(2) .* A, A)
  @test B == map_entries(x -> k(2) * x(3) * x(1), A)

  @test_throws ErrorException [1, 2] .* A
  @test_throws ErrorException A .* 2 .* 2 .* [1, 2]
end

@testset "broadcasting evaluate" begin
  Qa, (k1, k2, k3, k4) = rational_function_field(QQ, ["k1", "k2", "k3", "k4"])
  A = matrix(Qa, [k1 k2; k3 k4])
  @test evaluate.(A, Ref([1, 2, 3, 4])) == QQ[1 2; 3 4]
end

