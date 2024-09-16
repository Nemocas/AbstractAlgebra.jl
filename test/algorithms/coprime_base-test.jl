@testset "coprime base" begin
  c = BigInt[6, 10]
  d = coprime_base(c)
  @test issetequal(d, BigInt[2, 3, 5])

  c = BigInt[1]
  d = coprime_base(c)
  @test issetequal(d, BigInt[1])

  c = BigInt[6, 10, 1]
  d = coprime_base(c)
  @test issetequal(d, BigInt[2, 3, 5])
  coprime_base_push!(d, BigInt(1))
  @test issetequal(d, BigInt[2, 3, 5])
  coprime_base_push!(d, BigInt(14))
  @test issetequal(d, BigInt[2, 3, 5, 7])
end
