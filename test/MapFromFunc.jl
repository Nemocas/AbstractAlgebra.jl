@testset "MapFromFunc with inverse" begin
  f = map_from_func(ZZ, ZZ, x -> x+1, x -> x-1)

  @test image(f, ZZ(1)) == 2
  @test preimage(s, ZZ(1)) == 0

  @test image_fn(f)(ZZ(1)) == 2
  @test inverse_fn(f)(ZZ(1)) == 0

  @test preimage(f, image(f, ZZ(1))) == ZZ(1)

  finv = inv(f)
  y = finv(ZZ(1))
  @test preimage(finv, y) == ZZ(1)
end
