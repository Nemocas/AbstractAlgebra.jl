@testset "clear_cache" begin
  @test Generic._numpart(1000) == 24_061_467_864_032_622_473_692_149_727_991
  d = AbstractAlgebra.register_cache(Dict(1 => 2))
  v = AbstractAlgebra.register_cache([1, 2])
  cleared = Ref(false)
  clear = () -> cleared[] = true
  @test AbstractAlgebra.register_cache_clearer(clear) === clear

  @test AbstractAlgebra.clear_cache() === nothing
  @test isempty(d)
  @test isempty(v)
  @test cleared[]

  # registered while AbstractAlgebra precompiled
  @test isempty(Generic._numPartsTable)
  @test isempty(Generic._numPartsTableBig)
  @test [Generic._numpart(i) for i in 0:10] == [1,1,2,3,5,7,11,15,22,30,42]
  @test Generic._numpart(1000) == 24_061_467_864_032_622_473_692_149_727_991
end
