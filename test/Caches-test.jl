module CachesTest
  using AbstractAlgebra
  AbstractAlgebra.@register_cache const d = Dict(1 => 2)
  AbstractAlgebra.@register_cache const v = [1, 2]
  AbstractAlgebra.@register_cache const s = Set([1, 2, 3])
end

@testset "clear_cache" begin
  @test_throws ArgumentError @macroexpand AbstractAlgebra.@register_cache d = Dict()
  @test_throws ArgumentError @macroexpand AbstractAlgebra.@register_cache const d::Dict = Dict()

  @test Generic._numpart(1000) == 24_061_467_864_032_622_473_692_149_727_991
  io = IOBuffer()
  AbstractAlgebra.print_caches(io)
  s = String(take!(io))
  @test occursin("Main.CachesTest.d: length 1, $(@__FILE__):3\n", s)
  @test occursin("Main.CachesTest.v: length 2, $(@__FILE__):4\n", s)
  @test occursin("Main.CachesTest.s: length 3, $(@__FILE__):5\n", s)
  @test occursin("AbstractAlgebra.Generic._numPartsTableBig: length ", s)

  @test AbstractAlgebra.clear_cache() === nothing
  @test isempty(CachesTest.d)
  @test isempty(CachesTest.v)
  @test isempty(CachesTest.s)

  # registered while AbstractAlgebra precompiled
  @test isempty(Generic._numPartsTable)
  @test isempty(Generic._numPartsTableBig)
  @test [Generic._numpart(i) for i in 0:10] == [1,1,2,3,5,7,11,15,22,30,42]
  @test Generic._numpart(1000) == 24_061_467_864_032_622_473_692_149_727_991
end
