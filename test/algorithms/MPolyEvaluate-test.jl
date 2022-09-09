@testset "MPolyEvaluate.horner" begin
  R, (x, y, z, t) = ZZ["x", "y", "z", "t"]
  @test 1 == AbstractAlgebra.evaluate_horner(one(R), [x, y, z, t])
  for i in 1:100
    f = rand(R, 1:20, 1:10, 1:10)
    g = AbstractAlgebra.evaluate_horner(f, [x, y, z, t])
    @test f == g
  end
end
