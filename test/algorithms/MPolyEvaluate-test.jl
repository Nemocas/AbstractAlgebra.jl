@testset "MPolyEvaluate.horner" begin
  @testset "ZZ(ZZ[x])" begin
    R, (x, y, z, t) = ZZ["x", "y", "z", "t"]
    @test iszero(AbstractAlgebra.evaluate_horner(zero(R), [x, y, z, t]))
    @test isone(AbstractAlgebra.evaluate_horner(one(R), [x, y, z, t]))
    for i in 1:100
      f = rand(R, 1:20, 1:10, 1:10)
      g = AbstractAlgebra.evaluate_horner(f, [x, y, z, t])
      @test f == g
    end
    for i in 1:10
      f = rand(R, 1:15, 1:3, 1:3)
      v = elem_type(R)[rand(R, 1:5, 1:10, 1:10) for i in 1:nvars(R)]
      a = AbstractAlgebra.evaluate_horner(f, v)
      b = evaluate(f, v)
      @test a == b
    end
  end

  @testset "ZZ(QQ)" begin
    R, (x, y, z, t) = ZZ["x", "y", "z", "t"]
    @test iszero(AbstractAlgebra.evaluate_horner(zero(R), [QQ(), QQ(), QQ(), QQ()]))
    @test isone(AbstractAlgebra.evaluate_horner(one(R), [QQ(), QQ(), QQ(), QQ()]))
    for i in 1:100
      f = rand(R, 1:20, 1:10, 1:10)
      v = elem_type(QQ)[rand(QQ, 1:20) for i in 1:nvars(R)]
      a = AbstractAlgebra.evaluate_horner(f, v)
      b = evaluate(f, v)
      @test a == b
    end
  end

  @testset "QQ(ZZ)" begin
    R, (x, y, z, t) = QQ["x", "y", "z", "t"]
    @test iszero(AbstractAlgebra.evaluate_horner(zero(R), [QQ(), QQ(), QQ(), QQ()]))
    @test isone(AbstractAlgebra.evaluate_horner(one(R), [QQ(), QQ(), QQ(), QQ()]))
    for i in 1:100
      f = rand(R, 1:20, 1:10, 1:10)
      v = elem_type(ZZ)[rand(ZZ, 1:20) for i in 1:nvars(R)]
      a = AbstractAlgebra.evaluate_horner(f, v)
      b = evaluate(f, v)
      @test a == b
    end
  end
end
