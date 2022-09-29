function test_evaluate(f, v, q = evaluate(f, v))
  e1 = AbstractAlgebra.evaluate_horner(f, v)
  @test typeof(q) == typeof(e1)
  @test q == e1

  e2 = AbstractAlgebra.evaluate_log(f, v)
  @test typeof(q) == typeof(e2)
  @test q == e2

  e3 = AbstractAlgebra.evaluate_log([f, f], v)
  @test typeof(q) == typeof(e3[1]) == typeof(e3[2])
  @test q == e3[1]
  @test q == e3[2]

  e4 = AbstractAlgebra.evaluate_horner_lex(f, v)
  @test typeof(q) == typeof(e4)
  @test q == e4
end

@testset "MPolyEvaluate.ZZ(ZZ[x])" begin
  R, (x, y, z, t) = ZZ["x", "y", "z", "t"]
  test_evaluate(zero(R), [x, y, z, t])
  test_evaluate(one(R), [x, y, z, t])
  for i in 1:100
    f = rand(R, 1:20, 1:10, 1:10)
    test_evaluate(f, [x, y, z, t], f)
  end
  for i in 1:10
    f = rand(R, 1:15, 1:3, 1:3)
    v = elem_type(R)[rand(R, 1:5, 1:10, 1:10) for i in 1:nvars(R)]
    test_evaluate(f, v)
  end
end

@testset "MPolyEvaluate.ZZ(QQ)" begin
  R, (x, y, z, t) = PolynomialRing(ZZ, ["x", "y", "z", "t"], ordering = :degrevlex)
  test_evaluate(zero(R), [QQ(), QQ(), QQ(), QQ()])
  test_evaluate(one(R), [QQ(), QQ(), QQ(), QQ()])
  for i in 1:100
    f = rand(R, 1:20, 1:10, 1:10)
    v = elem_type(ZZ)[rand(ZZ, 1:5) for i in 1:nvars(R)]
    test_evaluate(f, v)
  end
end

@testset "MPolyEvaluate.QQ(ZZ)" begin
  R, (x, y, z, t) = QQ["x", "y", "z", "t"]
  test_evaluate(zero(R), [ZZ(), ZZ(), ZZ(), ZZ()])
  test_evaluate(one(R), [ZZ(), ZZ(), ZZ(), ZZ()])
  for i in 1:100
    f = rand(R, 1:20, 1:10, 1:10)
    v = elem_type(ZZ)[rand(ZZ, 1:20) for i in 1:nvars(R)]
    test_evaluate(f, v)
  end
end

