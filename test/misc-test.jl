@testset "evaluate" begin
  Qx, (x, y) = QQ["x", "y"]
  @test evaluate(:((x^2 + 1)//(y + 1//2)), Qx) == (x^2 + 1)//(y + 1//2)

  Qx, x = QQ["x"]
  @test evaluate(:(x - 1), Qx) == x - 1
  @test evaluate(:(-x), Qx) == -x
  @test evaluate(:(x*x*x*x + 1), Qx) == x^4 + 1

  @test evaluate(:((x^2 + 1)//(y + 1//2)), Dict(:x => 1, :y => QQ(2))) == 4//5

  @test_throws ErrorException evaluate(:(sin(x)), Dict(:x => 1))
  @test_throws ErrorException evaluate(:(x + x^2 + f(x)), Dict(:x => 1))
  @test_throws ErrorException evaluate((Expr(:x)), Dict(:y => 1))
end
