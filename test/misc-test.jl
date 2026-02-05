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

@testset "rising_factorial" begin
  @test_throws OverflowError rising_factorial(10,30)
  @test_throws OverflowError rising_factorial2(10,30)
end

### Test options to @varnames_interface ###
# For tests under default options, see doctests and test/generic/MPoly-test.jl

# note: macro definition not allowed inside a local scope
# hence using a module instead:
module VarNamesTest
  using AbstractAlgebra
  using Test

  f(a, s::Vector{Symbol}) = a, string.((a,), s)
  projective(a, s::Vector{Symbol}) = f(a, s)
  uses_n(n, s::Vector{Symbol}) = f(n, s)
  no_n_variant(a, s::Vector{Symbol}) = f(a, s)
  no_macros(a, s::Vector{Symbol}) = f(a, s)
  AbstractAlgebra.@varnames_interface f(a, s)
  AbstractAlgebra.@varnames_interface projective(a, s) range=0:n
  AbstractAlgebra.@varnames_interface uses_n(n, s) n=m range=1:m
  AbstractAlgebra.@varnames_interface no_n_variant(a, s) n=:no
  AbstractAlgebra.@varnames_interface no_macros(a, s) macros=:no

  @testset "VarNames" begin
    @test f("A", 3) == ("A", ["Ax1", "Ax2", "Ax3"])

    @f("A", :y => 1:3)
    @f("B", 3, "z")

    @test [y1, y2, y3] == ["Ay1", "Ay2", "Ay3"]
    @test [z1, z2, z3] == ["Bz1", "Bz2", "Bz3"]
    @test ! @isdefined x1
    @test ! @isdefined y4
    @test ! @isdefined z4

    @testset "VarNames.options" begin
      @test uses_n("A", 3) == ("A", ["Ax1", "Ax2", "Ax3"])
      @test projective("A", 3) == ("A", ["Ax0", "Ax1", "Ax2", "Ax3"])

      @projective("A", 3)
      @test [x0, x1, x2, x3] == ["Ax0", "Ax1", "Ax2", "Ax3"]

      @uses_n("B", 3, 'y')
      @test [y1, y2, y3] == ["By1", "By2", "By3"]

      @test ! @isdefined var"@no_macros"
      @test_throws MethodError no_n_variant("A", 3)
      @test_throws ArgumentError @macroexpand @no_n_variant("A", 3, :q)
      @test ! @isdefined q1
    end
  end
end

@testset "is_pairwise" begin
  @test is_pairwise(!is_associated, [2, -3, 5])
  @test !is_pairwise(!is_associated, [2, -3, -2])
end
