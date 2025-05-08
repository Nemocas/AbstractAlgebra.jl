@testset "getindex" begin
  Zxy, (x, y) = ZZ[:x, :y]
  v = Zxy[x^3, 3 * x^2 * y + 1, -3, ZZ(12)]

  @test v == [x^3, 3 * x^2 * y + 1, Zxy(-3), Zxy(12)]
end

include("error-test.jl")
include("AliasMacro-test.jl")
include("Attributes-test.jl")
include("WeakKeyIdDict-test.jl")
include("WeakValueDict-test.jl")
include("Groups-test.jl")
include("NCRings-test.jl")
include("Rings-test.jl")
include("Fields-test.jl")
include("Ideals-test.jl")
include("Matrix-test.jl")
include("Modules-test.jl")
include("Maps-test.jl")
include("Benchmark-test.jl")
include("Factor-test.jl")
include("PrettyPrinting-test.jl")
include("PrintHelper-test.jl")
include("misc-test.jl")
include("Solve-test.jl")
include("MatrixNormalForms-test.jl")
include("KnownProperties.jl")
