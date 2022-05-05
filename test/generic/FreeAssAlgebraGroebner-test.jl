@testset "Generic.FreeAssAlgebra.groebner" begin

   R, (x, y, u, v, t, s) = FreeAssociativeAlgebra(GF(2), ["x", "y", "u", "v", "t", "s"])
   g = AbstractAlgebra.groebner_basis([u*(x*y)^3 + u*(x*y)^2 + u + v,
                                      (y*x)^3*t + (y*x)^2*t + t + s])
   @test length(g) >= 5

   # Example 6.1 Kreuzer & Xiu
   R, (a, b) = FreeAssociativeAlgebra(QQ, ["a", "b"])

   g = AbstractAlgebra.groebner_basis([a^2 - 1, b^3 - 1, (a*b*a*b^2)^2 - 1])
   AbstractAlgebra.interreduce!(g)
   @test length(g) == 5

   g = AbstractAlgebra.groebner_basis([a^2 - 1, b^3 - 1, (a*b*a*b*a*b^2)^2 - 1])
   AbstractAlgebra.interreduce!(g)
   @test length(g) == 15
end
