@testset "known properties" begin
  R, (x, y) = ZZ[:x, :y]
  S, (u, v) = QQ[:u, :v]
  AbstractAlgebra.is_known(::typeof(dim), R::MPolyRing) = AbstractAlgebra.is_known(dim, coefficient_ring(R))
  AbstractAlgebra.is_known(::typeof(dim), ::MPolyRing{<:Rational}) = true
  @test AbstractAlgebra.is_known(dim, S)
  @test_throws MethodError AbstractAlgebra.is_known(dim, R)
  @test_throws MethodError AbstractAlgebra.is_known(dim, 5, 7, 8; dummy=true)
end
