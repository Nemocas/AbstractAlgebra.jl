@testset "known properties" begin
  R, (x, y) = ZZ[:x, :y]
  S, (u, v) = QQ[:u, :v]
  AbstractAlgebra.is_known(R::MPolyRing, ::typeof(dim)) = AbstractAlgebra.is_known(coefficient_ring(R), dim)
  AbstractAlgebra.is_known(::MPolyRing{<:Rational}, ::typeof(dim)) = true
  @test AbstractAlgebra.is_known(S, dim)
  @test_throws ErrorException AbstractAlgebra.is_known(R, dim)
  @test_throws ErrorException AbstractAlgebra.is_known(5, dim, 7, 8; dummy=true)
end
