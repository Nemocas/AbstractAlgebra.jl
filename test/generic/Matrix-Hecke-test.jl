
@testset "Additional generic.mat constructors..." begin
    R, t = PolynomialRing(QQ, "t")
    S = MatrixSpace(R, 3, 3)

    @test_throws DomainError matrix(AbstractAlgebra.GFElem{Int64}[])

    @test_throws DomainError matrix(fill(one(R), 0, 5))

    @test_throws DomainError matrix(fill(one(R), 0, 5))

end
