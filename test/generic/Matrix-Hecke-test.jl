
@testset "Additional generic.mat constructors..." begin
    R, t = PolynomialRing(QQ, "t")
    S = MatrixSpace(R, 3, 3)

    @test_throws DomainError matrix(AbstractAlgebra.GFElem{Int64}[])

    @test_throws DomainError matrix(fill(one(R), 0, 5))

    @test_throws DomainError matrix(fill(one(R), 0, 5))

    @test_throws DomainError diagonal_matrix(fill(one(R),0))

    @test isdiagonal(diagonal_matrix(fill(one(R), 5)))

    @test isdiagonal(diagonal_matrix(fill(zero(R), 5)))

    @test isdiagonal(matrix(fill(one(R), 5, 5))) == false

    @test isdiagonal(matrix(R,0,0,fill(one(R),0,0)))

    
end
