
@testset "Addtional Matrix functionality from Hecke..." begin
    
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

    @testset "Additional generic.mat AbstractArray behaviour..." begin

        A = diagonal_matrix(QQ, [1,2,3])

        @test (stride(A,-1), stride(A,1), stride(A,2), stride(A,3)) == (1,1,3,9)

        @test strides(A) == (1,3)

        @test size([x for x in A]) == size(A)

        # TODO: Add some getindex setindex tests for the 1-index accessor.
        
    end

    @testset "Additional generic.mat char/minpoly behaviour..." begin

        # Check to see if the functions exist. Beyond that, behaviour is covered
        # by other tests.

        A = matrix(QQ, fill(1, 5, 5))

        @test isa(charpoly(A), Any) && isa(minpoly(A), Any)

    end
end

nothing
