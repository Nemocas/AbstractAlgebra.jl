@testset "Matrix.is_lower_triangular" begin
    A = [ 1 0 0 ; 0 0 1 ; 1 0 1 ]
    @test !is_lower_triangular(A)
    @test !is_lower_triangular(matrix(ZZ, A))

    A = [ 1 0 0 ; 0 0 1 ; 0 0 1 ]
    @test !is_lower_triangular(A)
    @test !is_lower_triangular(matrix(ZZ, A))

    A = [ 1 0 0 ; 0 0 1 ]
    @test !is_lower_triangular(A)
    @test !is_lower_triangular(matrix(ZZ, A))

    A = [ 1 0 0 ; 0 1 0 ]
    @test is_lower_triangular(A)
    @test is_lower_triangular(matrix(ZZ, A))
end

@testset "Matrix.is_upper_triangular" begin
    A = [ 1 0 0 ; 0 0 1 ; 1 0 1 ]
    @test !is_upper_triangular(A)
    @test !is_upper_triangular(matrix(ZZ, A))

    A = [ 1 0 0 ; 0 0 1 ; 0 0 1 ]
    @test is_upper_triangular(A)
    @test is_upper_triangular(matrix(ZZ, A))

    A = [ 1 0 0 ; 0 0 1 ]
    @test is_upper_triangular(A)
    @test is_upper_triangular(matrix(ZZ, A))

    A = [ 1 0 0 ; 0 1 0 ]
    @test is_upper_triangular(A)
    @test is_upper_triangular(matrix(ZZ, A))
end
