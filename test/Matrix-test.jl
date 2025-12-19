@testset "Matrix.is_diagonal" begin
    A = [1 0 0; 0 0 1; 1 0 1]
    @test !is_diagonal(A)
    @test !is_diagonal(matrix(ZZ, A))

    A = [1 0 0; 0 0 1; 0 0 1]
    @test !is_diagonal(A)
    @test !is_diagonal(matrix(ZZ, A))

    A = [1 0 0; 0 0 1]
    @test !is_diagonal(A)
    @test !is_diagonal(matrix(ZZ, A))

    A = [1 0 0; 0 1 0]
    @test is_diagonal(A)
    @test is_diagonal(matrix(ZZ, A))

    A = [1 0 0; 0 2 0; 0 0 3]
    @test is_diagonal(A)
    @test is_diagonal(matrix(ZZ, A))
end

@testset "Matrix.is_lower_triangular" begin
    A = [1 0 0; 0 0 1; 1 0 1]
    @test !is_lower_triangular(A)
    @test !is_lower_triangular(matrix(ZZ, A))

    A = [1 0 0; 0 0 1; 0 0 1]
    @test !is_lower_triangular(A)
    @test !is_lower_triangular(matrix(ZZ, A))

    A = [1 0 0; 0 0 1]
    @test !is_lower_triangular(A)
    @test !is_lower_triangular(matrix(ZZ, A))

    A = [1 0 0; 0 1 0]
    @test is_lower_triangular(A)
    @test is_lower_triangular(matrix(ZZ, A))

    A = [1 0 0; 0 2 0; 0 0 3]
    @test is_lower_triangular(A)
    @test is_lower_triangular(matrix(ZZ, A))
end

@testset "Matrix.is_upper_triangular" begin
    A = [1 0 0; 0 0 1; 1 0 1]
    @test !is_upper_triangular(A)
    @test !is_upper_triangular(matrix(ZZ, A))

    A = [1 0 0; 0 0 1; 0 0 1]
    @test is_upper_triangular(A)
    @test is_upper_triangular(matrix(ZZ, A))

    A = [1 0 0; 0 0 1]
    @test is_upper_triangular(A)
    @test is_upper_triangular(matrix(ZZ, A))

    A = [1 0 0; 0 1 0]
    @test is_upper_triangular(A)
    @test is_upper_triangular(matrix(ZZ, A))

    A = [1 0 0; 0 2 0; 0 0 3]
    @test is_upper_triangular(A)
    @test is_upper_triangular(matrix(ZZ, A))
end

@testset "Matrix.is_nilpotent" begin
  @test is_nilpotent(zero_matrix(QQ, 3, 3))
  @test is_nilpotent(upper_triangular_matrix(QQ.([0,1,1,0,1,0])))
  @test !is_nilpotent(identity_matrix(QQ, 3))
  @test !is_nilpotent(diagonal_matrix(QQ, [1,-1,0]))
end

@testset "Matrix.concat" begin
    for R in [ZZ, QQ, polynomial_ring(QQ, [:x, :y])[1]]
        S = matrix_space(R, 3, 3)
        T = matrix_space(R, 3, 6)
        U = matrix_space(R, 6, 3)

        A = S([2 3 5; 1 4 7; 9 6 3])
        B = S([1 4 7; 9 6 7; 4 3 3])
        C = matrix(R, 2, 2, [1, 2, 3, 4])

        @test hcat(A, B) == T([2 3 5 1 4 7; 1 4 7 9 6 7; 9 6 3 4 3 3])

        @test vcat(A, B) == U([2 3 5; 1 4 7; 9 6 3; 1 4 7; 9 6 7; 4 3 3])

        @test [A B] == hcat(A, B)
        @test [A A A] == hcat(A, hcat(A, A))
        @test_throws ErrorException [A C]
        @test_throws ErrorException [A A C]
        @test [A; B] == vcat(A, B)
        @test [A; B; A] == vcat(A, vcat(B, A))
        @test_throws ErrorException [A; A; C]
        @test cat(A, A; dims=(1, 2)) == block_diagonal_matrix([A, A])
        @test cat(A, A; dims=1) == hcat(A, A)
        @test cat(A, A; dims=2) == vcat(A, A)
        @test_throws ErrorException cat(A, A; dims=3)
        @test hvcat((2, 3), A, B, C, C, C) == vcat(hcat(A, B), hcat(C, C, C))
    end

    @test_throws ErrorException vcat(ZZ[1 0; 0 1], QQ[1 0; 0 1])
    @test_throws ErrorException hcat(ZZ[1 0; 0 1], QQ[1 0; 0 1])
    @test_throws ErrorException cat(ZZ[1 0; 0 1], QQ[1 0; 0 1], dims=2)
    @test_throws ErrorException hvcat(2, ZZ[1 0; 0 1], QQ[1 0; 0 1])
end

@testset "Matrix.scalar_matrix" begin
    S = scalar_matrix(3, QQ(2 // 3))
    @test S == matrix(QQ, [2//3 0 0; 0 2//3 0; 0 0 2//3])
    T = scalar_matrix(QQ, 3, 42)
    @test T == matrix(QQ, [42 0 0; 0 42 0; 0 0 42])
end

@testset "Matrix.vector_space_dim" begin
  F, = residue_field(ZZ, 2)
  S1 = matrix_space(F, 3, 4)
  S2 = matrix_ring(QQ, 2)
  @test vector_space_dim(S1) == 12
  @test vector_space_dim(S2) == 4
end

@testset "Matrix.conversion" begin
  U, t = polynomial_ring(QQ, "t")

  R = matrix_ring(U, 2)
  S = matrix_space(U, 2, 2)
  a = U.([1 2; t^2 (t-1)]) # Matrix
  Ra = R(a) # MatRingElem
  Sa = S(a) # MatElem
  @test Ra == R(Ra)
  @test Ra == R(Sa)
  @test Sa == S(Sa)
  @test Sa == S(Ra)
  @test matrix(Ra) == Sa
  @test matrix(U, Ra) == Ra
  @test matrix(Sa) == Sa
  @test matrix(U, Sa) == Sa
  @test Matrix(Ra) == a
  @test Matrix(Sa) == a
end

@testset "Matrix.keys and pairs" begin
  a = matrix(ZZ, 2, 3, [6, 3, 0, 10, 12, 14])
  @test keys(a) == CartesianIndices((2, 3))
  @test issetequal(
    keys(a),
    [CartesianIndex(1, 1), CartesianIndex(1, 2), CartesianIndex(1, 3),
     CartesianIndex(2, 1), CartesianIndex(2, 2), CartesianIndex(2, 3)],
  )
  @test issetequal(
    pairs(a),
    [CartesianIndex(1, 1) => 6, CartesianIndex(1, 2) => 3, CartesianIndex(1, 3) => 0,
     CartesianIndex(2, 1) => 10, CartesianIndex(2, 2) => 12, CartesianIndex(2, 3) => 14],
  )
end

@testset "Strassen" begin
   S = matrix(QQ, rand(-10:10, 100, 100))
   T = S*S
   TT = Strassen.mul(S, S; cutoff = 50)
   @test T == TT

   P1 = Perm(100)
   S1 = deepcopy(S)
   r1 = lu!(P1, S1)
   P = Perm(100)
   r2 = Strassen.lu!(P, S; cutoff = 50)
   @test r1 == r2
   @test S1 == S

   a = matrix(randmat_triu(matrix_ring(ZZ, 10), -100:100))
   b = matrix(randmat_triu(matrix_ring(ZZ, 10), -100:100))
   c = identity_matrix(ZZ, 10)

   Strassen.mul_tt!(c, a, b; cutoff = 5)
   @test c == a*b
   Strassen.mul_tt!(c, a, b; cutoff = 2)
   @test c == a*b
end

@testset "Promotion" begin
  M = matrix(ZZ, 1, 1, [1])
  N = matrix(QQ, 1, 1, [2])

  L = @inferred M + N
  @test base_ring(L) === QQ
  @test L == change_base_ring(QQ, M) + N
  L = @inferred M - N
  @test base_ring(L) === QQ
  @test L == change_base_ring(QQ, M) - N
  L = @inferred M * N
  @test base_ring(L) === QQ
  @test L == change_base_ring(QQ, M) * N
  L = @inferred N + M
  @test base_ring(L) === QQ
  @test L == N + change_base_ring(QQ, M)
  L = @inferred N - M
  @test base_ring(L) === QQ
  @test L == N - change_base_ring(QQ, M)
  L = @inferred N * M
  @test base_ring(L) === QQ
  @test L == N * change_base_ring(QQ, M)

  @test M * QQ[1;] == QQ[1;]
  @test M * ZZ[1;] == ZZ[1;]
  @test N * QQ[1;] == QQ[2;]
  @test N * ZZ[1;] == QQ[2;]
  @test QQ[1;] * M == QQ[1;]
  @test ZZ[1;] * M == ZZ[1;]
  @test QQ[1;] * N == QQ[2;]
  @test ZZ[1;] * N == QQ[2;]

  @test M * QQ(1) == QQ(1) * M == QQ.(M)
  @test N * ZZ(1) == ZZ(1) * N == QQ.(N)
end

@testset "Rank Interpolation" begin
  Fx, x = polynomial_ring(GF(29), :x)
  Qy, y = polynomial_ring(QQ, 3, :y)
  Qz, z = rational_function_field(QQ, :z)
  
  A = rand(matrix_space(Fx, 2, 2), 0:5)
  B = rand(matrix_space(Qy, 2, 2), 1:3, 1:4, 1:2)
  C = rand(matrix_space(Qz, 2, 2), 0:5, 1:5)

  @test rank(A) == AbstractAlgebra.rank_interpolation(A)
  @test rank(A) >= AbstractAlgebra.rank_interpolation_mc(A)
  @test rank(B) == AbstractAlgebra.rank_interpolation(B)
  @test rank(B) >= AbstractAlgebra.rank_interpolation_mc(B)
  @test rank(C) == AbstractAlgebra.rank_interpolation(C)
  @test rank(C) >= AbstractAlgebra.rank_interpolation_mc(C)

  multiply_row!(A, 0, 1)
  multiply_row!(B, 0, 1)
  multiply_row!(C, 0, 1)

  @test rank(A) == AbstractAlgebra.rank_interpolation(A)
  @test rank(A) >= AbstractAlgebra.rank_interpolation_mc(A)
  @test rank(B) == AbstractAlgebra.rank_interpolation(B)
  @test rank(B) >= AbstractAlgebra.rank_interpolation_mc(B)
  @test rank(C) == AbstractAlgebra.rank_interpolation(C)
  @test rank(C) >= AbstractAlgebra.rank_interpolation_mc(C)

  A = rand(matrix_space(Fx, 0, 2), 0:5, 1:5)
  B = rand(matrix_space(Qy, 0, 2), 1:3, 1:4, 1:2)
  C = rand(matrix_space(Qz, 2, 0), 0:5, 1:5)

  @test rank(A) == AbstractAlgebra.rank_interpolation(A)
  @test rank(A) >= AbstractAlgebra.rank_interpolation_mc(A)
  @test rank(B) == AbstractAlgebra.rank_interpolation(B)
  @test rank(B) >= AbstractAlgebra.rank_interpolation_mc(B)
  @test rank(C) == AbstractAlgebra.rank_interpolation(C)
  @test rank(C) >= AbstractAlgebra.rank_interpolation_mc(C)

  A = rand(matrix_space(Fx, 10, 10), 0:5, 1:5)
  # @test_throws AbstractAlgebra.rank_interpolation(A)
  # @test_throws AbstractAlgebra.rank_interpolation_mc(A)
end
