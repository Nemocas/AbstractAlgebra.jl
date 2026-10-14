@testset "Julia.Matrix.manipulation" begin
  r, c = 1, 2
  A = zero_matrix(Int, r, c)
  A[1, :] = [10, -2]
  B = Int[10 -2]

  @test A == B
  @test number_of_rows(A) == r
  @test number_of_columns(A) == c
end

@testset "Julia.Matrix.array_creation" begin
  M = matrix_ring(ZZ, 2)
  m1 = M(ZZ[1 2; 3 4])
  m2 = M(ZZ[1 0; 0 1])
  arr = Array(M, 2)
  arr[1] = deepcopy(m1)
  arr[2] = deepcopy(m2)
  @test arr[1] == m1
  @test arr[2] == m2
end

@testset "Julia.Matrix.permutations" begin
  M = matrix_ring(ZZ, 2)
  m0 = M([1 2; 3 4])
  m1 = M([2 3; 4 5])
  m2 = M([3 4; 5 6])
  m3 = M([4 5; 6 7])
  m = elem_type(M)[m0 m1; m2 m3]
  Rm = elem_type(M)[m2 m3; m0 m1]
  Cm = elem_type(M)[m1 m0; m3 m2]
  rm = swap_rows(m, 2, 1)
  cm = swap_cols(m, 1, 2)
  @test rm == Rm
  @test cm == Cm
end
