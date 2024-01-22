@testset "Echelon form" begin
  M = matrix(QQ, [ 1 2 3 4 ; 5 6 7 8 ; 9 10 11 12 ])
  R = echelon_form(M)
  @test is_rref(R)
  @test nrows(R) == 3
  S = echelon_form(M, shape = :lower)
  R = reverse_rows(reverse_cols(S))
  @test is_rref(R)
  @test nrows(S) == 3
  R = echelon_form(M, trim = true)
  @test is_rref(R)
  @test nrows(R) == 2
  S = echelon_form(M, shape = :lower, trim = true)
  R = reverse_rows(reverse_cols(S))
  @test is_rref(R)
  @test nrows(S) == 2

  R, U = echelon_form_with_transformation(M)
  @test is_rref(R)
  @test U*M == R
  @test is_invertible(U)
  S, U = echelon_form_with_transformation(M, shape = :lower)
  R = reverse_rows(reverse_cols(S))
  @test is_rref(R)
  @test U*M == S
end

@testset "Hermite form" begin
  M = matrix(ZZ, [ 1 2 3 4 ; 5 6 7 8 ; 9 10 11 12 ])
  H = hermite_form(M)
  @test is_hnf(H)
  @test nrows(H) == 3
  Hl = hermite_form(M, shape = :lower)
  H = reverse_rows(reverse_cols(Hl))
  @test is_hnf(H)
  @test nrows(Hl) == 3
  H = hermite_form(M, trim = true)
  @test is_hnf(H)
  @test nrows(H) == 2
  Hl = hermite_form(M, shape = :lower, trim = true)
  H = reverse_rows(reverse_cols(Hl))
  @test is_hnf(H)
  @test nrows(Hl) == 2

  H, U = hermite_form_with_transformation(M)
  @test is_hnf(H)
  @test U*M == H
  @test is_invertible(U)
  Hl, U = hermite_form_with_transformation(M, shape = :lower)
  H = reverse_rows(reverse_cols(Hl))
  @test is_hnf(H)
  @test U*M == Hl
end
