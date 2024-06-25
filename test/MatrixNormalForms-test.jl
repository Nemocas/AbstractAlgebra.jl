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
  @test is_invertible(U)
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
  @test is_invertible(U)
end

@testset "Howell form" begin
  R, _ = residue_ring(ZZ, 30)
  M = R[2 3 4 5; 6 9 12 15; 10 15 20 30]
  H = howell_form(M)
  @test H == R[2 3 4 0; 0 15 0 0; 0 0 0 5; 0 0 0 0]
  Hl = howell_form(M, shape = :lower)
  @test Hl == R[0 0 0 0; 0 15 0 0; 16 9 2 0; 0 0 0 5]
  H = howell_form(M, trim = true)
  @test H == R[2 3 4 0; 0 15 0 0; 0 0 0 5]
  Hl = howell_form(M, shape = :lower, trim = true)
  @test Hl == R[0 15 0 0; 16 9 2 0; 0 0 0 5]

  H, U = howell_form_with_transformation(M)
  @test H == R[2 3 4 0; 0 15 0 0; 0 0 0 5; 0 0 0 0]
  @test U*M == H
  Hl, U = howell_form_with_transformation(M, shape = :lower)
  @test Hl == R[0 0 0 0; 0 15 0 0; 16 9 2 0; 0 0 0 5]
  @test U*M == Hl

  R, _ = residue_ring(ZZ, 12)
  M = R[4 1 0; 3 0 0; 0 0 5]
  @test howell_form(M) == R[1 1 0; 0 3 0; 0 0 1]

  a = R[4 1 0; 0 0 5; 0 0 0]
  b = R[8 5 5; 0 9 8; 0 0 10]
  c = R[4 1 0; 0 3 0; 0 0 1]
  @test howell_form(a) == c
  @test howell_form(b) == c

  Qt, t = QQ["t"]
  R, _ = residue_ring(Qt, t^5)
  M = R[t^2 t^2 - 1 1; t t - 1 0]
  H = howell_form(M)
  @test H == R[t 0 -1; 0 1 -t^3 - t^2 - t - 1; 0 0 t^4]
  H = howell_form(M, shape = :lower)
  @test H == R[0 0 0; -t^4 - t^3 - t^2 - t 1 0; -t 0 1]
  H = howell_form(M, trim = true)
  @test H == R[t 0 -1; 0 1 -t^3 - t^2 - t - 1; 0 0 t^4]
  H = howell_form(M, shape = :lower, trim = true)
  @test H == R[-t^4 - t^3 - t^2 - t 1 0; -t 0 1]
  H, U = howell_form_with_transformation(M)
  @test H == R[t 0 -1; 0 1 -t^3 - t^2 - t - 1; 0 0 t^4]
  @test U*M == H
  H, U = howell_form_with_transformation(M, shape = :lower)
  @test H == R[0 0 0; -t^4 - t^3 - t^2 - t 1 0; -t 0 1]
  @test U*M == H
end
