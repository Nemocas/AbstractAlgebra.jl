# We test everything over the following rings
# * GF(101) -> uses RREFTrait or LUTrait
# * QQ -> uses FFLUTrait
# * ZZ -> uses HermiteFormTrait
# * residue_ring(ZZ, 16) -> uses HowellFormTrait
# * fraction_field(QQ["x"][1]) -> uses MatrixInterpolateTrait

@testset "Linear solving over $R" for R in [ QQ, ZZ, GF(101), fraction_field(QQ["x"][1]), residue_ring(ZZ, 16)[1] ]
  M = matrix(R, [1 2 3 4 5; 0 0 8 9 10; 0 0 0 14 15])

  @test_throws ErrorException AbstractAlgebra.Solve.solve(M, [ R(1) ])
  @test_throws ErrorException AbstractAlgebra.Solve.solve(M, [ R(1) ], side = :right)
  @test_throws ErrorException AbstractAlgebra.Solve.solve(M, matrix(R, 1, 1, [ R(1) ]))
  @test_throws ErrorException AbstractAlgebra.Solve.solve(M, matrix(R, 1, 1, [ R(1) ]), side = :right)
  @test_throws ArgumentError AbstractAlgebra.Solve.solve(M, [ R(1), R(2), R(3) ], side = :test)
  @test_throws ArgumentError AbstractAlgebra.Solve.solve(M, matrix(R, 3, 1, [ R(1), R(2), R(3) ]), side = :test)

  if R isa AbstractAlgebra.GFField
    RR = GF(2)
    @test_throws AssertionError solve(M, [ RR(1), RR(2), RR(3) ])
    @test_throws AssertionError solve(M, matrix(RR, 3, 1, [1, 2, 3]))
  end

  for b in [ [ R(1), R(2), R(3) ],
             matrix(R, 3, 1, [ R(1), R(2), R(3) ]),
             matrix(R, 3, 2, [ R(1), R(2), R(3), R(4), R(5), R(6) ]) ]
    @test @inferred AbstractAlgebra.Solve.can_solve(M, b, side = :right)
    x = @inferred AbstractAlgebra.Solve.solve(M, b, side = :right)
    @test M*x == b
    fl, x = @inferred AbstractAlgebra.Solve.can_solve_with_solution(M, b, side = :right)
    @test fl
    @test M*x == b
    fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(M, b, side = :right)
    @test fl
    @test M*x == b
    @test is_zero(M*K)
    @test ncols(K) == 2
    K = @inferred AbstractAlgebra.Solve.kernel(M, side = :right)
    @test is_zero(M*K)
    @test ncols(K) == 2
  end

  for b in [ [ R(1), R(1), R(1), R(1), R(1) ],
             matrix(R, 1, 5, [ R(1), R(1), R(1), R(1), R(1) ]),
             matrix(R, 2, 5, [ R(1), R(1), R(1), R(1), R(1),
                               R(1), R(1), R(1), R(1), R(1) ]) ]
    @test_throws ArgumentError AbstractAlgebra.Solve.solve(M, b)
    @test @inferred !AbstractAlgebra.Solve.can_solve(M, b)
    fl, x = @inferred AbstractAlgebra.Solve.can_solve_with_solution(M, b)
    @test !fl
    fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(M, b)
    @test !fl
  end

  for b in [ [ R(1), R(2), R(3), R(4), R(5) ],
             matrix(R, 1, 5, [ R(1), R(2), R(3), R(4), R(5)]),
             matrix(R, 2, 5, [ R(1), R(2), R(3), R(4), R(5),
                               R(0), R(0), R(8), R(9), R(10) ]) ]
    @test @inferred AbstractAlgebra.Solve.can_solve(M, b)
    x = @inferred AbstractAlgebra.Solve.solve(M, b)
    @test x*M == b
    fl, x = @inferred AbstractAlgebra.Solve.can_solve_with_solution(M, b)
    @test fl
    @test x*M == b
    fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(M, b)
    @test fl
    @test x*M == b
    @test is_zero(K*M)
    @test nrows(K) == 0
    K = @inferred AbstractAlgebra.Solve.kernel(M)
    @test is_zero(K*M)
    @test nrows(K) == 0
  end

  N = zero_matrix(R, 2, 1)
  b = [zero(R), zero(R)]
  fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(N, b, side = :right)
  @test fl
  @test N*x == b
  @test K == identity_matrix(R, 1)
  K = @inferred AbstractAlgebra.Solve.kernel(N, side = :right)
  @test K == identity_matrix(R, 1)

  N = zero_matrix(R, 1, 2)
  b = [zero(R)]
  fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(N, b, side = :right)
  @test fl
  @test N*x == b
  @test K == identity_matrix(R, 2) || K == swap_cols!(identity_matrix(R, 2), 1, 2)
  K = @inferred AbstractAlgebra.Solve.kernel(N, side = :right)
  @test K == identity_matrix(R, 2) || K == swap_cols!(identity_matrix(R, 2), 1, 2)
end

@testset "Linear solving context over $R with $NFTrait" for (R, NFTrait, is_default) in [
    (QQ, AbstractAlgebra.Solve.FFLUTrait, true),
    (ZZ, AbstractAlgebra.Solve.HermiteFormTrait, true),
    (GF(101), AbstractAlgebra.Solve.LUTrait, true),
    (GF(101), AbstractAlgebra.Solve.RREFTrait, false),
    (fraction_field(QQ["x"][1]), AbstractAlgebra.Solve.FFLUTrait, true),
    (residue_ring(ZZ, 16)[1], AbstractAlgebra.Solve.HowellFormTrait, true)
   ]
  M = matrix(R, [1 2 3 4 5; 0 0 8 9 10; 0 0 0 14 15])

  if is_default
    C = solve_init(M)
    @test AbstractAlgebra.Solve.matrix_normal_form_type(typeof(C)) === NFTrait()
    @test AbstractAlgebra.Solve.matrix_normal_form_type(C) === NFTrait()
    @test C isa AbstractAlgebra.solve_context_type(typeof(R))
    @test C isa AbstractAlgebra.solve_context_type(R)
    @test C isa AbstractAlgebra.solve_context_type(typeof(M))
    @test C isa AbstractAlgebra.solve_context_type(M)
  end

  C = solve_init(NFTrait(), M)

  @test AbstractAlgebra.Solve.matrix_normal_form_type(typeof(C)) === NFTrait()
  @test AbstractAlgebra.Solve.matrix_normal_form_type(C) === NFTrait()
  @test C isa AbstractAlgebra.solve_context_type(NFTrait(), elem_type(R))
  @test C isa AbstractAlgebra.solve_context_type(NFTrait(), R(1))
  @test C isa AbstractAlgebra.solve_context_type(NFTrait(), typeof(R))
  @test C isa AbstractAlgebra.solve_context_type(NFTrait(), R)
  @test C isa AbstractAlgebra.solve_context_type(NFTrait(), typeof(M))
  @test C isa AbstractAlgebra.solve_context_type(NFTrait(), M)

  @test_throws ErrorException AbstractAlgebra.Solve.solve(C, [ R(1) ])
  @test_throws ErrorException AbstractAlgebra.Solve.solve(C, [ R(1) ], side = :right)
  @test_throws ErrorException AbstractAlgebra.Solve.solve(C, matrix(R, 1, 1, [ R(1) ]))
  @test_throws ErrorException AbstractAlgebra.Solve.solve(C, matrix(R, 1, 1, [ R(1) ]), side = :right)
  @test_throws ArgumentError AbstractAlgebra.Solve.solve(C, [ R(1), R(2), R(3) ], side = :test)
  @test_throws ArgumentError AbstractAlgebra.Solve.solve(C, matrix(R, 3, 1, [ R(1), R(2), R(3) ]), side = :test)

  if R isa AbstractAlgebra.GFField
    RR = GF(2)
    @test_throws AssertionError solve(C, [ RR(1), RR(2), RR(3) ])
    @test_throws AssertionError solve(C, matrix(RR, 3, 1, [1, 2, 3]))
  end

  for b in [ [ R(1), R(2), R(3) ],
             matrix(R, 3, 1, [ R(1), R(2), R(3) ]),
             matrix(R, 3, 2, [ R(1), R(2), R(3), R(4), R(5), R(6) ]) ]
    @test @inferred AbstractAlgebra.Solve.can_solve(C, b, side = :right)
    x = @inferred AbstractAlgebra.Solve.solve(C, b, side = :right)
    @test M*x == b
    fl, x = @inferred AbstractAlgebra.Solve.can_solve_with_solution(C, b, side = :right)
    @test fl
    @test M*x == b
    fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(C, b, side = :right)
    @test fl
    @test M*x == b
    @test is_zero(M*K)
    @test ncols(K) == 2
    K = @inferred AbstractAlgebra.Solve.kernel(C, side = :right)
    @test is_zero(M*K)
    @test ncols(K) == 2
  end

  for b in [ [ R(1), R(1), R(1), R(1), R(1) ],
             matrix(R, 1, 5, [ R(1), R(1), R(1), R(1), R(1) ]),
             matrix(R, 2, 5, [ R(1), R(1), R(1), R(1), R(1),
                               R(1), R(1), R(1), R(1), R(1) ]) ]
    @test_throws ArgumentError AbstractAlgebra.Solve.solve(C, b)
    @test @inferred !AbstractAlgebra.Solve.can_solve(C, b)
    fl, x = @inferred AbstractAlgebra.Solve.can_solve_with_solution(C, b)
    @test !fl
    fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(C, b)
    @test !fl
  end

  for b in [ [ R(1), R(2), R(3), R(4), R(5) ],
             matrix(R, 1, 5, [ R(1), R(2), R(3), R(4), R(5)]),
             matrix(R, 2, 5, [ R(1), R(2), R(3), R(4), R(5),
                               R(0), R(0), R(8), R(9), R(10) ]) ]
    @test @inferred AbstractAlgebra.Solve.can_solve(C, b)
    x = @inferred AbstractAlgebra.Solve.solve(C, b)
    @test x*M == b
    fl, x = @inferred AbstractAlgebra.Solve.can_solve_with_solution(C, b)
    @test fl
    @test x*M == b
    fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(C, b)
    @test fl
    @test x*M == b
    @test is_zero(K*M)
    @test nrows(K) == 0
    K = @inferred AbstractAlgebra.Solve.kernel(C)
    @test is_zero(K*M)
    @test nrows(K) == 0
  end

  N = zero_matrix(R, 2, 1)
  C = AbstractAlgebra.Solve.solve_init(N)
  b = [zero(R), zero(R)]
  fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(C, b, side = :right)
  @test fl
  @test N*x == b
  @test K == identity_matrix(R, 1)
  K = @inferred AbstractAlgebra.Solve.kernel(C, side = :right)
  @test K == identity_matrix(R, 1)

  N = zero_matrix(R, 1, 2)
  C = AbstractAlgebra.Solve.solve_init(N)
  b = [zero(R)]
  fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(C, b, side = :right)
  @test fl
  @test N*x == b
  @test K == identity_matrix(R, 2) || K == swap_cols!(identity_matrix(R, 2), 1, 2)
  K = @inferred AbstractAlgebra.Solve.kernel(C, side = :right)
  @test K == identity_matrix(R, 2) || K == swap_cols!(identity_matrix(R, 2), 1, 2)
end

@testset "Linear solving context with non-default normal form" begin
  R = GF(101)
  M = matrix(R, [1 2 3 4 5; 0 0 8 9 10; 0 0 0 14 15])
  C = solve_init(M)
  @test_throws ErrorException solve(AbstractAlgebra.Solve.RREFTrait(), C, [R(1), R(2), R(3), R(4), R(5)])
  @test_throws ErrorException solve(AbstractAlgebra.Solve.RREFTrait(), C, matrix(R, 1, 5, [1, 2, 3, 4, 5]))

  # This should work because we use RREFTrait for the kernel if C uses LUTrait
  # for solving
  K = kernel(AbstractAlgebra.Solve.RREFTrait(), C)
  @test nrows(K) == 0

  # This should not work if C uses LUTrait
  @test_throws ErrorException kernel(AbstractAlgebra.Solve.HermiteFormTrait(), C)
end

@testset "Lazy transpose" begin
  M = matrix(QQ, [1 2 3 4 5; 0 0 8 9 10; 0 0 0 14 15])
  MT = AbstractAlgebra.Solve.lazy_transpose(M)

  @test M === @inferred AbstractAlgebra.Solve.data(MT)
  @test M === @inferred AbstractAlgebra.Solve.lazy_transpose(MT)
  @test transpose(M) == MT

  @test nrows(MT) == 5
  @test ncols(MT) == 3

  V = @inferred view(MT, 1, 1:3)
  @test V == [QQ(1), QQ(0), QQ(0)]
  V = @inferred view(MT, 1:5, 1)
  @test V == [QQ(1), QQ(2), QQ(3), QQ(4), QQ(5)]
  V = @inferred view(MT, 3:4, 1:2)
  @test V == QQ[3 8; 4 9]

  @test MT[2, 1] == QQ(2)
  MT[2, 1] = QQ(100)
  @test MT[2, 1] == QQ(100)
  @test M[1, 2] == QQ(100)

  @test base_ring(MT) == QQ

  @test (@inferred zero(MT)) == AbstractAlgebra.Solve.lazy_transpose(zero_matrix(QQ, 3, 5))
  @test (@inferred zero(MT, 2, 3)) == AbstractAlgebra.Solve.lazy_transpose(zero_matrix(QQ, 3, 2))

  S = @inferred similar(MT)
  @test S isa AbstractAlgebra.Solve.LazyTransposeMatElem
  @test nrows(S) == 5
  @test ncols(S) == 3
  @test base_ring(S) == QQ

  S = @inferred similar(MT, 2, 3)
  @test S isa AbstractAlgebra.Solve.LazyTransposeMatElem
  @test nrows(S) == 2
  @test ncols(S) == 3
  @test base_ring(S) == QQ
end

@testset "solve_triu" begin
  A = matrix(ZZ, 10, 10, [i<=j ? i+j-1 : 0 for i=1:10 for j=1:10])
  x = matrix(ZZ, rand(-10:10, 10, 10))
  @test AbstractAlgebra._solve_triu(A, A*x; side = :right) == x
  @test AbstractAlgebra._solve_triu(A, x*A; side = :left) == x

  A = matrix(ZZ, 20, 20, [i<=j ? i+j-1 : 0 for i=1:20 for j=1:20])
  x = matrix(ZZ, rand(-10:10, 20, 20))
  @test AbstractAlgebra.Strassen._solve_triu(A, A*x; cutoff = 10, side = :right) == x
  @test AbstractAlgebra.Strassen._solve_triu(A, x*A; cutoff = 10, side = :left) == x
end
