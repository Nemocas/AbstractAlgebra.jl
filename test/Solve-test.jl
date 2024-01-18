@testset "Linear solving" begin
  for R in [ QQ, ZZ ]
    M = matrix(R, [1 2 3 4 5; 0 0 8 9 10; 0 0 0 14 15])

    @test_throws ErrorException AbstractAlgebra.Solve.solve(M, [ R(1) ])
    @test_throws ErrorException AbstractAlgebra.Solve.solve(M, [ R(1) ], side = :left)
    @test_throws ErrorException AbstractAlgebra.Solve.solve(M, matrix(R, 1, 1, [ R(1) ]))
    @test_throws ArgumentError AbstractAlgebra.Solve.solve(M, [ R(1), R(2), R(3) ], side = :test)
    @test_throws ArgumentError AbstractAlgebra.Solve.solve(M, matrix(R, 3, 1, [ R(1), R(2), R(3) ]), side = :test)

    for b in [ [ R(1), R(2), R(3) ],
               matrix(R, 3, 1, [ R(1), R(2), R(3) ]),
               matrix(R, 3, 2, [ R(1), R(2), R(3), R(4), R(5), R(6) ]) ]
      @test @inferred AbstractAlgebra.Solve.can_solve(M, b)
      x = @inferred AbstractAlgebra.Solve.solve(M, b)
      @test M*x == b
      fl, x = @inferred AbstractAlgebra.Solve.can_solve_with_solution(M, b)
      @test fl
      @test M*x == b
      fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(M, b)
      @test fl
      @test M*x == b
      @test is_zero(M*K)
      @test ncols(K) == 2
    end

    for b in [ [ R(1), R(1), R(1), R(1), R(1) ],
               matrix(R, 1, 5, [ R(1), R(1), R(1), R(1), R(1) ]),
               matrix(R, 2, 5, [ R(1), R(1), R(1), R(1), R(1),
                                 R(1), R(1), R(1), R(1), R(1) ]) ]
      @test_throws ArgumentError AbstractAlgebra.Solve.solve(M, b, side = :left)
      @test @inferred !AbstractAlgebra.Solve.can_solve(M, b, side = :left)
      fl, x = @inferred AbstractAlgebra.Solve.can_solve_with_solution(M, b, side = :left)
      @test !fl
      fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(M, b, side = :left)
      @test !fl
    end

    for b in [ [ R(1), R(2), R(3), R(4), R(5) ],
               matrix(R, 1, 5, [ R(1), R(2), R(3), R(4), R(5)]),
               matrix(R, 2, 5, [ R(1), R(2), R(3), R(4), R(5),
                                 R(0), R(0), R(8), R(9), R(10) ]) ]
      @test @inferred AbstractAlgebra.Solve.can_solve(M, b, side = :left)
      x = @inferred AbstractAlgebra.Solve.solve(M, b, side = :left)
      @test x*M == b
      fl, x = @inferred AbstractAlgebra.Solve.can_solve_with_solution(M, b, side = :left)
      @test fl
      @test x*M == b
      fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(M, b, side = :left)
      @test fl
      @test x*M == b
      @test is_zero(K*M)
      @test nrows(K) == 0
    end

    N = zero_matrix(R, 2, 1)
    b = zeros(R, 2)
    fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(N, b)
    @test fl
    @test N*x == b
    @test K == identity_matrix(R, 1)

    N = zero_matrix(R, 1, 2)
    b = zeros(R, 1)
    fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(N, b)
    @test fl
    @test N*x == b
    @test K == identity_matrix(R, 2) || K == swap_cols!(identity_matrix(R, 2), 1, 2)
  end
end

@testset "Linear solving context" begin
  for R in [ QQ, ZZ ]
    M = matrix(R, [1 2 3 4 5; 0 0 8 9 10; 0 0 0 14 15])
    C = AbstractAlgebra.Solve.solve_init(M)

    @test_throws ErrorException AbstractAlgebra.Solve.solve(C, [ R(1) ])
    @test_throws ErrorException AbstractAlgebra.Solve.solve(C, [ R(1) ], side = :left)
    @test_throws ErrorException AbstractAlgebra.Solve.solve(C, matrix(R, 1, 1, [ R(1) ]))
    @test_throws ArgumentError AbstractAlgebra.Solve.solve(C, [ R(1), R(2), R(3) ], side = :test)
    @test_throws ArgumentError AbstractAlgebra.Solve.solve(C, matrix(R, 3, 1, [ R(1), R(2), R(3) ]), side = :test)

    for b in [ [ R(1), R(2), R(3) ],
               matrix(R, 3, 1, [ R(1), R(2), R(3) ]),
               matrix(R, 3, 2, [ R(1), R(2), R(3), R(4), R(5), R(6) ]) ]
      @test @inferred AbstractAlgebra.Solve.can_solve(C, b)
      x = @inferred AbstractAlgebra.Solve.solve(C, b)
      @test M*x == b
      fl, x = @inferred AbstractAlgebra.Solve.can_solve_with_solution(C, b)
      @test fl
      @test M*x == b
      fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(C, b)
      @test fl
      @test M*x == b
      @test is_zero(M*K)
      @test ncols(K) == 2
    end

    for b in [ [ R(1), R(1), R(1), R(1), R(1) ],
               matrix(R, 1, 5, [ R(1), R(1), R(1), R(1), R(1) ]),
               matrix(R, 2, 5, [ R(1), R(1), R(1), R(1), R(1),
                                 R(1), R(1), R(1), R(1), R(1) ]) ]
      @test_throws ArgumentError AbstractAlgebra.Solve.solve(C, b, side = :left)
      @test @inferred !AbstractAlgebra.Solve.can_solve(C, b, side = :left)
      fl, x = @inferred AbstractAlgebra.Solve.can_solve_with_solution(C, b, side = :left)
      @test !fl
      fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(C, b, side = :left)
      @test !fl
    end

    for b in [ [ R(1), R(2), R(3), R(4), R(5) ],
               matrix(R, 1, 5, [ R(1), R(2), R(3), R(4), R(5)]),
               matrix(R, 2, 5, [ R(1), R(2), R(3), R(4), R(5),
                                 R(0), R(0), R(8), R(9), R(10) ]) ]
      @test @inferred AbstractAlgebra.Solve.can_solve(C, b, side = :left)
      x = @inferred AbstractAlgebra.Solve.solve(C, b, side = :left)
      @test x*M == b
      fl, x = @inferred AbstractAlgebra.Solve.can_solve_with_solution(C, b, side = :left)
      @test fl
      @test x*M == b
      fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(C, b, side = :left)
      @test fl
      @test x*M == b
      @test is_zero(K*M)
      @test nrows(K) == 0
    end

    N = zero_matrix(R, 2, 1)
    C = AbstractAlgebra.Solve.solve_init(N)
    b = zeros(R, 2)
    fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(C, b)
    @test fl
    @test N*x == b
    @test K == identity_matrix(R, 1)

    N = zero_matrix(R, 1, 2)
    C = AbstractAlgebra.Solve.solve_init(N)
    b = zeros(R, 1)
    fl, x, K = @inferred AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(C, b)
    @test fl
    @test N*x == b
    @test K == identity_matrix(R, 2) || K == swap_cols!(identity_matrix(R, 2), 1, 2)
  end
end

@testset "Lazy transpose" begin
  M = matrix(QQ, [1 2 3 4 5; 0 0 8 9 10; 0 0 0 14 15])
  MT = AbstractAlgebra.Solve.lazy_transpose(M)

  @test M === @inferred AbstractAlgebra.Solve.data(MT)
  @test M === @inferred AbstractAlgebra.Solve.lazy_transpose(MT)
  @test transpose(M) == MT

  @test nrows(MT) == 5
  @test ncols(MT) == 3

  @test MT[2, 1] == QQ(2)
  MT[2, 1] = QQ(100)
  @test MT[2, 1] == QQ(100)
  @test M[1, 2] == QQ(100)

  @test base_ring(MT) == QQ

  @test @inferred zero(MT) == AbstractAlgebra.Solve.lazy_transpose(zero_matrix(QQ, 3, 5))
  @test @inferred zero(MT, 2, 3) == AbstractAlgebra.Solve.lazy_transpose(zero_matrix(QQ, 3, 2))

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
