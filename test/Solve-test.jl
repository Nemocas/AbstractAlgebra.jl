@testset "Linear solving" begin
  M = matrix(QQ, [1 2 3 4 5; 0 0 8 9 10; 0 0 0 14 15])

  @test_throws ErrorException AbstractAlgebra.Solve.solve(M, [ QQ(1) ])
  @test_throws ErrorException AbstractAlgebra.Solve.solve(M, [ QQ(1) ], side = :left)
  @test_throws ErrorException AbstractAlgebra.Solve.solve(M, matrix(QQ, 1, 1, [ QQ(1) ]))
  @test_throws ArgumentError AbstractAlgebra.Solve.solve(M, [ QQ(1), QQ(2), QQ(3) ], side = :test)
  @test_throws ArgumentError AbstractAlgebra.Solve.solve(M, matrix(QQ, 3, 1, [ QQ(1), QQ(2), QQ(3) ]), side = :test)

  for b in [ [ QQ(1), QQ(2), QQ(3) ],
             matrix(QQ, 3, 1, [ QQ(1), QQ(2), QQ(3) ]),
             matrix(QQ, 3, 2, [ QQ(1), QQ(2), QQ(3), QQ(4), QQ(5), QQ(6) ]) ]
    @test AbstractAlgebra.Solve.can_solve(M, b)
    x = AbstractAlgebra.Solve.solve(M, b)
    @test M*x == b
    fl, x = AbstractAlgebra.Solve.can_solve_with_solution(M, b)
    @test fl
    @test M*x == b
    fl, x, K = AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(M, b)
    @test fl
    @test M*x == b
    @test is_zero(M*K)
    @test ncols(K) == 2
  end

  for b in [ [ QQ(1), QQ(1), QQ(1), QQ(1), QQ(1) ],
             matrix(QQ, 1, 5, [ QQ(1), QQ(1), QQ(1), QQ(1), QQ(1) ]),
             matrix(QQ, 2, 5, [ QQ(1), QQ(1), QQ(1), QQ(1), QQ(1),
                               QQ(1), QQ(1), QQ(1), QQ(1), QQ(1) ]) ]
    @test_throws ArgumentError AbstractAlgebra.Solve.solve(M, b, side = :left)
    @test !AbstractAlgebra.Solve.can_solve(M, b, side = :left)
    fl, x = AbstractAlgebra.Solve.can_solve_with_solution(M, b, side = :left)
    @test !fl
    fl, x, K = AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(M, b, side = :left)
    @test !fl
  end

  for b in [ [ QQ(1), QQ(2), QQ(3), QQ(4), QQ(5) ],
             matrix(QQ, 1, 5, [ QQ(1), QQ(2), QQ(3), QQ(4), QQ(5)]),
             matrix(QQ, 2, 5, [ QQ(1), QQ(2), QQ(3), QQ(4), QQ(5),
                               QQ(0), QQ(0), QQ(8), QQ(9), QQ(10) ]) ]
    @test AbstractAlgebra.Solve.can_solve(M, b, side = :left)
    x = AbstractAlgebra.Solve.solve(M, b, side = :left)
    @test x*M == b
    fl, x = AbstractAlgebra.Solve.can_solve_with_solution(M, b, side = :left)
    @test fl
    @test x*M == b
    fl, x, K = AbstractAlgebra.Solve.can_solve_with_solution_and_kernel(M, b, side = :left)
    @test fl
    @test x*M == b
    @test is_zero(K*M)
    @test nrows(K) == 0
  end
end
