pretty_sort = AbstractAlgebra.pretty_sort

@testset "PrettyPrinting: pretty sorting things" begin
  v = [2, 3, 1]
  sort!(v; by = pretty_sort)
  @test v == [1, 2, 3]

  v = [[2, 1], [1, 3], [-1, 5], [-1]]
  sort!(v; by = pretty_sort)
  @test v == [[-1], [-1, 5], [1, 3], [2, 1]]

  v = [(2, 1), (1, 3), (-1, 5), (-1,)]
  sort!(v; by = pretty_sort)
  @test v == [(-1,), (-1, 5), (1, 3), (2, 1)]

  for F in (GF(5), QQ)
    Fx, x = F[:x]
    v = [F(2), F(3), F(1)]
    sort!(v; by = pretty_sort)
    @test v == [1, 2, 3]

    v = [F(0), x^2 + 2, x + 2, x + 1, Fx(1)]
    sort!(v; by = pretty_sort)
    @test v == [F(0), F(1), x + 1, x + 2, x^2 + 2]

    v = [1//x, 1//x^2, (x + 1)//x^2]
    sort!(v; by = pretty_sort)
    @test v ==  [1//x^2, (x + 1)//x^2, 1//x]
  end
end
