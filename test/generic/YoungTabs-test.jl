function test_partition_type()
   print("youngtabs.partition_type...")

   @test Partition([4,3,1]) isa Generic.Partition
   @test Partition([4,3,1]) isa AbstractVector{Int}
   @test Partition([3,4,1]) isa Generic.Partition
   @test_throws ArgumentError Partition([4,3,0])

   p = Partition([4,3,1])

   @test size(p) == (3,)
   @test length(p) == 3

   @test_throws BoundsError p[0]
   @test p[1] == 4
   @test p[2] == 3
   @test p[3] == 1
   @test_throws BoundsError p[4]

   @test_throws ArgumentError p[1] = 2
   @test_throws ArgumentError p[2] = 5
   @test_throws ArgumentError p[2] = 0

   r = deepcopy(p)

   p[3] = 2
   @test p[3] == 2

   @test sum(p) == 9
   @test sum(r) == 8

   @test p == Partition([4,3,2])
   @test p != r

   println("PASS")
end

function test_partition_iter()
   print("youngtabs.partition_iter...")

   @test [Generic._numpart(i) for i in 0:10] == [1,1,2,3,5,7,11,15,22,30,42]
   @test Generic._numpart(100) == 190_569_292
   @test Generic._numpart(1000) == 24_061_467_864_032_622_473_692_149_727_991

   @test collect(AllParts(2)) == [Partition([1,1]), Partition([2])]
   @test collect(AllParts(1)) == [Partition([1])]
   @test collect(AllParts(0)) == [Partition(Int[])]

   println("PASS")
end

function test_ytabs_type()
   print("youngtabs.youngtableau_type...")

   lambda = [4,3,1]
   @test YoungTableau(Partition(lambda)) isa Generic.YoungTableau
   @test YoungTableau(Partition(lambda)) isa AbstractArray{Int, 2}
   @test YoungTableau(lambda) isa Generic.YoungTableau

   Y = YoungTableau(lambda)
   Z = YoungTableau(Partition(lambda))
   @test Y == Z
   @test length(Y) == 12
   @test size(Y) == (3,4)

   Y.part[2] = 2
   @test Y != Z

   @test YoungTableau([1,3,4,1]) == YoungTableau([4,3,1,1])

   Y = YoungTableau([4,3,1])

   @test Generic.inyoungtab((4,3), Y) == false
   @test Generic.inyoungtab((1,1), Y) == true
   @test Generic.inyoungtab((2,3), Y) == true
   @test Generic.inyoungtab((2,4), Y) == false
   @test Generic.inyoungtab((3,1), Y) == true
   @test Generic.inyoungtab((3,2), Y) == false

   @test Y[2,3] == 7
   @test Y[3,3] == 0
   @test Y[3,1] == 8
   @test_throws BoundsError Y[4,1]
   @test_throws BoundsError Y[1,5]

   for i in 1:sum(Y.part)
     @test Y[something(findfirst(isequal(i), Y), 0)] == i
   end

   println("PASS")
end

function test_conjugation()
   print("youngtabs.conjugation...")

   @test conj(Partition([1,1,1,1])) == Partition([4])
   @test conj(Partition([2,1,1])) == Partition([3,1])

   part = Partition([4,3,1])
   c, w = conj(part, collect(1:sum(part)))
   @test c == Partition([3,2,2,1])
   @test w == [1,5,8,2,6,3,7,4]

   part = Partition([3,1,1])
   c, w = conj(part, collect(1:sum(part)))
   @test c == Partition([3,1,1])
   @test w == [1,4,5,2,3]

   Y = YoungTableau([4,3,1])
   a,b = size(Y)
   @test size(conj(Y)) == (b,a)
   @test Y[1,2] == conj(Y)[2,1]
   @test Y[1,3] == conj(Y)[3,1]
   @test Y[2,2] == conj(Y)[2,2]

   @test Y[3,1] == conj(Y)[1,3]

   println("PASS")
end


function test_ytabs_dim()
   print("youngtabs.dim...")

   y = YoungTableau([5,4,3,3,1,1])
   # ┌───┬───┬───┬───┬───┐
   # │ 1 │ 2 │ 3 │ 4 │ 5 │
   # ├───┼───┼───┼───┼───┘
   # │ 6 │ 7 │ 8 │ 9 │
   # ├───┼───┼───┼───┘
   # │10 │11 │12 │
   # ├───┼───┼───┤
   # │13 │14 │15 │
   # ├───┼───┴───┘
   # │16 │
   # ├───┤
   # │17 │
   # └───┘

   @test Generic.collength(y, 1,1) == 5
   @test Generic.rowlength(y, 1,1) == 4
   @test hooklength(y, 1, 1) == 10

   @test Generic.collength(y, 1,2) == 3
   @test Generic.rowlength(y, 1,2) == 3
   @test hooklength(y, 1, 2) == 7

   @test Generic.collength(y, 2,3) == 2
   @test Generic.rowlength(y, 2,3) == 1
   @test hooklength(y, 2, 3) == 4

   @test Generic.collength(y, 4,2) == 0
   @test Generic.rowlength(y, 4,2) == 1
   @test hooklength(y, 4, 2) == 2

   @test Generic.collength(y, 4,3) == 0
   @test Generic.rowlength(y, 4,3) == 0
   @test hooklength(y, 4, 3) == 1

   @test Generic.collength(y, 4,4) == 0
   @test Generic.rowlength(y, 4,4) == 0
   @test hooklength(y, 4, 4) == 0

   Y = YoungTableau([5,3,3,1])

   @test [hooklength(Y,i,j) for i in 1:size(Y,1), j in 1:size(Y,2)] ==
      [8 6 5 2 1; 5 3 2 0 0; 4 2 1 0 0; 1 0 0 0 0]

   Y = YoungTableau([4,3,1])
   @test dim(Y) isa BigInt
   @test dim(Y) == 70

   @test dim(YoungTableau([2,2])) == 2
   @test dim(YoungTableau([3,1])) == 3
   @test dim(YoungTableau([2,1,1])) == 3
   @test dim(YoungTableau([1,1,1,1])) == 1
   @test dim(YoungTableau([4])) == 1

   @test dim(Int, YoungTableau([2,2])) == 2
   @test dim(Int, YoungTableau([3,1])) == 3
   @test dim(Int, YoungTableau([2,1,1])) == 3
   @test dim(Int, YoungTableau([1,1,1,1])) == 1
   @test dim(Int, YoungTableau([4])) == 1

   @test dim(YoungTableau([5,4,1])) == 288
   @test dim(YoungTableau([4,3,1,1])) == 216

   @test dim(YoungTableau(collect(10:-1:1))) == 44261486084874072183645699204710400

   println("PASS")
end

function test_skewdiags()
   print("youngtabs.skewdiags...")

   l = Partition([5,3,2,2,2,1,1])
   m = Partition([2,2,1])
   xi = l/m
   @test isa(xi, Generic.SkewDiagram)
   psi = SkewDiagram(l,m)
   @test isa(psi, Generic.SkewDiagram)
   @test xi == psi

   # 7×5 AbstractAlgebra.Generic.SkewDiagram:
   #  ⋅  ⋅  1  1  1
   #  ⋅  ⋅  1
   #  ⋅  1
   #  1  1
   #  1  1
   #  1
   #  1

   @test ( 1, 5) in xi
   @test ( 1, 3) in xi
   @test ( 2, 3) in xi
   @test !(( 2, 2) in xi)
   @test ( 5, 2) in xi
   @test !(( 6, 2) in xi)

   @test !(( 8, 9) in xi)
   @test !((-1, 1) in xi)
   @test !(( 1,-1) in xi)
   @test !((-1,-1) in xi)

   @test in((1, 4), SkewDiagram([4,3,1], [2]))

   @test matrix_repr(Partition([1], false)/Partition(Int[], false)) == ones(Int, 1,1)
   @test matrix_repr(xi) == [((i,j) in xi ? 1 : 0) for i in 1:size(xi.lam,1), j in 1:maximum(xi.lam)]

   @test has_left_neighbor(xi, 1, 5) == true
   @test has_left_neighbor(xi, 1, 3) == false
   @test has_left_neighbor(xi, 2, 3) == false
   @test has_left_neighbor(xi, 3, 3) == false
   @test has_left_neighbor(xi, 4, 2) == true
   @test has_left_neighbor(xi, 5, 2) == true
   @test has_left_neighbor(xi, 6, 1) == false
   @test has_left_neighbor(xi, 7, 1) == false

   @test has_bottom_neighbor(xi, 1, 5) == false
   @test has_bottom_neighbor(xi, 1, 3) == true
   @test has_bottom_neighbor(xi, 2, 3) == false
   @test has_bottom_neighbor(xi, 3, 3) == false
   @test has_bottom_neighbor(xi, 4, 2) == true
   @test has_bottom_neighbor(xi, 5, 2) == false
   @test has_bottom_neighbor(xi, 6, 1) == true
   @test has_bottom_neighbor(xi, 7, 1) == false

   println("PASS")
end

function test_rimhooks()
   print("youngtabs.rimhooks...")

   xi = Partition([2,1])/Partition(Int[], false)
   @test isrimhook(xi) == true

   xi = SkewDiagram([2,1], [1])
   @test isrimhook(xi) == false

   xi = Partition([2,2])/Partition(Int[], false)
   @test isrimhook(xi) == false

   xi = SkewDiagram([2,2], [1])
   @test isrimhook(xi) == true

   xi = SkewDiagram([4,3,1], [2,2])
   @test isrimhook(xi) == false

   xi = SkewDiagram([4,3,1], [2,1])
   @test isrimhook(xi) == false

   xi = SkewDiagram([4,3,1], [1])
   @test isrimhook(xi) == false

   xi = SkewDiagram([4,3,1], [2])
   @test isrimhook(xi) == true

   xi = Partition([4,3,2,1])/Partition([2,2,2,1])
   @test isrimhook(xi) == true

   xi = Partition([4,3,2,1])/Partition([3,3,2,1])
   @test isrimhook(xi) == true

   xi = Partition([4,3,2,1])/Partition([2,2,1,1])
   @test isrimhook(xi) == false

   xi = Partition([4,3,2,1])/Partition([3,2,2,1])
   @test isrimhook(xi) == false

   xi = Partition([4,3,2,1])/Partition([3,3,1,1])
   @test isrimhook(xi) == false

   xi = Partition([4,3,2,1])/Partition([4,3,2,1])
   @test isrimhook(xi) == true

   lambda = Partition([5,3,2,2,1])
   mu = Partition([2,2,1])

   @test isrimhook(lambda/mu) == false # is disconnected

   println("PASS")
end

function test_partitionseqs()
   print("youngtabs.partitionseqs...")
   @test partitionseq(Partition([1])) == BitVector([true, false])
   @test partitionseq([1]) == BitVector([true, false])
   @test partitionseq(Partition([1,1])) == BitVector([true, false, false])
   @test partitionseq(Partition([2,1])) == BitVector([true, false, true, false])

   t = partitionseq(Partition([5,4,2,1]))
   @test t == partitionseq([5,4,2,1])
   @test length(t) == 9
   @test t == BitVector([true, false, true, false, true, true, false, true, false])
   R = partitionseq(t)
   R[1] = false
   R[end] = true
   @test partitionseq(R) == BitVector([true,false,true,true,false])

   println("PASS")
end

function test_youngtabs()
   test_partition_type()
   test_partition_iter()
   test_ytabs_type()
   test_conjugation()
   test_ytabs_dim()
   test_skewdiags()
   test_rimhooks()
   test_partitionseqs()

   println("")
end
