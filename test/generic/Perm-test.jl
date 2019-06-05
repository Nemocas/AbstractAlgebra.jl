function test_perm_abstract_types()
   print("perm.abstract_types...")

   @test Generic.perm <: GroupElem

   @test Generic.PermGroup <: AbstractAlgebra.Group

   println("PASS")
end

function test_perm_constructors(types)
   print(rpad("perm.constructors...", 30))

   for T in types
      print("$T ")
      @test elem_type(Generic.PermGroup{T}) == Generic.perm{T}
      @test parent_type(Generic.perm{T}) == Generic.PermGroup{T}

      @test PermutationGroup(T(10)) isa Generic.PermGroup{T}
      G = PermutationGroup(T(10))
      @test elem_type(G) == Generic.perm{T}

      @test G() isa GroupElem
      @test G() isa Generic.perm{T}
      a = G()
      @test parent_type(typeof(a)) == Generic.PermGroup{T}
      @test parent(a) == G

      z = T[2, 3, 5, 4, 6, 7, 1, 9, 10, 8]

      @test G(z) isa GroupElem
      @test G(z) isa Generic.perm{T}
      b = G(z)
      @test typeof(b) == Generic.perm{T}
      @test parent_type(typeof(b)) == Generic.PermGroup{T}
      @test parent(b) == G

      @test rand(G) isa Generic.perm{T}
      g = rand(G)
      @test parent(g) == G
      @test parent(g) == PermutationGroup(T(10))

      @test convert(Vector{T}, G(z)) == z

      if T != Int
         @test parent(g) != PermutationGroup(10)
      end

   end

   println("PASS")
end

function test_perm_parsingGAP(types)
   print(rpad("perm.parsingGAP...", 30))

   @test Generic.parse_cycles("()") == (Int[], [1])
   @test Generic.parse_cycles("(1)(2)(3)") == ([1,2,3], [1,2,3,4])
   @test Generic.parse_cycles("Cycle decomposition: (1)(2,3)") == ([1,2,3], [1,2,4])
   @test Generic.parse_cycles("(1)(\n2, 3)") == ([1,2,3], [1,2,4])
   @test Generic.parse_cycles("(3,2,1)(4,5)") == ([3,2,1,4,5], [1,4,6])

   s = """
( 1, 22,73,64,78,81,  24 ,89,90,54,51,82,91,53, 18
 ,38,19,52,44,77,62,95,94,50,43,42,
10,67,87,60,36,12)(2,57,34,88)(3,92,76,17,99,96,30,55,45,41,98)(4,56,59,97,49,
21,15,9,26,86,83,29,27,66,6,58,28,5,68,40,72,7,84,93,39,79,23,46,63,32,61,100,
11)(8,80,71,75,35,14,85,25,20,70,65,16,48,47,37,74,33,13,31,69)
"""

s2 ="""
(1,22,73,64,78,81,24,89,90,54,51,82,91,53,18,38,19,52,44,77,62,95,94,50,43,42,\n10,67,87,60,36,12)(2,57,34,88)(3,92,76,17,99,96,30,55,45,41,98)(4,56,59,97,49,\n21,15,9,26,86,83,29,27,66,6,58,28,5,68,40,72,7,84,93,39,79,23,46,63,32,61,100,\n11)(8,80,71,75,35,14,85,25,20,70,65,16,48,47,37,74,33,13,31,69)
"""
   @test Generic.parse_cycles(s) == Generic.parse_cycles(s2)

   @test perm"(1,2,3)(5)(10)" isa GroupElem
   @test perm"(1,2,3)(5)(10)" isa Generic.perm
   @test perm"(1,2,3)(5)(10)" isa Generic.perm{Int}
   @test parent(perm"(1,2,3)(5)(10)") == PermGroup(10)
   @test parent(perm"(1,2,3)(5)(6)") == PermGroup(6)
   @test perm"(1,2,3,4,5)" == perm([2,3,4,5,1])
   @test perm"(3,2,1)(4,5)" == perm([3,1,2,5,4])

@test perm"""
( 1, 22,73,64,78,81,  24 ,89,90,54,51,82,91,53,18,38,19,52,44,77,62,95,94,50,43,42,
10,67,87,60,36,12)(2,57,34,88)(3,92,76,17,99,96,30,55,45,41,98)(4,56,59,97,49,
21,15,9,26,86,83,29,27,66,6,58,28,5,68,40,72,7,84,93,39,79,23,46,63,32,61,100,
11)(8,80,71,75,35,14,85,25,20,70,65,16,48,47,37,74,33,13,31,69)
""" == PermGroup(100)(s2)

   for T in types
      print("$T ")
      G = PermutationGroup(T(5))
      @test G("()") == G()
      @test G("()()") == G()
      @test G("(1)(2,3)") == perm(T[1,3,2,4,5])
      @test G("(2,3)") == G("(1)(2,3)")
      @test G("(3,2)") == G("(2,3)")
      @test G("(3,2,1)") == perm(T[3,1,2,4,5])
      @test G("(3,2,1)") == G("(1,3,2)") == G("(2,1,3)")
      @test G("(1,2,3,4,5)") == perm(T[2,3,4,5,1])
      @test G("(3,2,1)(4,5)") == perm(T[3,1,2,5,4])
      gg = G("(3,2,1)(4,5)")
      @test isdefined(gg, :cycles)

      @test_throws String G("(1,2)(3,4,5,6)")
      @test_throws String G("(1,2)(3,4,5,2)")

      @test all(elt == G(string(cycles(elt))) for elt in G)
   end
   println("PASS")
end

function test_perm_printing(types)
   print(rpad("perm.printing...", 30))

   for T in types
      print("$T ")
      G = PermutationGroup(T(10))

      b = G(T[2, 3, 5, 4, 6, 7, 1, 9, 10, 8])

      Generic.setpermstyle(:array);
      @test string(b) == "[2, 3, 5, 4, 6, 7, 1, 9, 10, 8]"

      Generic.setpermstyle(:cycles);
      @test string(b) == "(1,2,3,5,6,7)(8,9,10)"

      @test string(perm(T[1,2,3])) == "()"
      @test string(perm(T[3,2,1])) == "(1,3)"
   end

   println("PASS")
end

function test_perm_basic_manipulation(types)
   print(rpad("perm.basic_manipulation...", 30))

   for T in types
      print("$T ")
      G = PermutationGroup(T(10))

      a = G()
      b = deepcopy(a)
      c = G(T[2, 3, 5, 4, 6, 7, 1, 9, 10, 8])

      @test a == b

      @test parity(a) == 0
      @test parity(c) == 1

      @test length(unique([c,deepcopy(c)])) == 1

      @test setindex!(a, 5, 2) isa perm
      @test a[2] == T(5)

      a[1] = T(5)
      @test a[1] == T(5)
   end

   println("PASS")
end

function test_perm_iteration(types)
   print(rpad("perm.iteration...", 30))
   for T in types
      print("$T ")
      G = PermutationGroup(T(6))
      @test length(AllPerms(T(6))) == 720
      @test length(unique([deepcopy(p) for p in AllPerms(T(6))])) == 720
      @test order(G) isa BigInt
      @test order(T, G) isa promote_type(T, Int)
      @test order(G) == 720

      @test collect(G) isa Vector{perm{T}}

      elts = collect(G)

      @test elts[1] == G()
      @test length(elts) == 720
      @test length(unique(elts)) == 720
      @test length(G) == 720
      @test length(unique(collect(G))) == 720
   end

   println("PASS")
end

function test_perm_binary_ops(types)
   print(rpad("perm.binary_ops...", 30))

   for T in types
      print("$T ")
      G = PermutationGroup(T(3))

      a = perm(T[2,1,3]) # (1,2)
      b = perm(T[2,3,1]) # (1,2,3)

      @test a*b == G(T[3,2,1]) # (1,2)*(1,2,3) == (1,3)
      @test b*a == G(T[1,3,2]) # (1,2,3)*(1,2) == (2,3)
      @test a*a == G()
      @test b*b*b == G()

      # (1,2,3)*(2,3,4) == (1,3)(2,4)
      @test perm(T[2,3,1,4])*perm(T[1,3,4,2]) == perm(T[3,4,1,2])

      @test parity(G()) == 0
      p = parity(a)
      @test p == 1
      cycles(a)
      @test parity(a) == p

      z = G()

      for a in G, b in G
         @test parity(a*b) == (parity(b)+parity(a)) % 2
         @test AbstractAlgebra.mul!(z, a, b) == a*b
      end

      G = PermutationGroup(T(10))
      p = G(T[9,5,4,7,3,8,2,10,1,6]) # (1,9)(2,5,3,4,7)(6,8,10)
      @test p^0 == G()
      @test p^1 == p
      @test p^-1 == inv(p)
      @test p^5 == p*p*p*p*p
      @test p^-4 == inv(p)*inv(p)*inv(p)*inv(p)
      @test p^2 * p^-2 == G()

   end
   println("PASS")
end

function test_perm_mixed_binary_ops(types)
   print(rpad("perm.mixed_binary_ops...", 30))
      G = PermutationGroup(6)
      for T in types
         print("$T ")
         H = PermutationGroup(T(6))

         @test G(H()) == G()
         @test H(G()) == H()
         @test G() == H()
         @test rand(G)*rand(H) isa perm{promote_type(Int, T)}
         @test rand(H)*rand(G) isa perm{promote_type(Int, T)}
         @test G(rand(H)) isa perm{Int}
         @test H(rand(G)) isa perm{T}
      end

   println("PASS")
end

function test_perm_inversion(types)
   print(rpad("perm.inversion...", 30))
   for T in types
      print("$T ")
      G = PermutationGroup(T(10))

      a = G()
      b = G(T[2, 3, 5, 4, 6, 7, 1, 9, 10, 8])

      @test a == inv(a)
      @test a == b*inv(b)

      G = PermutationGroup(3)
      for a in G, b in G
         @test inv(a*b) == inv(b)*inv(a)
      end
   end
   println("PASS")
end

function test_misc_functions(types)
   print(rpad("perm.misc...", 30))

   for T in types
      print("$T ")
      G = PermutationGroup(T(10))
      a = G([2, 3, 5, 4, 6, 7, 1, 9, 10, 8])

      @test cycles(G()) isa Generic.CycleDec{T}
      @test collect(cycles(G())) == [T[i] for i in 1:10]
      @test order(G()) isa BigInt
      @test order(T, G()) isa promote_type(Int, T)
      @test order(G()) == 1

      @test collect(cycles(a)) == [T[1,2,3,5,6,7], T[4], T[8,9,10]]
      @test Generic.permtype(a) isa Vector{Int}
      @test Generic.permtype(a) == [6,3,1]
      @test cycles(a)[1] isa Vector{T}
      @test cycles(a)[1] == T[1,2,3,5,6,7]
      @test cycles(a)[2] == T[4]
      @test cycles(a)[3] == T[8,9,10]
      @test cycles(a)[1:3] == [T[1,2,3,5,6,7], T[4], T[8,9,10]]

      @test order(a) isa BigInt
      @test order(T, a) isa promote_type(Int, T)
      @test order(a) == 6
      @test a^6 == G()

      p = G([9,5,4,7,3,8,2,10,1,6])

      @test collect(cycles(p)) == [T[1,9],T[2,5,3,4,7],T[6,8,10]]
      @test Generic.permtype(p) == [5, 3, 2]
      @test order(p) == 30
      @test collect(cycles(p^2)) == [T[1], T[2,3,7,5,4], T[6,10,8], T[9]]
      @test Generic.permtype(p^2) == [5, 3, 1, 1]
      @test order(p^2) == 15
      @test collect(cycles(p^3)) == [T[1,9], T[2,4,5,7,3], T[6], T[8], T[10]]
      @test Generic.permtype(p^3) == [5, 2, 1, 1, 1]
      @test order(p^3) == 10
      @test collect(cycles(p^5)) == [T[1,9], T[2], T[3], T[4], T[5], T[6,10,8], T[7]]
      @test Generic.permtype(p^5) == [3, 2, 1, 1, 1, 1, 1]
      @test order(p^5) == 6

      @test matrix_repr(a) isa AbstractArray{T,2}
      @test matrix_repr(a) isa SparseMatrixCSC{T,T}
      M = matrix_repr(a)
      for (idx, val) in enumerate(a.d)
         @test M[idx, val] == 1
      end
   end

   println("PASS")
end

function test_characters(types)
   print(rpad("perm.characters...",30))

   for T in types
      print("$T ")
      N = T(7)
      G = PermutationGroup(N)

      @test all(character(p)(G()) == dim(YoungTableau(p)) for p in AllParts(N))

      N = T(3)
      G = PermutationGroup(N)
      ps = Partition.([[1,1,1], [2,1], [3]])
      l = Partition(T[1,1,1])

      @test typeof(character(l, ps[1])) == BigInt
      TT = (T !== BigInt ? Int : BigInt)
      @test typeof(character(T, l, rand(G))) == TT
      @test typeof(character(T, l, ps[1])) == TT

      @test [character(l, m) for m in ps] == [ 1,-1, 1]
      l = Partition([2,1])
      @test [character(l, m) for m in ps] == [ 2, 0,-1]
      l = Partition([3])
      @test [character(l, m) for m in ps] == [ 1, 1, 1]

      k_big = character(Partition([1,1,1]), G([2,3,1]))
      @test typeof(k_big) == BigInt

      k_int = character(Int, Partition([1,1,1]), G([2,3,1]))
      @test typeof(k_int) == Int

      @test k_big == k_int
   end

   @test character(Partition([2,2,2,2]), Partition([8])) == 0

   N = 4
   G = PermutationGroup(N)

   ps = Partition.([[1,1,1,1], [2,1,1], [2,2], [3,1], [4]])
   @test Set(AllParts(N)) == Set(ps)

   l = Partition([1,1,1,1])

   @test [character(l, m) for m in ps] == [ 1,-1, 1, 1,-1]
   l = Partition([2,1,1])
   @test [character(l, m) for m in ps] == [ 3,-1,-1, 0, 1]
   l = Partition([2,2])
   @test [character(l, m) for m in ps] == [ 2, 0, 2,-1, 0]
   l = Partition([3,1])
   @test [character(l, m) for m in ps] == [ 3, 1,-1, 0,-1]
   l = Partition([4])
   @test [character(l, m) for m in ps] == [ 1, 1, 1, 1, 1]

   # values taken from GAP; note that we specify the order of partitions to be
   # compatible with GAP numbering of conjugacy classes. This is NOT the order
   # of partitions given by AllParts.
   N = 5
   G = PermutationGroup(N)
   ps = Partition.([[1,1,1,1,1], [2,1,1,1], [2,2,1], [3,1,1], [3,2], [4,1], [5]])
   @test Set(AllParts(N)) == Set(ps)

   l = Partition([1,1,1,1,1])
   @test [character(l, m) for m in ps] == [   1,  -1,   1,   1,  -1,  -1,   1 ]
   l = Partition([2,1,1,1])
   @test [character(l, m) for m in ps] == [   4,  -2,   0,   1,   1,   0,  -1 ]
   l = Partition([2,2,1])
   @test [character(l, m) for m in ps] == [   5,  -1,   1,  -1,  -1,   1,   0 ]
   l = Partition([3,1,1])
   @test [character(l, m) for m in ps] == [   6,   0,  -2,   0,   0,   0,   1 ]
   l = Partition([3,2])
   @test [character(l, m) for m in ps] == [   5,   1,   1,  -1,   1,  -1,   0 ]
   l = Partition([4,1])
   @test [character(l, m) for m in ps] == [   4,   2,   0,   1,  -1,   0,  -1 ]
   l = Partition([5])
   @test [character(l, m) for m in ps] == [   1,   1,   1,   1,   1,   1,   1 ]

   # test for overflow

   p = Partition(collect(10:-1:1))
   @test character(p, PermutationGroup(big(55))()) == 44261486084874072183645699204710400
   @test dim(YoungTableau(p)) == 44261486084874072183645699204710400

   println("PASS")
end

function test_perm()
   IntTypes = [Int8, Int16, Int32, Int, UInt8, UInt16, UInt32, UInt, BigInt]
   test_perm_abstract_types()
   test_perm_constructors(IntTypes)
   test_perm_parsingGAP(IntTypes)
   test_perm_printing(IntTypes)
   test_perm_basic_manipulation(IntTypes)
   test_perm_iteration(IntTypes)
   test_perm_binary_ops(IntTypes)
   test_perm_mixed_binary_ops(IntTypes)
   test_perm_inversion(IntTypes)
   test_misc_functions(IntTypes)
   test_characters(IntTypes)
   println("")
end
