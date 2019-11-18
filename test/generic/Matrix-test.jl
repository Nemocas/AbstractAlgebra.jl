primes100 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311,
313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401,
409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
499, 503, 509, 521, 523, 541]


function randprime(n::Int)
   if n > 100 || n < 1
      throw(DomainError())
   end
   return primes100[rand(1:n)]
end

function istriu(A::Generic.Mat)
   m = nrows(A)
   n = ncols(A)
   d = 0
   for c = 1:n
      for r = m:-1:1
         if !iszero(A[r,c])
            if r < d
               return false
            end
            d = r
            break
         end
      end
   end
   return true
end

function is_snf(A::Generic.Mat)
   m = nrows(A)
   n = ncols(A)
   a = A[1,1]
   for i = 2:min(m,n)
      q, r = divrem(A[i,i], a)
      if !iszero(r)
         return false
      end
      a = A[i,i]
   end
   for i = 1:n
      for j = 1:m
         if i == j
            continue
         end
         if !iszero(A[j,i])
            return false
         end
      end
   end
   return true
end

function is_weak_popov(P::Generic.Mat, rank::Int)
   zero_rows = 0
   pivots = zeros(ncols(P))
   for r = 1:nrows(P)
      p = AbstractAlgebra.find_pivot_popov(P, r)
      if P[r,p] == 0
         zero_rows += 1
         continue
      end
      if pivots[p] != 0
         return false
      end
      pivots[p] = r
   end
   if zero_rows != nrows(P)-rank
      return false
   end
   return true
end

# Simulate user matrix type belonging to AbstractArray
# with getindex but no setindex!
struct MyTestMatrix{T} <: AbstractArray{T, 2}
   d::T
   dim::Int
end

Base.getindex(a::MyTestMatrix{T}, r::Int, c::Int) where T = a.d

Base.size(a::MyTestMatrix{T}) where T = a.dim, a.dim

@testset "Generic.Mat.constructors..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   @test elem_type(S) == Generic.MatSpaceElem{elem_type(R)}
   @test elem_type(Generic.MatSpace{elem_type(R)}) == Generic.MatSpaceElem{elem_type(R)}
   @test parent_type(Generic.MatSpaceElem{elem_type(R)}) == Generic.MatSpace{elem_type(R)}
   @test base_ring(S) == R
   @test nrows(S) == 3
   @test ncols(S) == 3

   @test typeof(S) <: Generic.MatSpace

   f = S(t^2 + 1)

   @test isa(f, MatElem)
   @test parent_type(f) == typeof(S)

   g = S(2)

   @test isa(g, MatElem)

   h = S(BigInt(23))

   @test isa(h, MatElem)

   k = S([t t + 2 t^2 + 3t + 1; 2t R(2) t + 1; t^2 + 2 t + 1 R(0)])

   @test isa(k, MatElem)

   l = S(k)

   @test isa(l, MatElem)

   m = S()

   @test isa(m, MatElem)

   @test_throws ErrorConstrDimMismatch S([t t^2 ; t^3 t^4])
   @test_throws ErrorConstrDimMismatch S([t t^2 t^3 ; t^4 t^5 t^6 ; t^7 t^8 t^9 ; t t^2 t^3])
   @test_throws ErrorConstrDimMismatch S([t, t^2])
   @test_throws ErrorConstrDimMismatch S([t, t^2, t^3, t^4, t^5, t^6, t^7, t^8, t^9, t^10])

   # test literal construction
   for T = (R, ZZ)
      m = T[1 2]
      @test m isa MatElem{elem_type(T)}
      @test size(m) == (1, 2)
      @test m[1, 1] == T(1)
      @test m[1, 2] == T(2)

      m = T[1; 2; 3]
      @test m isa MatElem{elem_type(T)}
      @test size(m) == (3, 1)
      @test m[1, 1] == T(1)
      @test m[2, 1] == T(2)
      @test m[3, 1] == T(3)

      m = T[1 2; 3 4; 5 6]
      @test m isa MatElem{elem_type(T)}
      @test size(m) == (3, 2)
      @test m[1, 1] == T(1)
      @test m[1, 2] == T(2)
      @test m[2, 1] == T(3)
      @test m[2, 2] == T(4)
      @test m[3, 1] == T(5)
      @test m[3, 2] == T(6)

      @test_throws ArgumentError T[1; 2 3]
   end

   arr = [1 2; 3 4]
   arr2 = [1, 2, 3, 4, 5, 6]

   for T in [R, Int, BigInt, Rational{Int}, Rational{BigInt}]
      for M in (matrix(R, map(T, arr)),
                matrix(R, 2, 2, map(T, arr)))
         @test isa(M, Generic.MatSpaceElem{elem_type(R)})
         @test M.base_ring == R
         @test nrows(M) == 2
         @test ncols(M) == 2
      end
      @test_throws ErrorConstrDimMismatch matrix(R, 2, 3, map(T, arr))

      M2 = matrix(R, 2, 3, map(T, arr2))
      @test isa(M2, Generic.MatSpaceElem{elem_type(R)})
      @test M2.base_ring == R
      @test nrows(M2) == 2
      @test ncols(M2) == 3
      @test_throws ErrorConstrDimMismatch matrix(R, 2, 2, map(T, arr2))
      @test_throws ErrorConstrDimMismatch matrix(R, 2, 4, map(T, arr2))
   end

   M = matrix(R, arr')
   @test isa(M, Generic.MatSpaceElem{elem_type(R)})

   M = matrix(R, 2, 3, view(arr2', 1:6))
   @test isa(M, Generic.MatSpaceElem{elem_type(R)})

   M3 = zero_matrix(R, 2, 3)

   @test isa(M3, Generic.MatSpaceElem{elem_type(R)})
   @test M3.base_ring == R

   M4 = identity_matrix(R, 3)

   @test isa(M4, Generic.MatSpaceElem{elem_type(R)})
   @test M4.base_ring == R

   M5 = identity_matrix(R, 2, 3)
   M6 = identity_matrix(M5)

   @test isa(M5, Generic.MatSpaceElem{elem_type(R)})
   @test M5.base_ring == R
   @test M5 == M6

   M7 = identity_matrix(R, 3, 2)
   M8 = identity_matrix(M7)

   @test isa(M7, Generic.MatSpaceElem{elem_type(R)})
   @test M7.base_ring == R
   @test M7 == M8

   x = zero_matrix(R, 2, 2)
   y = zero_matrix(ZZ, 2, 3)

   @test x in [x, y]
   @test x in [y, x]
   @test !(x in [y])

   @test x in keys(Dict(x => 1))
   @test !(y in keys(Dict(x => 1)))

   # Test creation from AbstractArray without setindex!
   A = MyTestMatrix(BigInt(3), 2)
   S = MatrixSpace(ZZ, 2, 2)

   @test isa(S(A), Generic.MatSpaceElem{BigInt})

   # Test original matrix not modified
   S = MatrixSpace(QQ, 2, 2)
   a = Rational{BigInt}(1)
   A = [a a; a a]
   M = S(A)

   @test A[1, 1] === a

   a = BigInt(1)
   A = [a a; a a]
   M = S(A)

   @test A[1, 1] === a
end

@testset "Generic.Mat.size/axes..." begin
   A = matrix(QQ, [1 2 3; 4 5 6; 7 8 9])
   B = matrix(QQ, [1 2 3 4; 5 6 7 8])

   @test size(A) == (3,3)
   @test size(A, 1) == 3
   @test size(A, 2) == 3
   @test size(A, rand(3:99)) == 1
   @test_throws BoundsError size(A, 0)
   @test_throws BoundsError size(A, -rand(1:99))

   @test axes(A) == (1:3, 1:3)
   @test axes(A, 1) == 1:3
   @test axes(A, 2) == 1:3
   @test axes(A, rand(3:99)) == 1:1
   @test_throws BoundsError axes(A, 0)
   @test_throws BoundsError axes(A, -rand(1:99))

   @test issquare(A)

   @test size(B) == (2,4)
   @test size(B, 1) == 2
   @test size(B, 2) == 4
   @test size(B, rand(3:99)) == 1
   @test_throws BoundsError size(B, 0)
   @test_throws BoundsError size(B, -rand(1:99))

   @test axes(B) == (1:2, 1:4)
   @test axes(B, 1) == 1:2
   @test axes(B, 2) == 1:4
   @test axes(B, rand(3:99)) == 1:1
   @test_throws BoundsError axes(A, 0)
   @test_throws BoundsError axes(A, -rand(1:99))

   @test !issquare(B)
end

@testset "Generic.Mat.manipulation..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])

   @test dense_matrix_type(R) == elem_type(S)

   let ET = AbstractAlgebra.Generic.Poly{Rational{BigInt}}
      @test eltype(typeof(A)) == eltype(A) == elem_type(base_ring(A)) == ET
      @test eltype(typeof(B)) == eltype(B) == elem_type(base_ring(B)) == ET
   end

   @test iszero(zero(S))
   @test isone(one(S))

   @test zero(A) == zero(S)
   @test one(A) == one(S)

   B[1, 1] = R(3)
   @test B[1, 1] == R(3)

   B[1, 1] = 4
   @test B[1, 1] == R(4)

   B[1, 1] = BigInt(5)
   @test B[1, 1] == R(5)

   @test nrows(B) == 3
   @test ncols(B) == 3

   let AA = deepcopy(A)
      @test AA == A
      @test AA[1,1] !== A[1,1]
   end

   let AA = copy(A)
      @test AA == A
      @test AA[1,1] === A[1,1]
   end

   C = S([t + 1 R(0) R(1); t^2 R(0) t; R(0) R(0) R(0)])

   @test iszero_row(C, 3)
   @test !iszero_row(C, 1)
   @test iszero_column(C, 2)
   @test !iszero_column(C, 1)

   @test length(A) == length(B) == length(C) == 9
   @test !any(isempty, (A, B, C))

   @test length(matrix(R, zeros(Int, 2, 3))) == 6

   n = matrix(R, zeros(Int, 3, 2))
   @test length(n) == 6
   let ET = AbstractAlgebra.Generic.Poly{Rational{BigInt}}
      @test eltype(typeof(n)) == ET
      @test eltype(n) == ET
      @test elem_type(base_ring(n)) == ET
   end

   for n = (matrix(R, zeros(Int, 2, 0)),
            matrix(R, zeros(Int, 0, 2)))
      @test length(n) == 0
      @test isempty(n)
   end

   M3 = MatrixAlgebra(R, 3)
   for m3 in [rand(M3, 0:9, -9:9),
              rand(rng, M3, 0:9, -9:9),
              randmat_triu(M3, 0:9, -9:9),
              randmat_triu(rng, M3, 0:9, -9:9),
              randmat_with_rank(M3, 2, 0:9, -9:9),
              randmat_with_rank(rng, M3, 2, 0:9, -9:9)]
      @test length(m3) == 9
      ET = AbstractAlgebra.Generic.Poly{Rational{BigInt}}
      @test eltype(typeof(m3)) == ET
      @test eltype(m3) == ET
      @test elem_type(base_ring(m3)) == ET
      @test !isempty(m3)
      @test !iszero(m3)
      @test m3 isa Generic.MatAlgElem
      @test parent(m3) == M3
   end

   M0 = MatrixAlgebra(R, 0)
   m0 = rand(M0, 0:9, -9, 9)
   @test length(m0) == 0
   @test isempty(m0)

   M45 = MatrixSpace(R, 4, 5)
   for m45 in [rand(M45, 0:9, -9:9),
               rand(rng, M45, 0:9, -9:9),
               randmat_triu(M45, 0:9, -9:9),
               randmat_triu(rng, M45, 0:9, -9:9),
               randmat_with_rank(M45, 3, 0:9, -9:9),
               randmat_with_rank(rng, M45, 3, 0:9, -9:9)]
      @test length(m45) == 20
      @test !iszero(m45)
      @test m45 isa Generic.MatSpaceElem
      @test parent(m45) == M45
      @test_throws DomainError one(m45)
   end

   @test_throws DomainError one(M45)

   let m = matrix(ZZ, 2, 3, 1:6)
      @test typeof(m[1, 1]) == BigInt # not in AbstractAlgebra's hierarchy
      @test eltype(m) == BigInt
      @test eltype(typeof(m)) == BigInt
   end
end

@testset "Generic.Mat.unary_ops..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([-t - 1 (-t) -R(1); -t^2 (-t) (-t); -R(-2) (-t - 2) (-t^2 - t - 1)])

   @test -A == B

   # Exact ring
   S = MatrixSpace(ZZ, rand(0:9), rand(0:9))
   A = rand(S, -1000:1000)

   @test iszero(A + (-A))
   @test A == -(-A)
   @test -A == S(-A.entries)

   # Exact field
   S = MatrixSpace(GF(7), rand(0:9), rand(0:9))
   A = rand(S)

   @test iszero(A + (-A))
   @test A == -(-A)
   @test -A == S(-A.entries)

   # Inexact ring
   S = MatrixSpace(RealField["t"][1], rand(0:9), rand(0:9))
   A = rand(S, 0:200, -1000:1000)

   @test iszero(A + (-A))
   @test A == -(-A)
   @test -A == S(-A.entries)

   # Inexact field
   S = MatrixSpace(RealField, rand(0:9), rand(0:9))
   A = rand(S, -1000:1000)

   @test iszero(A + (-A))
   @test A == -(-A)
   @test -A == S(-A.entries)

   # Non-integral domain
   S = MatrixSpace(ResidueRing(ZZ, 6), rand(0:9), rand(0:9))
   A = rand(S, 0:5)

   @test iszero(A + (-A))
   @test A == -(-A)
   @test -A == S(-A.entries)
end

@testset "Generic.Mat.getindex..." begin
   S = MatrixSpace(ZZ, 3, 3)

   A = S([1 2 3; 4 5 6; 7 8 9])

   B = @inferred A[1:2, 1:2]

   @test typeof(B) == typeof(A)
   @test B == MatrixSpace(ZZ, 2, 2)([1 2; 4 5])

   B[1, 1] = 10
   @test A == S([1 2 3; 4 5 6; 7 8 9])

   C = @inferred B[1:2, 1:2]

   @test typeof(C) == typeof(A)
   @test C == MatrixSpace(ZZ, 2, 2)([10 2; 4 5])

   C[1, 1] = 20
   @test B == MatrixSpace(ZZ, 2, 2)([10 2; 4 5])
   @test A == S([1 2 3; 4 5 6; 7 8 9])

   D = @inferred A[:, 2:3]

   @test D == matrix(ZZ, 3, 2, [2, 3, 5, 6, 8, 9])

   @test A == @inferred A[:, :]
   @test B == @inferred B[:, :]
   @test C == @inferred C[:, :]
   @test D == @inferred D[:, :]

   # bounds check
   S = MatrixSpace(ZZ, rand(1:9), 0)
   A = S()
   @test isempty(A)
   rows = UnitRange(extrema(rand(axes(A)[1], 2))...)
   # A[rows, :] must be a valid indexing
   @test size(A[rows, :]) == (length(rows), 0)
   @test_throws BoundsError A[1:10, :]

   S = MatrixSpace(ZZ, 0, rand(1:9))
   A = S()
   @test isempty(A)
   cols = UnitRange(extrema(rand(axes(A)[2], 2))...)
   # A[:, cols] must be a valid indexing
   @test size(A[:, cols]) == (0, length(cols))
   @test_throws BoundsError A[:, 1:10]

   S = MatrixSpace(ZZ, 0, 0)
   A = S()
   @test isempty(A)
   # A[:, :] must be a valid indexing
   @test size(A[:, :]) == (0, 0)
   @test_throws BoundsError A[2:3, 1:10]

   # Exact ring
   R = ZZ
   S = MatrixSpace(R, rand(1:9), rand(1:9))

   A = rand(S, -1000:1000)
   ((i, j), (k, l)) = extrema.(rand.(axes(A), 2))

   @test A[i:j, k:l] == matrix(R, A.entries[i:j, k:l])
   @test A[:, k:l] == matrix(R, A.entries[:, k:l])
   @test A[i:j, :] == matrix(R, A.entries[i:j, :])   
   @test A[[i:j;], [k:l;]] == matrix(R, A.entries[[i:j;], [k:l;]])
   @test A[i:2:j, k:2:l] == matrix(R, A.entries[i:2:j, k:2:l])
   
   rows, cols = randsubseq.(axes(A), rand(2))
   @test A[rows, cols] == matrix(R, A.entries[rows, cols])

   # Exact field
   R = GF(7)
   S = MatrixSpace(R, rand(1:9), rand(1:9))

   A = rand(S)
   ((i, j), (k, l)) = extrema.(rand.(axes(A), 2))

   @test A[i:j, k:l] == matrix(R, A.entries[i:j, k:l])
   @test A[:, k:l] == matrix(R, A.entries[:, k:l])
   @test A[i:j, :] == matrix(R, A.entries[i:j, :])   
   @test A[[i:j;], [k:l;]] == matrix(R, A.entries[[i:j;], [k:l;]])
   @test A[i:2:j, k:2:l] == matrix(R, A.entries[i:2:j, k:2:l])
   
   rows, cols = randsubseq.(axes(A), rand(2))
   @test A[rows, cols] == matrix(R, A.entries[rows, cols])

   # Inexact ring
   R = RealField["t"][1]
   S = MatrixSpace(R, rand(1:9), rand(1:9))

   A = rand(S, 0:200, -1000:1000)
   ((i, j), (k, l)) = extrema.(rand.(axes(A), 2))

   @test A[i:j, k:l] == matrix(R, A.entries[i:j, k:l])
   @test A[:, k:l] == matrix(R, A.entries[:, k:l])
   @test A[i:j, :] == matrix(R, A.entries[i:j, :])   
   @test A[[i:j;], [k:l;]] == matrix(R, A.entries[[i:j;], [k:l;]])
   @test A[i:2:j, k:2:l] == matrix(R, A.entries[i:2:j, k:2:l])
   
   rows, cols = randsubseq.(axes(A), rand(2))
   @test A[rows, cols] == matrix(R, A.entries[rows, cols])

   # Inexact field
   R = RealField
   S = MatrixSpace(R, rand(1:9), rand(1:9))

   A = rand(S, -1000:1000)
   ((i, j), (k, l)) = extrema.(rand.(axes(A), 2))

   @test A[i:j, k:l] == matrix(R, A.entries[i:j, k:l])
   @test A[:, k:l] == matrix(R, A.entries[:, k:l])
   @test A[i:j, :] == matrix(R, A.entries[i:j, :])   
   @test A[[i:j;], [k:l;]] == matrix(R, A.entries[[i:j;], [k:l;]])
   @test A[i:2:j, k:2:l] == matrix(R, A.entries[i:2:j, k:2:l])
   
   rows, cols = randsubseq.(axes(A), rand(2))
   @test A[rows, cols] == matrix(R, A.entries[rows, cols])

   # Non-integral domain
   R = ResidueRing(ZZ, 6)
   S = MatrixSpace(R, rand(1:9), rand(1:9))

   A = rand(S, 0:5)
   ((i, j), (k, l)) = extrema.(rand.(axes(A), 2))

   @test A[i:j, k:l] == matrix(R, A.entries[i:j, k:l])
   @test A[:, k:l] == matrix(R, A.entries[:, k:l])
   @test A[i:j, :] == matrix(R, A.entries[i:j, :])   
   @test A[[i:j;], [k:l;]] == matrix(R, A.entries[[i:j;], [k:l;]])
   @test A[i:2:j, k:2:l] == matrix(R, A.entries[i:2:j, k:2:l])
   
   rows, cols = randsubseq.(axes(A), rand(2))
   @test A[rows, cols] == matrix(R, A.entries[rows, cols])

   # Fraction field
   R = QQ
   S = MatrixSpace(R, rand(1:9), rand(1:9))

   A = rand(S, -1000:1000)
   ((i, j), (k, l)) = extrema.(rand.(axes(A), 2))

   @test A[i:j, k:l] == matrix(R, A.entries[i:j, k:l])
   @test A[:, k:l] == matrix(R, A.entries[:, k:l])
   @test A[i:j, :] == matrix(R, A.entries[i:j, :])   
   @test A[[i:j;], [k:l;]] == matrix(R, A.entries[[i:j;], [k:l;]])
   @test A[i:2:j, k:2:l] == matrix(R, A.entries[i:2:j, k:2:l])
   
   rows, cols = randsubseq.(axes(A), rand(2))
   @test A[rows, cols] == matrix(R, A.entries[rows, cols])
end

@testset "Generic.Mat.block_replacement..." begin
   S = MatrixSpace(ZZ, 9, 9)
   (r, c) = (rand(1:9), rand(1:9))
   T = MatrixSpace(ZZ, r, c)
   a = rand(S, -100:100)
   b = rand(T, -100:100)
   startr = rand(1:(9-r+1))
   endr = startr + r - 1
   startc = rand(1:(9-c+1))
   endc = startc + c - 1
   a[startr:endr, startc:endc] = b
   @test a[startr:endr, startc:endc] == b
end

@testset "Generic.Mat.binary_ops..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])

   @test A + B == S([t+3 t+3 R(2); t^2 + t 2*t+1 2*t+2; R(-3) t^2 + t + 2 t^3 + 1*t^2 + t + 1])

   @test A - B == S([t-1 t-3 R(0); t^2 - t R(-1) R(-2); R(-1) (-t^2 + t + 2) (-t^3 + t^2 + t + 1)])

   @test A*B == S([t^2 + 2*t + 1 2*t^2 + 4*t + 3 t^3 + t^2 + 3*t + 1; 3*t^2 - t (t^3 + 4*t^2 + t) t^4 + 2*t^2 + 2*t; t-5 t^4 + t^3 + 2*t^2 + 3*t - 4 t^5 + 1*t^4 + t^3 + t^2 + 4*t + 2])

   # Exact ring
   R = ZZ

   for S in (MatrixSpace(R, rand(1:9), rand(1:9)),
             let n = rand(1:9)
                MatrixSpace(R, n, n)
             end)

      A = rand(S, -1000:1000)
      B = rand(S, -1000:1000)

      @test A + B == S(A.entries + B.entries)
      @test A - B == S(A.entries - B.entries)
      @test A + B == A - (-B)
      if issquare(A)
         @test A * B == S(A.entries * B.entries)
      end
   end

   # Exact field
   R = GF(7)

   for S in (MatrixSpace(R, rand(1:9), rand(1:9)),
             let n = rand(1:9)
                MatrixSpace(R, n, n)
             end)

      A = rand(S)
      B = rand(S)

      @test A + B == S(A.entries + B.entries)
      @test A - B == S(A.entries - B.entries)
      @test A + B == A - (-B)
      if issquare(A)
         @test A * B == S(A.entries * B.entries)
      end
   end

   # Inexact ring
   R = RealField["t"][1]

   for S in (MatrixSpace(R, rand(1:9), rand(1:9)),
             let n = rand(1:9)
                MatrixSpace(R, n, n)
             end)

      A = rand(S, 0:20, -100:100)
      B = rand(S, 0:20, -100:100)

      @test A + B == S(A.entries + B.entries)
      @test A - B == S(A.entries - B.entries)
      @test A + B == A - (-B)
      if issquare(A)
         @test A * B == S(A.entries * B.entries)
      end
   end

   # Inexact field
   R = RealField

   for S in (MatrixSpace(R, rand(1:9), rand(1:9)),
             let n = rand(1:9)
                MatrixSpace(R, n, n)
             end)

      A = rand(S, -1000:1000)
      B = rand(S, -1000:1000)

      @test A + B == S(A.entries + B.entries)
      @test A - B == S(A.entries - B.entries)
      @test A + B == A - (-B)
      if issquare(A)
         @test A * B == S(A.entries * B.entries)
      end
   end

   # Non-integral domain
   R = ResidueRing(ZZ, 6)

   for S in (MatrixSpace(R, rand(1:9), rand(1:9)),
             let n = rand(1:9)
                MatrixSpace(R, n, n)
             end)

      A = rand(S, 0:5)
      B = rand(S, 0:5)

      @test A + B == S(A.entries + B.entries)
      @test A - B == S(A.entries - B.entries)
      @test A + B == A - (-B)
      if issquare(A)
         @test A * B == S(A.entries * B.entries)
      end
   end

   # Fraction field
   R = QQ

   for S in (MatrixSpace(R, rand(1:9), rand(1:9)),
             let n = rand(1:9)
                MatrixSpace(R, n, n)
             end)

      A = rand(S, -1000:1000)
      B = rand(S, -1000:1000)

      @test A + B == S(A.entries + B.entries)
      @test A - B == S(A.entries - B.entries)
      @test A + B == A - (-B)
      if issquare(A)
         @test A * B == S(A.entries * B.entries)
      end
   end
end

# add x to all the elements of the main diagonal of a copy of M
add_diag(M::Matrix, x) = [i != j ? M[i, j] : M[i, j] + x for (i, j) in Tuple.(CartesianIndices(M))]

@testset "Generic.Mat.adhoc_binary..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test 12 + A == A + 12
   @test BigInt(11) + A == A + BigInt(11)
   @test Rational{BigInt}(11) + A == A + Rational{BigInt}(11)
   @test (t + 1) + A == A + (t + 1)
   @test A - (t + 1) == -((t + 1) - A)
   @test A - 3 == -(3 - A)
   @test A - BigInt(7) == -(BigInt(7) - A)
   @test A - Rational{BigInt}(7) == -(Rational{BigInt}(7) - A)
   @test 3*A == A*3
   @test BigInt(3)*A == A*BigInt(3)
   @test Rational{BigInt}(3)*A == A*Rational{BigInt}(3)
   @test (t - 1)*A == A*(t - 1)

   # Exact ring
   R = ZZ
   S = MatrixSpace(R, rand(1:9), rand(1:9))

   A = rand(S, -1000:1000)

   for t in Any[rand(-1000:1000), big(rand(-1000:1000)), rand(R, -1000:1000)]
      @test A + t == t + A
      @test A + t == S(add_diag(A.entries, t))
      @test A - t == -(t - A)
      @test A - t == S(add_diag(A.entries, -t))
      @test A * t == t * A
      @test A * t == S(A.entries .* t)
   end

   # Exact field
   R = GF(7)
   S = MatrixSpace(R, rand(1:9), rand(1:9))

   A = rand(S)

   for t in Any[rand(-1000:1000), big(rand(-1000:1000)), rand(R)]
      @test A + t == t + A
      @test A + t == S(add_diag(A.entries, t))
      @test A - t == -(t - A)
      @test A - t == S(add_diag(A.entries, -t))
      @test A * t == t * A
      @test A * t == S(A.entries .* t)
   end

   # Inexact ring
   R = RealField["t"][1]
   S = MatrixSpace(R, rand(1:9), rand(1:9))

   A = rand(S, 0:200, -1000:1000)

   for t in Any[rand(-1000:1000), big(rand(-1000:1000)), rand(R, 0:200, -1000:1000)]
      @test A + t == t + A
      @test A + t == S(add_diag(A.entries, t))
      @test A - t == -(t - A)
      @test A - t == S(add_diag(A.entries, -t))
      @test A * t == t * A
      @test A * t == S(A.entries .* t)
   end

   # Inexact field
   R = RealField
   S = MatrixSpace(R, rand(1:9), rand(1:9))

   A = rand(S, -1000:1000)

   for t in Any[rand(-1000:1000), big(rand(-1000:1000)), rand(R, -1000:1000)]
      @test A + t == t + A
      @test A + t == S(add_diag(A.entries, t))
      @test A - t == -(t - A)
      @test A - t == S(add_diag(A.entries, -t))
      @test A * t == t * A
      @test A * t == S(A.entries .* t)
   end

   # Non-integral domain
   R = ResidueRing(ZZ, 6)
   S = MatrixSpace(R, rand(1:9), rand(1:9))

   A = rand(S, 0:5)

   for t in Any[rand(-1000:1000), big(rand(-1000:1000)), rand(R, 0:5)]
      @test A + t == t + A
      @test A + t == S(add_diag(A.entries, t))
      @test A - t == -(t - A)
      @test A - t == S(add_diag(A.entries, -t))
      @test A * t == t * A
      @test A * t == S(A.entries .* t)
   end

   # Fraction field
   R = QQ
   S = MatrixSpace(R, rand(1:9), rand(1:9))

   A = rand(S, -1000:1000)

   for t in Any[rand(-1000:1000), big(rand(-1000:1000)), rand(R, -1000:1000)]
      @test A + t == t + A
      @test A + t == S(add_diag(A.entries, t))
      @test A - t == -(t - A)
      @test A - t == S(add_diag(A.entries, -t))
      @test A * t == t * A
      @test A * t == S(A.entries .* t)
   end
end

@testset "Generic.Mat.permutation..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   T = PermutationGroup(3)
   P = T([2, 3, 1])

   @test A == inv(P)*(P*A)
end

@testset "Generic.Mat.comparison..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test A == B

   @test A != one(S)
end

@testset "Generic.Mat.adhoc_comparison..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test S(12) == 12
   @test S(5) == BigInt(5)
   @test S(5) == Rational{BigInt}(5)
   @test S(t + 1) == t + 1
   @test 12 == S(12)
   @test BigInt(5) == S(5)
   @test Rational{BigInt}(5) == S(5)
   @test t + 1 == S(t + 1)
   @test A != one(S)
   @test one(S) == one(S)
end

@testset "Generic.Mat.powering..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test A^5 == A^2*A^3

   @test A^0 == one(S)

   S = MatrixSpace(QQ, 2, 2)

   A = S(Rational{BigInt}[2 3; 7 -4])

   @test A^-1 == inv(A)
end

@testset "Generic.Mat.adhoc_exact_division..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test divexact(5*A, 5) == A
   @test divexact(12*A, BigInt(12)) == A
   @test divexact(12*A, Rational{BigInt}(12)) == A
   @test divexact((1 + t)*A, 1 + t) == A
end

@testset "Generic.Mat.transpose..." begin
   R, t = PolynomialRing(QQ, "t")
   arr = [t + 1 t R(1); t^2 t t]
   A = matrix(R, arr)
   B = matrix(R, permutedims(arr, [2, 1]))
   @test transpose(A) == B

   arr = [t + 1 t; t^2 t]
   A = matrix(R, arr)
   B = matrix(R, permutedims(arr, [2, 1]))
   @test transpose(A) == B
end

@testset "Generic.Mat.gram..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test gram(A) == S([2*t^2 + 2*t + 2 t^3 + 2*t^2 + t 2*t^2 + t - 1; t^3 + 2*t^2 + t t^4 + 2*t^2 t^3 + 3*t; 2*t^2 + t - 1 t^3 + 3*t t^4 + 2*t^3 + 4*t^2 + 6*t + 9])
end

@testset "Generic.Mat.tr..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test tr(A) == t^2 + 3t + 2
end

@testset "Generic.Mat.content..." begin
   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test content((1 + t)*A) == 1 + t
end

@testset "Generic.Mat.lu..." begin
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   S = MatrixSpace(K, 3, 3)

   A = S([a + 1 2a + 3 a^2 + 1; 2a^2 - 1 a - 1 2a; a^2 + 3a + 1 2a K(1)])

   r, P, L, U = lu(A)

   @test r == 3
   @test P*A == L*U

   A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])

   r, P, L, U = lu(A)

   @test r == 3
   @test P*A == L*U

   A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 - 2 a - 1 2a])

   r, P, L, U = lu(A)

   @test r == 2
   @test P*A == L*U

   R, z = PolynomialRing(ZZ, "z")
   F = FractionField(R)

   A = matrix(F, 3, 3, [0, 0, 11, 78*z^3-102*z^2+48*z+12, 92, -16*z^2+80*z-149, -377*z^3+493*z^2-232*z-58, -448, 80*z^2-385*z+719])

   r, P, L, U = lu(A)

   @test r == 3
   @test P*A == L*U
end

@testset "Generic.Mat.fflu..." begin
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   S = MatrixSpace(K, 3, 3)

   A = S([a + 1 2a + 3 a^2 + 1; 2a^2 - 1 a - 1 2a; a^2 + 3a + 1 2a K(1)])

   r, d, P, L, U = fflu(A)

   D = S()
   D[1, 1] = inv(U[1, 1])
   D[2, 2] = inv(U[1, 1]*U[2, 2])
   D[3, 3] = inv(U[2, 2])

   @test r == 3
   @test P*A == L*D*U

   A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])

   r, d, P, L, U = fflu(A)

   D = S()
   D[1, 1] = inv(U[1, 1])
   D[2, 2] = inv(U[1, 1]*U[2, 2])
   D[3, 3] = inv(U[2, 2])

   @test r == 3
   @test P*A == L*D*U

   A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 - 2 a - 1 2a])

   r, d, P, L, U = fflu(A)

   D = S()
   D[1, 1] = inv(U[1, 1])
   D[2, 2] = inv(U[1, 1]*U[2, 2])
   D[3, 3] = inv(U[2, 2])

   @test r == 2
   @test P*A == L*D*U

   A = matrix(QQ, 3, 3, [0, 0, 1, 12, 1, 11, 1, 0, 1])

   r, d, P, L, U, = fflu(A)

   D = zero_matrix(QQ, 3, 3)
   D[1, 1] = inv(U[1, 1])
   D[2, 2] = inv(U[1, 1]*U[2, 2])
   D[3, 3] = inv(U[2, 2])
   @test r == 3
   @test P*A == L*D*U
end

@testset "Generic.Mat.det..." begin
   S, x = PolynomialRing(ResidueRing(ZZ, 1009*2003), "x")

   for dim = 0:5
      R = MatrixSpace(S, dim, dim)

      M = rand(R, 0:5, -100:100)

      @test det(M) == AbstractAlgebra.det_clow(M)
   end

   S, z = PolynomialRing(ZZ, "z")

   for dim = 0:5
      R = MatrixSpace(S, dim, dim)

      M = rand(R, 0:3, -20:20)

      @test det(M) == AbstractAlgebra.det_clow(M)
   end

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")

   for dim = 0:7
      S = MatrixSpace(K, dim, dim)

      M = rand(S, 0:2, -100:100)

      @test det(M) == AbstractAlgebra.det_clow(M)
   end

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   for dim = 0:5
      T = MatrixSpace(S, dim, dim)
      M = rand(T, 0:2, 0:2, -10:10)

      @test det(M) == AbstractAlgebra.det_clow(M)
   end
end

@testset "Generic.Mat.minors..." begin
   S, z = PolynomialRing(ZZ,"z")
   n = 5
   R = MatrixSpace(S,n,n)
   for r = 0:n
      M = randmat_with_rank(R, r, 0:3, 0:3)
      @test [1] == minors(M, 0)
      for i = r + 1:n
         for m in minors(M, i)
            @test m == 0
         end
      end
      @test [det(M)] == minors(M, n)
      @test [] == minors(M, n + 1)
   end
end

@testset "Generic.Mat.rank..." begin
   S = ResidueRing(ZZ, 20011*10007)
   R = MatrixSpace(S, 5, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, -100:100)
      do_test = false
      r = 0

      try
         r = rank(M)
         do_test = true
      catch e
         if !(e isa ErrorException)
            rethrow(e)
         end
      end

      if do_test
         @test r == i
      end
   end

   S, z = PolynomialRing(ZZ, "z")
   R = MatrixSpace(S, 4, 4)

   M = R([S(-2) S(0) S(5) S(3); 5*z^2+5*z-5 S(0) S(-z^2+z) 5*z^2+5*z+1; 2*z-1 S(0) z^2+3*z+2 S(-4*z); 3*z-5 S(0) S(-5*z+5) S(1)])

   @test rank(M) == 3

   R = MatrixSpace(S, 5, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, 0:3, -20:20)

      @test rank(M) == i
   end

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   S = MatrixSpace(K, 3, 3)

   M = S([a a^2 + 2*a - 1 2*a^2 - 1*a; 2*a+2 2*a^2 + 2*a (-2*a^2 - 2*a); (-a) (-a^2) a^2])

   @test rank(M) == 2

   S = MatrixSpace(K, 5, 5)

   for i = 0:5
      M = randmat_with_rank(S, i, 0:2, -100:100)

      @test rank(M) == i
   end

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   T = MatrixSpace(S, 3, 3)

   M = T([(2*x^2)*y^2+(-2*x^2-2*x)*y+(-x^2+2*x) S(0) (-x^2-2)*y^2+(x^2+2*x+2)*y+(2*x^2-x-1);
    (-x)*y^2+(-x^2+x-1)*y+(x^2-2*x+2) S(0) (2*x^2+x-1)*y^2+(-2*x^2-2*x-2)*y+(x^2-x);
    (-x+2)*y^2+(x^2+x+1)*y+(-x^2+x-1) S(0) (-x^2-x+2)*y^2+(-x-1)*y+(-x-1)])

   @test rank(M) == 2

   T = MatrixSpace(S, 5, 5)

   for i = 0:5
      M = randmat_with_rank(T, i, 0:2, 0:2, -20:20)

      @test rank(M) == i
   end
end

@testset "Generic.Mat.solve_lu..." begin
   S = QQ

   for dim = 0:5
      R = MatrixSpace(S, dim, dim)
      U = MatrixSpace(S, dim, rand(1:5))

      M = randmat_with_rank(R, dim, -100:100)
      b = rand(U, -100:100)

      x = Generic.solve_lu(M, b)

      @test M*x == b
   end

   S, y = PolynomialRing(ZZ, "y")
   K = FractionField(S)

   for dim = 0:5
      R = MatrixSpace(S, dim, dim)
      U = MatrixSpace(S, dim, rand(1:5))

      M = randmat_with_rank(R, dim, 0:5, -100:100)
      b = rand(U, 0:5, -100:100);

      MK = matrix(K, elem_type(K)[ K(M[i, j]) for i in 1:nrows(M), j in 1:ncols(M) ])
      bK = matrix(K, elem_type(K)[ K(b[i, j]) for i in 1:nrows(b), j in 1:ncols(b) ])

      x = Generic.solve_lu(MK, bK)

      @test MK*x == bK
   end
end

@testset "Generic.Mat.solve_rational..." begin
   S = ResidueRing(ZZ, 20011*10007)

   for dim = 0:5
      R = MatrixSpace(S, dim, dim)
      U = MatrixSpace(S, dim, rand(1:5))

      M = randmat_with_rank(R, dim, -100:100)
      b = rand(U, -100:100)

      do_test = false
      try
         x, d = solve_rational(M, b)
         do_test = true
      catch e
         if !(e isa ErrorException)
             rethrow(e)
         end
      end

      if do_test
         @test M*x == d*b
      end
   end

   S, z = PolynomialRing(ZZ, "z")

   for dim = 0:5
      R = MatrixSpace(S, dim, dim)
      U = MatrixSpace(S, dim, rand(1:5))

      M = randmat_with_rank(R, dim, 0:3, -20:20)
      b = rand(U, 0:3, -20:20);

      x, d = solve_rational(M, b)

      @test M*x == d*b
   end

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")

   for dim = 0:5
      S = MatrixSpace(K, dim, dim)
      U = MatrixSpace(K, dim, rand(1:5))

      M = randmat_with_rank(S, dim, 0:2, -100:100)
      b = rand(U, 0:2, -100:100)

      x = solve(M, b)

      @test M*x == b
   end

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   for dim = 0:5
      T = MatrixSpace(S, dim, dim)
      U = MatrixSpace(S, dim, rand(1:5))

      M = randmat_with_rank(T, dim, 0:2, 0:2, -20:20)
      b = rand(U, 0:2, 0:2, -20:20)

      x, d = solve_rational(M, b)

      @test M*x == d*b
   end

   R, t = PolynomialRing(AbstractAlgebra.JuliaQQ, "t")
   K, a = NumberField(t^3 + 3t + 1, "a")
   S, y = PolynomialRing(K, "y")
   T = MatrixSpace(S, 3, 3)
   U = MatrixSpace(S, 3, 1)

   M = T([3y*a^2 + (y + 1)*a + 2y (5y+1)*a^2 + 2a + y - 1 a^2 + (-a) + 2y; (y + 1)*a^2 + 2y - 4 3y*a^2 + (2y - 1)*a + y (4y - 1)*a^2 + (y - 1)*a + 5; 2a + y + 1 (2y + 2)*a^2 + 3y*a + 3y a^2 + (-y-1)*a + (-y - 3)])
   b = U(permutedims([4y*a^2 + 4y*a + 2y + 1 5y*a^2 + (2y + 1)*a + 6y + 1 (y + 1)*a^2 + 3y*a + 2y + 4], [2, 1]))

   x, d = solve_rational(M, b)

   @test M*x == d*b
end

@testset "Generic.Mat.solve_left..." begin
   for R in [ZZ, QQ]
      for iter = 1:40
         for dim = 0:5
            r = rand(1:5)
            n = rand(1:5)
            c = rand(1:5)

            S = MatrixSpace(R, r, n)
            U = MatrixSpace(R, n, c)

            X1 = rand(S, -20:20)
            M = rand(U, -20:20)

            B = X1*M
            X = solve_left(M, X1*M)

            @test X*M == B
         end
      end
   end

   R, x = PolynomialRing(QQ, "x")

   for iter = 1:4
      for dim = 0:5
         r = rand(1:5)
         n = rand(1:5)
         c = rand(1:5)

         S = MatrixSpace(R, r, n)
         U = MatrixSpace(R, n, c)

         X1 = rand(S, 1:2, -10:10)
         M = rand(U, 1:2, -10:10)

         B = X1*M
         X = solve_left(M, X1*M)

         @test X*M == B
      end
   end
end

@testset "Generic.Mat.solve_triu..." begin
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")

   for dim = 0:10
      S = MatrixSpace(K, dim, dim)
      U = MatrixSpace(K, dim, rand(1:5))

      M = randmat_triu(S, 0:2, -100:100)
      b = rand(U, 0:2, -100:100)

      x = solve_triu(M, b, false)

      @test M*x == b
   end
end

@testset "Generic.Mat.solve_left_reduced_triu..." begin
   for iter = 1:40
      n = rand(1:6)
      m = rand(1:n)
      S = MatrixSpace(ZZ, m, n)
      U = MatrixSpace(ZZ, 1, n)

      M = randmat_with_rank(S, rand(1:m), -20:20)
      r = rand(U, -20:20)

      M = hnf(M)

      flag, x = can_solve_left_reduced_triu(r, M)

      @test flag == false || x*M == r
   end
end

@testset "Generic.Mat.rref..." begin
   S = ResidueRing(ZZ, 20011*10007)
   R = MatrixSpace(S, 5, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, -100:100)

      do_test = false
      r = 0
      d = S(0)
      A = M

      try
          r, d, A = rref(M)
          do_test = true
      catch e
         if !(e isa ErrorException)
            rethrow(e)
         end
      end

      if do_test
         @test r == i
         @test isrref(A)
      end
   end

   S, z = PolynomialRing(ZZ, "z")
   R = MatrixSpace(S, 5, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, 0:3, -20:20)

      r, d, A = rref(M)

      @test r == i
      @test isrref(A)
   end

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   S = MatrixSpace(K, 5, 5)

   for i = 0:5
      M = randmat_with_rank(S, i, 0:2, -100:100)

      r, A = rref(M)

      @test r == i
      @test isrref(A)
   end

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   T = MatrixSpace(S, 5, 5)

   for i = 0:5
      M = randmat_with_rank(T, i, 0:2, 0:2, -20:20)

      r, d, A = rref(M)

      @test r == i
      @test isrref(A)
   end
end

@testset "Generic.Mat.nullspace..." begin
   S = ResidueRing(ZZ, 20011*10007)
   R = MatrixSpace(S, 5, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, -100:100)

      do_test = false
      n = 0
      N = M
      r = 0

      try
         n, N = nullspace(M)
         r = rank(N)
         do_test = true
      catch e
         if !(e isa ErrorException)
            rethrow(e)
         end
      end

      if do_test
         @test n == 5 - i
         @test r == n
         @test iszero(M*N)
      end
   end

   S, z = PolynomialRing(ZZ, "z")
   R = MatrixSpace(S, 5, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, 0:3, -20:20)

      n, N = nullspace(M)

      @test n == 5 - i
      @test rank(N) == n
      @test iszero(M*N)
   end

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   S = MatrixSpace(K, 5, 5)

   for i = 0:5
      M = randmat_with_rank(S, i, 0:2, -100:100)

      n, N = nullspace(M)

      @test n == 5 - i
      @test rank(N) == n
      @test iszero(M*N)
   end

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   T = MatrixSpace(S, 5, 5)

   for i = 0:5
      M = randmat_with_rank(T, i, 0:2, 0:2, -20:20)

      n, N = nullspace(M)

      @test n == 5 - i
      @test rank(N) == n
      @test iszero(M*N)
   end
end

@testset "Generic.Mat.kernel..." begin
   R = MatrixSpace(ZZ, 5, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, -20:20)

      n, N = kernel(M)

      @test n == 5 - i
      @test rank(N) == n
      @test iszero(M*N)

      n, N = left_kernel(M)

      @test n == 5 - i
      @test rank(N) == n
      @test iszero(N*M)
   end

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   S = MatrixSpace(K, 5, 5)

   for i = 0:5
      M = randmat_with_rank(S, i, 0:2, -100:100)

      n, N = kernel(M)

      @test n == 5 - i
      @test rank(N) == n
      @test iszero(M*N)

      n, N = left_kernel(M)

      @test n == 5 - i
      @test rank(N) == n
      @test iszero(N*M)
   end

   R, x = PolynomialRing(QQ, "x")
   T = MatrixSpace(R, 5, 5)

   for i = 0:5
      M = randmat_with_rank(T, i, 0:2, -20:20)

      n, N = kernel(M)

      @test n == 5 - i
      @test rank(N) == n
      @test iszero(M*N)

      n, N = left_kernel(M)

      @test n == 5 - i
      @test rank(N) == n
      @test iszero(N*M)
   end
end

@testset "Generic.Mat.inversion..." begin
   for dim = 2:5
      R = MatrixSpace(ZZ, dim, dim)
      M = R(1)
      i = rand(1:dim-1)
      j = rand(i+1:dim)
      M[i,j] = 1 # E_{i,j} elementary matrix

      N, c = pseudo_inv(M)
      @test N isa elem_type(R)
      @test c isa eltype(M)

      @test isunit(c)
      @test N[i,j] == -1
      @test M*N == N*M == c*R(1)

      M[j,i] = -1
      NN, cc = pseudo_inv(M)
      @test NN[i,j] == -1
      @test NN[j,i] == 1

      @test M*NN == NN*M == cc*R(1)
   end

   S = ResidueRing(ZZ, 20011*10007)

   for dim = 1:5
      R = MatrixSpace(S, dim, dim)

      M = randmat_with_rank(R, dim, -100:100)

      do_test = false
      X = M
      d = R(0)

      try
          X, d = pseudo_inv(M)
          do_test = true
      catch e
         if !(e isa ErrorException)
            rethrow(e)
         end
      end

      if do_test
         @test M*X == d*one(R)
      end
   end

   S, z = PolynomialRing(ZZ, "z")

   for dim = 1:5
      R = MatrixSpace(S, dim, dim)

      M = randmat_with_rank(R, dim, 0:3, -20:20)

      X, d = pseudo_inv(M)

      @test M*X == d*one(R)
   end

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")

   for dim = 1:5
      S = MatrixSpace(K, dim, dim)

      M = randmat_with_rank(S, dim, 0:2, -100:100)

      X, d = pseudo_inv(M)

      @test M*X == d*one(S)
   end

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   for dim = 1:5
      T = MatrixSpace(S, dim, dim)

      M = randmat_with_rank(T, dim, 0:2, 0:2, -20:20)

      X, d = pseudo_inv(M)

      @test M*X == d*one(T)
   end
end

@testset "Generic.Mat.hessenberg..." begin
   R = ResidueRing(ZZ, 18446744073709551629)

   for dim = 0:5
      S = MatrixSpace(R, dim, dim)
      U, x = PolynomialRing(R, "x")

      for i = 1:10
         M = rand(S, -5:5)

         A = hessenberg(M)

         @test ishessenberg(A)
      end
   end
end

@testset "Generic.Mat.kronecker_product..." begin
   R = ResidueRing(ZZ, 18446744073709551629)
   S = MatrixSpace(R, 2, 3)
   S2 = MatrixSpace(R, 2, 2)
   S3 = MatrixSpace(R, 3, 3)

   A = S(R.([2 3 5; 9 6 3]))
   B = S2(R.([2 3; 1 4]))
   C = S3(R.([2 3 5; 1 4 7; 9 6 3]))

   @test size(kronecker_product(A, A)) == (4,9)
   @test kronecker_product(B*A,A*C) == kronecker_product(B,A) * kronecker_product(A,C)
end

@testset "Generic.Mat.charpoly..." begin
   R = ResidueRing(ZZ, 18446744073709551629)

   for dim = 0:5
      S = MatrixSpace(R, dim, dim)
      U, x = PolynomialRing(R, "x")

      for i = 1:10
         M = rand(S, -5:5)

         p1 = charpoly(U, M)
         p2 = charpoly_danilevsky!(U, M)

         @test p1 == p2
      end

      for i = 1:10
         M = rand(S, -5:5)

         p1 = charpoly(U, M)
         p2 = charpoly_danilevsky_ff!(U, M)

         @test p1 == p2
      end

      for i = 1:10
         M = rand(S, -5:5)

         p1 = charpoly(U, M)
         p2 = charpoly_hessenberg!(U, M)

         @test p1 == p2
      end
   end

   R, x = PolynomialRing(ZZ, "x")
   U, z = PolynomialRing(R, "z")
   T = MatrixSpace(R, 6, 6)

   M = T()
   for i = 1:3
      for j = 1:3
         M[i, j] = rand(R, 0:2, -10:10)
         M[i + 3, j + 3] = deepcopy(M[i, j])
      end
   end

   p1 = charpoly(U, M)

   for i = 1:10
      similarity!(M, rand(1:6), R(rand(R, 0:2, -3:3)))
   end

   p2 = charpoly(U, M)

   @test p1 == p2
end

@testset "Generic.Mat.minpoly..." begin
   R = GF(103)
   T, y = PolynomialRing(R, "y")

   M = R[92 97 8;
          0 5 13;
          0 16 2]

   @test minpoly(T, M) == y^2+96*y+8

   R = GF(3)
   T, y = PolynomialRing(R, "y")

   M = R[1 2 0 2;
         1 2 1 0;
         1 2 2 1;
         2 1 2 0]

   @test minpoly(T, M) == y^2 + 2y

   R = GF(13)
   T, y = PolynomialRing(R, "y")

   M = R[7 6 1;
         7 7 5;
         8 12 5]

   @test minpoly(T, M) == y^2+10*y

   M = R[4 0 9 5;
         1 0 1 9;
         0 0 7 6;
         0 0 3 10]

   @test minpoly(T, M) == y^2 + 9y

   M = R[2 7 0 0 0 0;
         1 0 0 0 0 0;
         0 0 2 7 0 0;
         0 0 1 0 0 0;
         0 0 0 0 4 3;
         0 0 0 0 1 0]

   @test minpoly(T, M) == (y^2+9*y+10)*(y^2+11*y+6)

   M = R[2 7 0 0 0 0;
         1 0 1 0 0 0;
         0 0 2 7 0 0;
         0 0 1 0 0 0;
         0 0 0 0 4 3;
         0 0 0 0 1 0]

   @test minpoly(T, M) == (y^2+9*y+10)*(y^2+11*y+6)^2

   S = MatrixSpace(R, 1, 1)
   M = S()

   @test minpoly(T, M) == y

   S = MatrixSpace(R, 0, 0)
   M = S()

   @test minpoly(T, M) == 1

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   U, z = PolynomialRing(S, "z")
   T = MatrixSpace(S, 6, 6)

   M = T()
   for i = 1:3
      for j = 1:3
         M[i, j] = rand(S, 0:3, 0:3, -10:10)
         M[i + 3, j + 3] = deepcopy(M[i, j])
      end
   end

   f = minpoly(U, M)

   @test degree(f) <= 3

   R, x = PolynomialRing(ZZ, "x")
   U, z = PolynomialRing(R, "z")
   T = MatrixSpace(R, 6, 6)

   M = T()
   for i = 1:3
      for j = 1:3
         M[i, j] = rand(R, 0:2, -10:10)
         M[i + 3, j + 3] = deepcopy(M[i, j])
      end
   end

   p1 = minpoly(U, M)

   for i = 1:10
      similarity!(M, rand(1:6), R(rand(R, 0:2, -3:3)))
   end

   p2 = minpoly(U, M)

   @test p1 == p2
end

@testset "Generic.Mat.row_col_swapping..." begin
   R, x = PolynomialRing(ZZ, "x")
   M = MatrixSpace(R, 3, 2)

   a = M(map(R, [1 2; 3 4; 5 6]))

   @test swap_rows(a, 1, 3) == M(map(R, [5 6; 3 4; 1 2]))

   swap_rows!(a, 2, 3)

   @test a == M(map(R, [1 2; 5 6; 3 4]))

   @test swap_cols(a, 1, 2) == matrix(R, [2 1; 6 5; 4 3])

   swap_cols!(a, 2, 1)

   @test a == matrix(R, [2 1; 6 5; 4 3])

   a = matrix(R, [1 2; 3 4])
   @test reverse_rows(a) == matrix(R, [3 4; 1 2])
   reverse_rows!(a)
   @test a == matrix(R, [3 4; 1 2])

   a = matrix(R, [1 2; 3 4])
   @test reverse_cols(a) == matrix(R, [2 1; 4 3])
   reverse_cols!(a)
   @test a == matrix(R, [2 1; 4 3])

   a = matrix(R, [1 2 3; 3 4 5; 5 6 7])

   @test reverse_rows(a) == matrix(R, [5 6 7; 3 4 5; 1 2 3])
   reverse_rows!(a)
   @test a == matrix(R, [5 6 7; 3 4 5; 1 2 3])

   a = matrix(R, [1 2 3; 3 4 5; 5 6 7])
   @test reverse_cols(a) == matrix(R, [3 2 1; 5 4 3; 7 6 5])
   reverse_cols!(a)
   @test a == matrix(R, [3 2 1; 5 4 3; 7 6 5])
end

@testset "Generic.Mat.gen_mat_elem_op..." begin
   R, x = PolynomialRing(ZZ, "x")
   for i in 1:10
      r = rand(1:50)
      c = rand(1:50)
      S = MatrixSpace(R, r, c)
      M = rand(S, 0:3, -100:100)
      c1, c2 = rand(1:c), rand(1:c)
      s = rand(-100:100)
      r1 = rand(1:r)
      r2 = rand(r1:r)

      # add column

      N = add_column(M, s, c1, c2, r1:r2)
      for ci in 1:c
         if ci == c2
            @test all(N[k, c2] == M[k, c2] + s * M[k, c1] for k in r1:r2)
            @test all(N[k, c2] == M[k, c2] for k in 1:r if !(k in r1:r2))
         else
            @test all(N[k, ci] == M[k, ci] for k in 1:r)
         end
      end

      MM = deepcopy(M)
      add_column!(MM, s, c1, c2, r1:r2)
      for ci in 1:c
         if ci == c2
            @test all(MM[k, c2] == M[k, c2] + s * M[k, c1] for k in r1:r2)
            @test all(MM[k, c2] == M[k, c2] for k in 1:r if !(k in r1:r2))
         else
            @test all(MM[k, ci] == M[k, ci] for k in 1:r)
         end
      end

      if c1 != c2
         @test add_column(add_column(M, s, c1, c2), -s, c1, c2) == M
      end

      # multiply column

      N = multiply_column(M, s, c1, r1:r2)
      for ci in 1:c
         if ci == c1
            @test all(N[k, c1] == s * M[k, c1] for k in r1:r2)
            @test all(N[k, c1] == M[k, c1] for k in 1:r if !(k in r1:r2))
         else
            @test all(N[k, ci] == M[k, ci] for k in 1:r)
         end
      end

      MM = deepcopy(M)
      multiply_column!(MM, s, c1, r1:r2)
      for ci in 1:c
         if ci == c1
            @test all(MM[k, c1] == s * M[k, c1] for k in r1:r2)
            @test all(MM[k, c1] == M[k, c1] for k in 1:r if !(k in r1:r2))
         else
            @test all(MM[k, ci] == M[k, ci] for k in 1:r)
         end
      end

      @test multiply_column(multiply_column(M, -one(R), c1), -one(R), c1) == M

      # add row

      r1, r2 = rand(1:r), rand(1:r)
      s = rand(-100:100)
      c1 = rand(1:c)
      c2 = rand(c1:c)

      N = add_row(M, s, r1, r2, c1:c2)
      for ci in 1:r
         if ci == r2
            @test all(N[r2, k] == M[r2, k] + s * M[r1, k] for k in c1:c2)
            @test all(N[r2, k] == M[r2, k] for k in 1:c if !(k in c1:c2))
         else
            @test all(N[ci, k] == M[ci, k] for k in 1:c)
         end
      end

      MM = deepcopy(M)
      add_row!(MM, s, r1, r2, c1:c2)
      for ci in 1:r
         if ci == r2
            @test all(MM[r2, k] == M[r2, k] + s * M[r1, k] for k in c1:c2)
            @test all(MM[r2, k] == M[r2, k] for k in 1:c if !(k in c1:c2))
         else
            @test all(MM[ci, k] == M[ci, k] for k in 1:c)
         end
      end

      if r1 != r2
         @test add_row(add_row(M, s, r1, r2), -s, r1, r2) == M
      end

      # multiply row

      N = multiply_row(M, s, r1, c1:c2)
      for ci in 1:r
         if ci == r1
            @test all(N[r1, k] == s * M[r1, k] for k in c1:c2)
            @test all(N[r1, k] == M[r1, k] for k in 1:c if !(k in c1:c2))
         else
            @test all(N[ci, k] == M[ci, k] for k in 1:c)
         end
      end

      MM = deepcopy(M)
      multiply_row!(MM, s, r1, c1:c2)
      for ci in 1:r
         if ci == r1
            @test all(MM[r1, k] == s * M[r1, k] for k in c1:c2)
            @test all(MM[r1, k] == M[r1, k] for k in 1:c if !(k in c1:c2))
         else
            @test all(MM[ci, k] == M[ci, k] for k in 1:c)
         end
      end

      @test multiply_row(multiply_row(M, -one(R), r1), -one(R), r1) == M
   end
end

@testset "Generic.Mat.concat..." begin
   R, x = PolynomialRing(ZZ, "x")

   for i = 1:10
      r = rand(0:10)
      c1 = rand(0:10)
      c2 = rand(0:10)

      S1 = MatrixSpace(R, r, c1)
      S2 = MatrixSpace(R, r, c2)

      M1 = rand(S1, 0:3, -100:100)
      M2 = rand(S2, 0:3, -100:100)

      @test vcat(transpose(M1), transpose(M2)) == transpose(hcat(M1, M2))
   end

   A = matrix(R, 2, 2, [1, 2, 3, 4])
   B = matrix(R, 4, 2, [1, 2, 3, 4, 0, 1, 0, 1])
   C = matrix(R, 4, 1, [0, 1, 0, 2])
   D = matrix(R, 2, 3, [1, 2, 3, 4, 5, 6])

   @test hcat(B, C) == matrix(R, [1 2 0;
                                  3 4 1;
                                  0 1 0;
                                  0 1 2;])
   @test hcat(B, C) == [B C]
   @test hcat(B, C, C, B) == reduce(hcat, [B, C, C, B])
   let BC = hcat([B, C])
      @test size(BC) == (2, 1)
      @test BC[1] == B
      @test BC[2] == C
   end

   @test vcat(A, B) == matrix(R, [1 2;
                                  3 4;
                                  1 2;
                                  3 4;
                                  0 1;
                                  0 1;])

   @test vcat(A, B) == [A; B]
   @test vcat(A, B, B, A) == reduce(vcat, [A, B, B, A])
   let AB = vcat([A, B])
      @test size(AB) == (2,)
      @test AB[1] == A
      @test AB[2] == B
   end

   @test [A D; B B C] == matrix(R, [1 2 1 2 3;
                                    3 4 4 5 6;
                                    1 2 1 2 0;
                                    3 4 3 4 1;
                                    0 1 0 1 0;
                                    0 1 0 1 2;])
end

@testset "Generic.Mat.hnf_minors..." begin
   R, x = PolynomialRing(QQ, "x")

   M = MatrixSpace(R, 4, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))

   H = hnf_minors(A)
   @test ishnf(H)

   H, U = hnf_minors_with_transform(A)
   @test ishnf(H)
   @test isunit(det(U))
   @test U*A == H

   # Fake up finite field of char 7, degree 2
   R, x = PolynomialRing(GF(7), "x")
   F = ResidueField(R, x^2 + 6x + 3)
   a = F(x)

   S, y = PolynomialRing(F, "y")

   N = MatrixSpace(S, 4, 4)

   B = N(map(S, Any[1 0 a 0; a*y^3 0 3*a^2 0; y^4+a 0 y^2+y 5; y 1 y 2]))

   H = hnf_minors(B)
   @test ishnf(H)

   H, U = hnf_minors_with_transform(B)
   @test ishnf(H)
   @test isunit(det(U))
   @test U*B == H
end

@testset "Generic.Mat.hnf_kb..." begin
   M = matrix(ZZ, BigInt[4 6 2; 0 0 10; 0 5 3])

   H, U = AbstractAlgebra.hnf_kb_with_transform(M)

   @test H == matrix(ZZ, BigInt[4 1 9; 0 5 3; 0 0 10])
   @test isunit(det(U))
   @test U*M == H

   R, x = PolynomialRing(QQ, "x")

   M = MatrixSpace(R, 4, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))

   H = AbstractAlgebra.hnf_kb(A)
   @test ishnf(H)

   H, U = AbstractAlgebra.hnf_kb_with_transform(A)
   @test ishnf(H)
   @test isunit(det(U))
   @test U*A == H

   # Fake up finite field of char 7, degree 2
   R, x = PolynomialRing(GF(7), "x")
   F = ResidueField(R, x^2 + 6x + 3)
   a = F(x)

   S, y = PolynomialRing(F, "y")

   N = MatrixSpace(S, 3, 4)

   B = N(map(S, Any[1 0 a 0; a*y^3 0 3*a^2 0; y^4+a 0 y^2+y 5]))

   H = AbstractAlgebra.hnf_kb(B)
   @test ishnf(H)

   H, U = AbstractAlgebra.hnf_kb_with_transform(B)
   @test ishnf(H)
   @test isunit(det(U))
   @test U*B == H
end

@testset "Generic.Mat.hnf_cohen..." begin
   R, x = PolynomialRing(QQ, "x")

   M = MatrixSpace(R, 4, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))

   H = AbstractAlgebra.hnf_cohen(A)
   @test ishnf(H)

   H, U = AbstractAlgebra.hnf_cohen_with_transform(A)
   @test ishnf(H)
   @test isunit(det(U))
   @test U*A == H

   # Fake up finite field of char 7, degree 2
   R, x = PolynomialRing(GF(7), "x")
   F = ResidueField(R, x^2 + 6x + 3)
   a = F(x)

   S, y = PolynomialRing(F, "y")

   N = MatrixSpace(S, 3, 4)

   B = N(map(S, Any[1 0 a 0; a*y^3 0 3*a^2 0; y^4+a 0 y^2+y 5]))

   H = AbstractAlgebra.hnf_cohen(B)
   @test ishnf(H)

   H, U = AbstractAlgebra.hnf_cohen_with_transform(B)
   @test ishnf(H)
   @test isunit(det(U))
   @test U*B == H
end

@testset "Generic.Mat.hnf..." begin
   R, x = PolynomialRing(QQ, "x")

   M = MatrixSpace(R, 4, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))

   H = hnf(A)
   @test ishnf(H)

   H, U = hnf_with_transform(A)
   @test ishnf(H)
   @test isunit(det(U))
   @test U*A == H

   # Fake up finite field of char 7, degree 2
   R, x = PolynomialRing(GF(7), "x")
   F = ResidueField(R, x^2 + 6x + 3)
   a = F(x)

   S, y = PolynomialRing(F, "y")

   N = MatrixSpace(S, 3, 4)

   B = N(map(S, Any[1 0 a 0; a*y^3 0 3*a^2 0; y^4+a 0 y^2+y 5]))

   H = hnf(B)
   @test ishnf(H)

   H, U = hnf_with_transform(B)
   @test ishnf(H)
   @test isunit(det(U))
   @test U*B == H
end

@testset "Generic.Mat.snf_kb..." begin
   R, x = PolynomialRing(QQ, "x")

   M = MatrixSpace(R, 4, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))

   T = AbstractAlgebra.snf_kb(A)
   @test is_snf(T)

   T, U, K = AbstractAlgebra.snf_kb_with_transform(A)
   @test is_snf(T)
   @test isunit(det(U))
   @test isunit(det(K))
   @test U*A*K == T

   # Fake up finite field of char 7, degree 2
   R, x = PolynomialRing(GF(7), "x")
   F = ResidueField(R, x^2 + 6x + 3)
   a = F(x)

   S, y = PolynomialRing(F, "y")

   N = MatrixSpace(S, 3, 4)

   B = N(map(S, Any[1 0 a 0; a*y^3 0 3*a^2 0; y^4+a 0 y^2+y 5]))

   T = AbstractAlgebra.snf_kb(B)
   @test is_snf(T)

   T, U, K = AbstractAlgebra.snf_kb_with_transform(B)
   @test is_snf(T)
   @test isunit(det(U))
   @test isunit(det(K))
   @test U*B*K == T
end

@testset "Generic.Mat.snf..." begin
   R, x = PolynomialRing(QQ, "x")

   M = MatrixSpace(R, 4, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))

   T = snf(A)
   @test is_snf(T)

   T, U, K = snf_with_transform(A)
   @test is_snf(T)
   @test isunit(det(U))
   @test isunit(det(K))
   @test U*A*K == T

   # Fake up finite field of char 7, degree 2
   R, x = PolynomialRing(GF(7), "x")
   F = ResidueField(R, x^2 + 6x + 3)
   a = F(x)

   S, y = PolynomialRing(F, "y")

   N = MatrixSpace(S, 3, 4)

   B = N(map(S, Any[1 0 a 0; a*y^3 0 3*a^2 0; y^4+a 0 y^2+y 5]))

   T = snf(B)
   @test is_snf(T)

   T, U, K = snf_with_transform(B)
   @test is_snf(T)
   @test isunit(det(U))
   @test isunit(det(K))
   @test U*B*K == T
end

@testset "Generic.Mat.weak_popov..." begin
   R, x = PolynomialRing(QQ, "x")

   A = matrix(R, map(R, Any[1 2 3 x; x 2*x 3*x x^2; x x^2+1 x^3+x^2 x^4+x^2+1]))
   r = rank(A)

   P = weak_popov(A)
   @test is_weak_popov(P, r)

   P, U = weak_popov_with_transform(A)
   @test is_weak_popov(P, r)
   @test U*A == P
   @test isunit(det(U))

   F = GF(7)

   S, y = PolynomialRing(F, "y")

   B = matrix(S, map(S, Any[ 4*y^2+3*y+5 4*y^2+3*y+4 6*y^2+1; 3*y+6 3*y+5 y+3; 6*y^2+4*y+2 6*y^2 2*y^2+y]))
   s = rank(B)

   P = weak_popov(B)
   @test is_weak_popov(P, s)

   P, U = weak_popov_with_transform(B)
   @test is_weak_popov(P, s)
   @test U*B == P
   @test isunit(det(U))

   # some random tests

   for i in 1:3
      M = MatrixSpace(PolynomialRing(QQ, "x")[1], rand(1:5), rand(1:5))
      A = rand(M, 0:5, -5:5)
      r = rank(A)
      P = weak_popov(A)
      @test is_weak_popov(P, r)

      P, U = weak_popov_with_transform(A)
      @test is_weak_popov(P, r)
      @test U*A == P
      @test isunit(det(U))
   end

   R = GF(randprime(100))

   M = MatrixSpace(PolynomialRing(R, "x")[1], rand(1:5), rand(1:5))

   for i in 1:2
      A = rand(M, 1:5)
      r = rank(A)
      P = weak_popov(A)
      @test is_weak_popov(P, r)

      P, U = weak_popov_with_transform(A)
      @test is_weak_popov(P, r)
      @test U*A == P
      @test isunit(det(U))
   end

   R = ResidueField(ZZ, randprime(100))

   M = MatrixSpace(PolynomialRing(R, "x")[1], rand(1:5), rand(1:5))

   for i in 1:2
      A = rand(M, 1:5, 0:100)
      r = rank(A)
      P = weak_popov(A)
      @test is_weak_popov(P, r)

      P, U = weak_popov_with_transform(A)
      @test is_weak_popov(P, r)
      @test U*A == P
      @test isunit(det(U))
   end
end

@testset "Generic.Mat.views..." begin
   M = matrix(ZZ, 3, 3, BigInt[1, 2, 3, 2, 3, 4, 3, 4, 5])
   M2 = deepcopy(M)

   N1 = @view M[:,1:2]
   N2 = @view M[1:2, :]
   N3 = @view M[:,:]

   @test isa(N1, Generic.MatSpaceView)
   @test isa(N2, Generic.MatSpaceView)
   @test isa(N3, Generic.MatSpaceView)

   @test N2*N1 == matrix(ZZ, 2, 2, BigInt[14, 20, 20, 29])
   @test N3*N3 == M*M

   @test fflu(N3) == fflu(M) # tests that deepcopy is correct
   @test M2 == M
end

@testset "Generic.Mat.change_base_ring..." begin
   for (P, Q, T) in ((MatrixSpace(ZZ, 2, 3), MatrixSpace(ZZ, 3, 2), MatElem),
                     (MatrixAlgebra(ZZ, 3), MatrixAlgebra(ZZ, 3), MatAlgElem))
      M = rand(P, -10:10)
      N = rand(Q, -10:10)
      for R in [QQ, ZZ, GF(2), GF(5)]
         MQ = change_base_ring(R, M)
         @test MQ isa T
         @test base_ring(MQ) == R
         NQ = change_base_ring(R, N)
         @test NQ isa T
         @test base_ring(NQ) == R
         MNQ = change_base_ring(R, M * N)
         @test MNQ isa T
         @test base_ring(MNQ) == R
         @test MQ * NQ == MNQ
      end
   end
end

@testset "Generic.Mat.map..." begin
   u, v = rand(0:9, 2)
   for (mat, algebra) = ((rand(1:9, u, v), false),
                         (rand(1:9, u, u), true))
      for R = [QQ, ZZ, GF(2), GF(7), PolynomialRing(GF(5), 'x')[1]]
         M = algebra ? MatrixAlgebra(R, u) : MatrixSpace(R, u, v)
         m0 = M(mat)
         for f0 = (x -> x + 1, x -> x*2, x -> one(R), x -> zero(R))
            for f = (f0, map_from_func(f0, R, R))
               m = deepcopy(m0)
               n0 = similar(m)
               n = map_entries!(f, n0, m)
               @test n === n0 # map! must return its argument
               if !isempty(mat)
                  # when empty, it may happen that the result of map below has Any
                  # as eltype, and calling M on it fails, cf. issue #423
                  @test n == M(map(f isa Function ? f : f.image_fn, mat))
               end
            end
         end
      end
      m0 = algebra ? MatrixAlgebra(ZZ, u)(mat) : MatrixSpace(ZZ, u, v)(mat)
      m = deepcopy(m0)
      for S = [QQ, ZZ, GF(2), GF(7), PolynomialRing(GF(5), 'x')[1]]
         for f0 = (x -> S(x), x -> S(x + 1))
            for f = (f0, map_from_func(f0, ZZ, S))
               n = map_entries(f, m)
               @test n !== m
               @test m == m0 # map's input must not be mutated
               M = algebra ? MatrixAlgebra(S, u) : MatrixSpace(S, u, v)
               if !isempty(mat)
                  @test n == M(map(f isa Function ? f : f.image_fn, mat))
               end
               @test n isa (algebra ? MatAlgElem : MatElem)
            end
         end
      end
   end
end

@testset "Generic.Mat.similar/zero..." begin
   for sim_zero in (similar, zero)
      test_zero = sim_zero === zero
      for R = (ZZ, GF(11))
         M = MatrixSpace(R, rand(0:9), rand(0:9))
         m = R == ZZ ? rand(M, -10:10) : rand(M)
         n = sim_zero(m)
         @test !test_zero || iszero(n)
         @test parent(n) == M
         @test size(n) == (nrows(M), ncols(M))
         r, c = rand(0:9, 2)
         n = sim_zero(m, r, c)
         @test !test_zero || iszero(n)
         @test parent(n) == MatrixSpace(R, r, c)
         @test size(n) == (r, c)
         for S = [QQ, ZZ, GF(2), GF(5)]
            n = sim_zero(m, S)
            @test !test_zero || iszero(n)
            @test parent(n) == MatrixSpace(S, size(n)...)
            @test size(n) == (nrows(M), ncols(M))
            r, c = rand(0:9, 2)
            n = sim_zero(m, S, r, c)
            @test !test_zero || iszero(n)
            @test parent(n) == MatrixSpace(S, r, c)
            @test size(n) == (r, c)
         end
      end
   end
end

@testset "Generic.Mat.printing..." begin
   @test string(matrix(ZZ, [3 1 2; 2 0 1])) == "[3  1  2]\n[2  0  1]"
   @test string(matrix(ZZ, [3 1 2; 2 0 1])) == "[3  1  2]\n[2  0  1]"
   @test string(matrix(ZZ, 2, 0, [])) == "2 by 0 matrix"
   @test string(matrix(ZZ, 0, 3, [])) == "0 by 3 matrix"
   S = MatrixAlgebra(QQ, 3)
   @test string(S([1 2 3; 4 5 6; 7 8 9])) ==
      "[1//1  2//1  3//1]\n[4//1  5//1  6//1]\n[7//1  8//1  9//1]"
   @test string(MatrixAlgebra(QQ, 0)()) == "0 by 0 matrix"
   @test string(similar(matrix(ZZ, [3 1 2; 2 0 1]))) ==
      "[#undef  #undef  #undef]\n[#undef  #undef  #undef]"
end
