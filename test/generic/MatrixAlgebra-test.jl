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

@testset "Generic.MatAlg.constructors" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

   @test S === MatrixAlgebra(R, 3)

   @test elem_type(S) == Generic.MatAlgElem{elem_type(R)}
   @test elem_type(Generic.MatAlgebra{elem_type(R)}) == Generic.MatAlgElem{elem_type(R)}
   @test parent_type(Generic.MatAlgElem{elem_type(R)}) == Generic.MatAlgebra{elem_type(R)}

   @test typeof(S) <: Generic.MatAlgebra

   f = S(t^2 + 1)

   @test isa(f, MatAlgElem)

   g = S(2)

   @test isa(g, MatAlgElem)

   h = S(BigInt(23))

   @test isa(h, MatAlgElem)

   k = S([t t + 2 t^2 + 3t + 1; 2t R(2) t + 1; t^2 + 2 t + 1 R(0)])

   @test isa(k, MatAlgElem)

   l = S(k)

   @test isa(l, MatAlgElem)

   m = S()

   @test isa(m, MatAlgElem)

   n = identity_matrix(m)

   @test isa(n, MatAlgElem)

   @test nrows(n) == degree(S)

   @test_throws ErrorConstrDimMismatch S([t t^2 ; t^3 t^4])
   @test_throws ErrorConstrDimMismatch S([t t^2 t^3 ; t^4 t^5 t^6 ; t^7 t^8 t^9 ; t t^2 t^3])
   @test_throws ErrorConstrDimMismatch S([t, t^2])
   @test_throws ErrorConstrDimMismatch S([t, t^2, t^3, t^4, t^5, t^6, t^7, t^8, t^9, t^10])

   # Test constructors over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   
   S = MatrixAlgebra(R, 2)

   @test isa(S, MatAlgebra)

   @test base_ring(S) == R

   @test elem_type(S) == Generic.MatAlgElem{elem_type(R)}
   @test elem_type(Generic.MatAlgebra{elem_type(R)}) == Generic.MatAlgElem{elem_type(R)}
   @test parent_type(Generic.MatAlgElem{elem_type(R)}) == Generic.MatAlgebra{elem_type(R)}

   @test is_exact_type(elem_type(S)) == true
   @test is_domain_type(elem_type(S)) == false

   @test isa(S(), MatAlgElem)
   @test isa(S(ZZ(1)), MatAlgElem)
   @test isa(S(one(R)), MatAlgElem)
   @test isa(S([1 2; 3 4]), MatAlgElem)
   @test isa(S([1, 2, 3, 4]), MatAlgElem)

   @test parent(S()) == S
end

@testset "Generic.MatAlg.printing" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)
   @test !occursin("\n", sprint(show, S))
end

@testset "Generic.MatAlg.manipulation" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])

   @test nrows(S) == 3
   @test ncols(S) == 3

   @test is_exact_type(typeof(A))

   @test iszero(zero(S))
   @test isone(one(S))

   @test zero(A) == zero(S)
   @test one(A) == one(S)

   @test characteristic(S) == 0

   U = MatrixAlgebra(QQ, 3)
   C = U([1 2 3; 5 6 7; 9 8 5])

   @test !is_unit(A)
   @test is_unit(C)
   @test !is_zero_divisor(A)
   @test !is_zero_divisor(C)

   B[1, 1] = R(3)
   @test B[1, 1] == R(3)

   B[1, 1] = 4
   @test B[1, 1] == R(4)

   B[1, 1] = BigInt(5)
   @test B[1, 1] == R(5)

   @test nrows(B) == 3
   @test ncols(B) == 3

   @test deepcopy(A) == A

   C = S([t + 1 R(0) R(1); t^2 R(0) t; R(0) R(0) R(0)])

   @test is_zero_row(C, 3)
   @test !is_zero_row(C, 1)
   @test is_zero_column(C, 2)
   @test !is_zero_column(C, 1)

   S = MatrixAlgebra(QQ, 3)
   A = S([1 2 3; 4 5 6; 7 8 9])

   @test nrows(S) == ncols(S) == 3
   @test degree(S) == 3

   @test nrows(A) == ncols(A) == 3
   @test degree(A) == 3

   # Tests over residue ring
   S = MatrixAlgebra(residue_ring(ZZ, 6)[1], 2)
   A = S([1 2; 3 4])
   B = S([0 0; 3 3])
   @test is_zero(A*B)
   @test !is_unit(A)
   @test !is_unit(B)
   @test is_zero_divisor(A)
   @test is_zero_divisor(B)
   @test !is_zero_divisor(one(S))

   # Tests over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   
   S = MatrixAlgebra(R, 2)

   M = rand(S, -10:10)

   @test isa(hash(M), UInt)
   @test nrows(M) == 2
   @test ncols(M) == 2
   @test length(M) == 4
   @test isempty(M) == false
   @test isassigned(M, 1, 1) == true

   @test iszero(zero(M, 3, 3))
   @test iszero(zero(M, QQ, 3, 3))
   @test iszero(zero(M, QQ))
   
   zero!(M)
   @test iszero(M)

   @test isone(one(R))

   @test is_zero_row(M, 1)
   @test is_zero_column(M, 1)

   @test degree(M) == 2
end

@testset "Generic.MatAlg.size/axes" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   
   @test is_square(A)
   @test size(A) == (3, 3)
   @test size(A, 1) == 3
   @test size(A, 2) == 3
   @test_throws BoundsError size(A, 0)
   @test_throws BoundsError size(A, -rand(1:99))
   @test axes(A) == (1:3, 1:3)
   @test axes(A, 1) == 1:3
   @test axes(A, 2) == 1:3
   @test axes(A, rand(3:99)) == 1:1
   @test_throws BoundsError axes(A, 0)
   @test_throws BoundsError axes(A, -rand(1:99))

   # test over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   
   S = MatrixAlgebra(R, 2)

   M = rand(S, -10:10)

   @test firstindex(M, 1) == 1
   @test lastindex(M, 1) == 2
   @test size(M) == (2, 2)
   @test size(M, 1) == 2
   @test axes(M) == (1:2, 1:2)
   @test axes(M, 1) == 1:2
   @test is_square(M)
end

@testset "Generic.MatAlg.unary_ops" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([-t - 1 (-t) -R(1); -t^2 (-t) (-t); -R(-2) (-t - 2) (-t^2 - t - 1)])

   @test -A == B

   # Tests over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   
   S = MatrixAlgebra(R, 2)

   M = rand(S, -10:10)

   @test -(-M) == M
end

@testset "Generic.MatAlg.binary_ops" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])

   @test A + B == S([t+3 t+3 R(2); t^2 + t 2*t+1 2*t+2; R(-3) t^2 + t + 2 t^3 + 1*t^2 + t + 1])

   @test A - B == S([t-1 t-3 R(0); t^2 - t R(-1) R(-2); R(-1) (-t^2 + t + 2) (-t^3 + t^2 + t + 1)])

   @test A*B == S([t^2 + 2*t + 1 2*t^2 + 4*t + 3 t^3 + t^2 + 3*t + 1; 3*t^2 - t (t^3 + 4*t^2 + t) t^4 + 2*t^2 + 2*t; t-5 t^4 + t^3 + 2*t^2 + 3*t - 4 t^5 + 1*t^4 + t^3 + t^2 + 4*t + 2])

   # Tests over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   
   S = MatrixAlgebra(R, 2)

   M = rand(S, -10:10)
   N = rand(S, -10:10)
   P = rand(S, -10:10)

   @test M + N == N + M
   @test M - N == M + (-N)
   @test M*(N + P) == M*N + M*P
end

@testset "Generic.MatAlg.adhoc_binary" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

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

   # Tests over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   
   S = MatrixAlgebra(R, 2)
   
   M = rand(S, -10:10)
   N = rand(S, -10:10)
   
   t1 = rand(ZZ, -10:10)
   t2 = rand(R, -10:10)

   @test t1*(M + N) == t1*M + t1*N
   @test t1*(M - N) == t1*M - t1*N
   @test (M + N)*t1 == M*t1 + N*t1
   @test (M - N)*t1 == M*t1 - N*t1

   @test t2*(M + N) == t2*M + t2*N
   @test t2*(M - N) == t2*M - t2*N
   @test (M + N)*t2 == M*t2 + N*t2
   @test (M - N)*t2 == M*t2 - N*t2

   @test M + t1 == M - (-t1)
   @test M + t2 == M - (-t2)

   @test t1 + M == t1 - (-M)
   @test t2 + M == t2 - (-M)

   r1 = rand(R, -10:10)
   r2 = rand(R, -10:10)

   @test (M + N)*[r1, r2] == M*[r1, r2] + N*[r1, r2]
   @test [r1, r2]*(M + N) == [r1, r2]*M + [r1, r2]*N
end

@testset "Generic.MatAlg.promotion" begin
   m = [1 2; 3 4]
   F, = residue_field(ZZ, 3)
   R, t = polynomial_ring(F, "t")
   A = MatrixAlgebra(R, 2)(m)
   B = MatrixAlgebra(F, 2)(m)

   @test typeof(A * B) == typeof(A)
   @test typeof(B * A) == typeof(A)
   @test A * B == A^2
   @test B * A == A^2
   
   @test typeof(A + B) == typeof(A)
   @test typeof(B + A) == typeof(A)
   @test A + B == A + A
   @test B + A == A + A

   @test typeof(F(2) * A) == typeof(A)
   @test typeof(A * F(2)) == typeof(A)
   @test F(2) * A == R(2) * A
   @test A * F(2) == A * R(2)

   @test typeof(F(2) + A) == typeof(A)
   @test typeof(A + F(2)) == typeof(A)
   @test F(2) + A == R(2) + A
   @test A + F(2) == A + R(2)

   @test one(F) == R[1 0; 0 1]
   @test R[1 0; 0 1] == one(R)

   # vector * matrix
   m = [1 2; 3 4]
   F, = residue_field(ZZ, 3)
   R, t = polynomial_ring(F, "t")
   A = MatrixAlgebra(R, 2)(m)
   B = MatrixAlgebra(F, 2)(m)
   v = [one(F), 2*one(F)]
   vv = [one(R), 2*one(R)]
   @test (@inferred A * v) == A * vv
   @test (@inferred v * A) == vv * A

   @test (@inferred B * vv) == A * vv
   @test (@inferred vv * B) == vv * A

   @test_throws ErrorException A * Rational{BigInt}[1 ,2]
   @test_throws ErrorException Rational{BigInt}[1 ,2] * A
end

@testset "Generic.MatAlg.permutation" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   T = SymmetricGroup(3)
   P = T([2, 3, 1])

   @test A == inv(P)*(P*A)
end

@testset "Generic.MatAlg.comparison" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test A == B

   @test A != one(S)


   # Tests over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   
   S = MatrixAlgebra(R, 2)
   
   M = rand(S, -10:10)
   N = deepcopy(M)
   
   @test M == M
   @test M == N
   @test M == copy(M)
   @test isequal(M, M)
   @test isequal(M, N)
   @test isequal(M, copy(M))
end

@testset "Generic.MatAlg.adhoc_comparison" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

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

   # Tests over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   
   S = MatrixAlgebra(R, 2)
   
   @test S(5) == 5
   @test 5 == S(5)
   @test S(BigInt(5)) == 5
   @test 5 == S(BigInt(5))

   m = rand(R, -10:10)

   @test S(m) == m
   @test m == S(m)
end

@testset "Generic.MatAlg.powering" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test A^5 == A^2*A^3

   @test A^0 == one(S)

   # Tests over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   
   S = MatrixAlgebra(R, 2)
   
   M = rand(S, -10:10)

   @test M^0 == one(S)
   @test M^1 == M
   @test M^2 == M*M
   @test M^3 == M*M*M
end

@testset "Generic.MatAlg.exact_division" begin
   S = MatrixAlgebra(QQ, 3)

   M = rand(S, -20:20)
   N = randmat_with_rank(S, 3, -20:20)

   @test divexact_right(M*N, N) == M
   @test divexact_left(N*M, N) == M
end

@testset "Generic.MatAlg.adhoc_exact_division" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test divexact(5*A, 5) == A
   @test divexact(12*A, BigInt(12)) == A
   @test divexact(12*A, Rational{BigInt}(12)) == A
   @test divexact((1 + t)*A, 1 + t) == A

   # Tests over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   U, x = polynomial_ring(R, "x")

   S = MatrixAlgebra(R, 2)
   T = MatrixAlgebra(U, 2)

   for i = 1:50
       M = rand(S, -10:10)

       @test divexact(5*M, 5) == M

       c = rand(R, 1:10)  # create a regular matrix
       while rank(c) < nrows(c)
          c = rand(R, 1:10)
       end

       @test divexact_left(c*M, c) == M
       @test divexact_right(M*c, c) == M

       N = rand(T, 0:5, -10:10)
       d = rand(U, 0:5, -10:10)
       while iszero(d) || rank(leading_coefficient(d)) != 2
          d = rand(U, 0:5, -10:10)
       end

       @test divexact_left(d*N, d) == N
       @test divexact_right(N*d, d) == N
   end
end

@testset "Generic.MatAlg.transpose" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)
   arr = [t + 1 t R(1); t^2 t t; t+1 t^2 R(-1)]
   A = S(arr)
   B = S(permutedims(arr, [2, 1]))
   @test transpose(A) == B

   # Tests over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   
   S = MatrixAlgebra(R, 2)
   
   M = rand(S, -10:10)

   @test is_symmetric(M + transpose(M))
end

@testset "Generic.MatAlg.gram" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test gram(A) == S([2*t^2 + 2*t + 2 t^3 + 2*t^2 + t 2*t^2 + t - 1; t^3 + 2*t^2 + t t^4 + 2*t^2 t^3 + 3*t; 2*t^2 + t - 1 t^3 + 3*t t^4 + 2*t^3 + 4*t^2 + 6*t + 9])
end

@testset "Generic.MatAlg.tr" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test tr(A) == t^2 + 3t + 2
end

@testset "Generic.MatAlg.content" begin
   R, t = polynomial_ring(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test content((1 + t)*A) == 1 + t
end

@testset "Generic.MatAlg.lu" begin
   R, x = polynomial_ring(QQ, "x")
   K, f = residue_field(R, x^3 + 3x + 1)
   a = f(x)
   S = MatrixAlgebra(K, 3)

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

   R, z = polynomial_ring(ZZ, "z")
   F = fraction_field(R)
   S = MatrixAlgebra(F, 3)

   A = S([F(0), F(0), F(11), 78*z^3-102*z^2+48*z+12, F(92), -16*z^2+80*z-149, -377*z^3+493*z^2-232*z-58, F(-448), 80*z^2-385*z+719])

   r, P, L, U = lu(A)

   @test r == 3
   @test P*A == L*U
end

@testset "Generic.MatAlg.fflu" begin
   R, x = polynomial_ring(QQ, "x")
   K, f = residue_field(R, x^3 + 3x + 1)
   a = f(x)
   S = MatrixAlgebra(K, 3)

   A = S([a + 1 2a + 3 a^2 + 1; 2a^2 - 1 a - 1 2a; a^2 + 3a + 1 2a K(1)])

   r, d, P, L, U = fflu(A)

   D = S()
   D[1, 1] = inv(L[1, 1])
   D[2, 2] = inv(L[1, 1]*L[2, 2])
   D[3, 3] = inv(L[2, 2]*L[3, 3])

   @test r == 3
   @test P*A == L*D*U

   A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])

   r, d, P, L, U = fflu(A)

   D = S()
   D[1, 1] = inv(L[1, 1])
   D[2, 2] = inv(L[1, 1]*L[2, 2])
   D[3, 3] = inv(L[2, 2]*L[3, 3])

   @test r == 3
   @test P*A == L*D*U

   A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 - 2 a - 1 2a])

   r, d, P, L, U = fflu(A)

   D = S()
   D[1, 1] = inv(L[1, 1])
   D[2, 2] = inv(L[1, 1]*L[2, 2])
   D[3, 3] = inv(L[2, 2]*L[3, 3])

   @test r == 2
   @test P*A == L*D*U

   S = MatrixAlgebra(QQ, 3)
   A = S([0, 0, 1, 12, 1, 11, 1, 0, 1])

   r, d, P, L, U, = fflu(A)

   D = S()
   D[1, 1] = inv(L[1, 1])
   D[2, 2] = inv(L[1, 1]*L[2, 2])
   D[3, 3] = inv(L[2, 2]*L[3, 3])
   @test r == 3
   @test P*A == L*D*U
end

@testset "Generic.MatAlg.det" begin
   S, x = polynomial_ring(residue_ring(ZZ, 1009*2003)[1], "x")

   for dim = 0:5
      R = MatrixAlgebra(S, dim)

      M = rand(R, -1:5, -100:100)

      @test det(M) == AbstractAlgebra.det_clow(M)
   end

   S, z = polynomial_ring(ZZ, "z")

   for dim = 0:5
      R = MatrixAlgebra(S, dim)

      M = rand(R, -1:3, -20:20)

      @test det(M) == AbstractAlgebra.det_clow(M)
   end

   R, x = polynomial_ring(QQ, "x")
   K, f = residue_field(R, x^3 + 3x + 1)
   a = f(x)

   for dim = 0:7
      S = MatrixAlgebra(K, dim)

      M = rand(S, -100:100)

      @test det(M) == AbstractAlgebra.det_clow(M)
   end

   R, x = polynomial_ring(ZZ, "x")
   S, y = polynomial_ring(R, "y")

   for dim = 0:5
      T = MatrixAlgebra(S, dim)
      M = rand(T, -1:2, -1:2, -10:10)

      @test det(M) == AbstractAlgebra.det_clow(M)
   end
end

@testset "Generic.MatAlg.block_replacement" begin
   _test_block_replacement = function(a, b, r, c)
      rr = r isa Colon ? (1:nrows(a)) : r
      cc = c isa Colon ? (1:ncols(a)) : c
      if (b isa Vector)
         return all([a[i1, j1] == b[i + j - 1] for (i, i1) in enumerate(rr) for (j, j1) in enumerate(cc)])
      else
         return all([a[i1, j1] == b[i, j] for (i, i1) in enumerate(rr) for (j, j1) in enumerate(cc)])
      end
   end

   S = MatrixAlgebra(ZZ, 9)
   (r, c) = (rand(1:9), rand(1:9))
   T = matrix_space(ZZ, r, c)
   a = rand(S, -100:100)
   b = rand(T, -100:100)
   startr = rand(1:(9-r+1))
   endr = startr + r - 1
   startc = rand(1:(9-c+1))
   endc = startc + c - 1
   a[startr:endr, startc:endc] = b
   @test _test_block_replacement(a, b, startr:endr, startc:endc)

   for i in 1:10
      n = rand(1:9)
      m = n
      S = MatrixAlgebra(ZZ, n)
      a = rand(S, -100:100)
      (r, c) = (rand(1:n), rand(1:m))
      T = matrix_space(zz, r, c)
      b = rand(T, -2:2)
      startr = rand(1:(n-r+1))
      endr = startr + r - 1
      startc = rand(1:(m-c+1))
      endc = startc + c - 1

      rrs = Int[]
      for j in 1:rand(1:n)
         rr = rand(1:n)
         while (rr in rrs)
            rr = rand(1:n)
         end
         push!(rrs, rr)
      end

      ccs = Int[]
      for j in 1:rand(1:m)
         cc = rand(1:m)
         while (cc in ccs)
            cc = rand(1:m)
         end
         push!(ccs, cc)
      end

      for r in [rand(1:n), Colon(), rrs, startr:endr]
         for c in [rand(1:m), Colon(), ccs, startc:endc]
            if c isa Int && r isa Int
               continue
            end
            for R in [zz, ZZ]
               lr = r isa Colon ? nrows(a) : length(r)
               lc = c isa Colon ? ncols(a) : length(c)
               T = matrix_space(zz, lr, lc)
               b = rand(T, -2:2)
               aa = deepcopy(a)
               a[r, c] = b
               @test _test_block_replacement(a, b, r, c)
               bb = Matrix(b)
               aa[r, c] = bb
               @test _test_block_replacement(aa, bb, r, c)
            end
         end
      end
      for r in [rand(1:n)]
         for c in [rand(1:m), Colon(), ccs, startc:endc]
            if c isa Int && r isa Int
               continue
            end
            for R in [zz, ZZ]
               lr = r isa Colon ? nrows(a) : length(r)
               lc = c isa Colon ? ncols(a) : length(c)
               T = matrix_space(zz, lr, lc)
               _b = rand(T, -2:2)
               b = vec(Matrix(_b))
               a[r, c] = b
               @test _test_block_replacement(a, b, r, c)
            end
         end
      end
      for r in [rand(1:n), Colon(), rrs, startr:endr]
         for c in [rand(1:m)]
            if c isa Int && r isa Int
               continue
            end
            for R in [zz, ZZ]
               lr = r isa Colon ? nrows(a) : length(r)
               lc = c isa Colon ? ncols(a) : length(c)
               T = matrix_space(zz, lr, lc)
               _b = rand(T, -2:2)
               b = vec(Matrix(_b))
               a[r, c] = b
               @test _test_block_replacement(a, b, r, c)
            end
         end
      end

   end
end

@testset "Generic.MatAlg.rank" begin
   S, = residue_ring(ZZ, 20011*10007)
   R = MatrixAlgebra(S, 5)

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

   S, z = polynomial_ring(ZZ, "z")
   R = MatrixAlgebra(S, 4)

   M = R([S(-2) S(0) S(5) S(3); 5*z^2+5*z-5 S(0) S(-z^2+z) 5*z^2+5*z+1; 2*z-1 S(0) z^2+3*z+2 S(-4*z); 3*z-5 S(0) S(-5*z+5) S(1)])

   @test rank(M) == 3

   R = matrix_space(S, 5, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, 0:3, -20:20)

      @test rank(M) == i
   end

   R, x = polynomial_ring(QQ, "x")
   K, f = residue_field(R, x^3 + 3x + 1)
   a = f(x)

   S = MatrixAlgebra(K, 3)

   M = S([a a^2 + 2*a - 1 2*a^2 - 1*a; 2*a+2 2*a^2 + 2*a (-2*a^2 - 2*a); (-a) (-a^2) a^2])

   @test rank(M) == 2

   S = MatrixAlgebra(K, 5)

   for i = 0:5
      M = randmat_with_rank(S, i, -100:100)

      @test rank(M) == i
   end

   R, x = polynomial_ring(ZZ, "x")
   S, y = polynomial_ring(R, "y")
   T = MatrixAlgebra(S, 3)

   M = T([(2*x^2)*y^2+(-2*x^2-2*x)*y+(-x^2+2*x) S(0) (-x^2-2)*y^2+(x^2+2*x+2)*y+(2*x^2-x-1);
    (-x)*y^2+(-x^2+x-1)*y+(x^2-2*x+2) S(0) (2*x^2+x-1)*y^2+(-2*x^2-2*x-2)*y+(x^2-x);
    (-x+2)*y^2+(x^2+x+1)*y+(-x^2+x-1) S(0) (-x^2-x+2)*y^2+(-x-1)*y+(-x-1)])

   @test rank(M) == 2

   T = MatrixAlgebra(S, 5)

   for i = 0:5
      M = randmat_with_rank(T, i, 0:2, 0:2, -20:20)

      @test rank(M) == i
   end
end

@testset "Generic.MatAlg.solve_lu" begin
   S = QQ

   for dim = 0:5
      R = MatrixAlgebra(S, dim)
      U = MatrixAlgebra(S, dim)

      M = randmat_with_rank(R, dim, -100:100)
      b = rand(U, -100:100)

      flag, x = Generic.can_solve_with_solution_lu(M, b)

      @test flag && M*x == b
   end

   S, y = polynomial_ring(ZZ, "y")
   K = fraction_field(S)

   for dim = 0:5
      R = MatrixAlgebra(S, dim)
      U = MatrixAlgebra(S, dim)
      T = MatrixAlgebra(K, dim)

      M = randmat_with_rank(R, dim, 0:5, -100:100)
      b = rand(U, 0:5, -100:100);

      MK = T(elem_type(K)[ K(M[i, j]) for i in 1:nrows(M), j in 1:ncols(M) ])
      bK = T(elem_type(K)[ K(b[i, j]) for i in 1:nrows(b), j in 1:ncols(b) ])

      flag, x = Generic.can_solve_with_solution_lu(MK, bK)

      @test flag && MK*x == bK
   end
end

@testset "Generic.MatAlg.rref" begin
   S, = residue_ring(ZZ, 20011*10007)
   R = MatrixAlgebra(S, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, -100:100)

      do_test = false
      r = 0
      A = M
      try
         r, A, d = rref_rational(M)
         do_test = true
      catch e
         if !(e isa ErrorException)
            rethrow(e)
         end
      end

      if do_test
         @test r == i
         @test is_rref(A)
      end
   end

   S, z = polynomial_ring(ZZ, "z")
   R = MatrixAlgebra(S, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, 0:3, -20:20)

      r, A, d = rref_rational(M)

      @test r == i
      @test is_rref(A)
   end

   R, x = polynomial_ring(QQ, "x")
   K, f = residue_field(R, x^3 + 3x + 1)
   a = f(x)
   S = MatrixAlgebra(K, 5)

   for i = 0:5
      M = randmat_with_rank(S, i, -100:100)

      r, A = rref(M)

      @test r == i
      @test is_rref(A)
   end

   R, x = polynomial_ring(ZZ, "x")
   S, y = polynomial_ring(R, "y")
   T = MatrixAlgebra(S, 5)

   for i = 0:5
      M = randmat_with_rank(T, i, 0:2, 0:2, -20:20)

      r, A, d = rref_rational(M)

      @test r == i
      @test is_rref(A)
   end
end

@testset "Generic.MatAlg.inversion" begin
   indexing(n) = [(i,j) for i in 1:n for j in 1:n if i !=j ]
   E(R,i,j, val=1) = (M=one(R); M[i,j] = val; return M)
   E(R::MatAlgebra; vals=[1,-1]) = [E(R, i,j,val) for (i,j) in indexing(R.n) for val in vals]
   random_product(S::Vector{<:NCRingElem}, len=10) = prod(i->S[i], rand(1:length(S), len))

   @testset "Matrix Algebra over ZZ" begin
   S = ZZ

   for dim = 2:5
      R = MatrixAlgebra(S, dim)
      M = R(1)
      i = rand(1:dim-1)
      j = rand(i+1:dim)
      M[i,j] = 1 # E_{i,j} elementary matrix

      @test inv(M) isa MatAlgElem
      N = inv(M)
      @test N[i,j] == -1
      @test M*N == N*M == R(1)

      M[j,i] = -1
      @test_throws DomainError inv(M) # we would need to invert 2
      M[i,i] = 0
      @test inv(M) isa MatAlgElem
      NN = inv(M)
      @test NN[i,j] == -1
      @test NN[j,i] == 1

      @test M*NN == NN*M == one(R)

      gens = E(R)
      random_matrices = [random_product(gens) for _ in 1:10]

      @test all(isone(m*inv(m)) for m in random_matrices)
   end
   end

   @testset "Matrix Algebra over residue_ring" begin
   S, = residue_ring(ZZ, 20011*10007)

   for dim = 1:5
      R = MatrixAlgebra(S, dim)

      M = randmat_with_rank(R, dim, -100:100)

      do_test = false
      X = M

      try
         X = inv(M)
         do_test = true
      catch e
         e isa DomainError ? nothing : rethrow(e)
      end

      do_test && @test isone(M*X)
   end
   end

   @testset "Matrix Algebra over ZZ[x]" begin
   S, x = polynomial_ring(ZZ, "x")

   for dim = 2:5
      R = MatrixAlgebra(S, dim)
      M = one(R)
      i = rand(1:dim-1)
      j = rand(i+1:dim)
      M[i,j] = 1
      @test inv(M) isa MatAlgElem
      X = inv(M)
      @test isone(M*X)

      M[i,j] = x
      @test inv(M) isa MatAlgElem
      X = inv(M)
      @test isone(M*X)

      M[j,i] = 1
      @test_throws DomainError inv(M)

      gens = E(R, vals=[1, x, -x^2])
      random_matrices = [random_product(gens) for _ in 1:10]

      @test all(isone(m*inv(m)) for m in random_matrices)
   end
   end

   @testset "Matrix Algebra over number_field over QQ" begin
   R, x = polynomial_ring(QQ, "x")
   S, f = residue_field(R, x^3 + 3x + 1)
   a = f(x)

   for dim = 1:5
      R = MatrixAlgebra(S, dim)

      M = randmat_with_rank(R, dim, -100:100)

      X = inv(M)

      @test isone(M*X)
   end
   end

   @testset "Matrix Algebra over (ZZ[x])[y]" begin
   R, x = polynomial_ring(ZZ, "x")
   S, y = polynomial_ring(R, "y")

   for dim = 2:5
      T = MatrixAlgebra(S, dim)
      M = one(T)
      i = rand(1:dim-1)
      j = rand(i+1:dim)
      M[i,j] = R(1)

      @test inv(M) isa MatAlgElem
      @test isone(M*inv(M))

      M[i,j] = x
      @test isone(M*inv(M))

      M[j,i] = 1
      @test_throws DomainError inv(M)
      M[j,i] = 0

      M[i,j] = y
      @test isone(M*inv(M))
      M[j,i] = x
      @test_throws DomainError inv(M)

      gens = E(T, vals=[1,x,y])
      random_matrices = [random_product(gens) for _ in 1:10]

      @test all(isone(m*inv(m)) for m in random_matrices)
   end
   end
end # of @testset "Generic.MatAlg.inversion"

@testset "Generic.MatAlg.hessenberg" begin
   R, = residue_ring(ZZ, 18446744073709551629)

   for dim = 0:5
      S = MatrixAlgebra(R, dim)

      for i = 1:10
         M = rand(S, -5:5)

         A = hessenberg(M)

         @test is_hessenberg(A)
      end
   end
end

@testset "Generic.MatAlg.charpoly" begin
   R, = residue_ring(ZZ, 18446744073709551629)

   for dim = 0:5
      S = matrix_space(R, dim, dim)
      U, x = polynomial_ring(R, "x")

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

   R, x = polynomial_ring(ZZ, "x")
   U, z = polynomial_ring(R, "z")
   T = matrix_space(R, 6, 6)

   M = T()
   for i = 1:3
      for j = 1:3
         M[i, j] = rand(R, -1:2, -10:10)
         M[i + 3, j + 3] = deepcopy(M[i, j])
      end
   end

   p1 = charpoly(U, M)

   for i = 1:10
      similarity!(M, rand(1:6), R(rand(R, 0:2, -3:3)))
   end

   p2 = charpoly(U, M)

   @test p1 == p2

   p3 = charpoly(M)

   @test p1 == evaluate(p3, gen(U))
end

@testset "Generic.MatAlg.minpoly" begin
   R = GF(103)
   T, y = polynomial_ring(R, "y")
   S = MatrixAlgebra(R, 3)

   M = S([92 97 8;
          0 5 13;
          0 16 2])

   @test minpoly(T, M) == y^2+96*y+8

   R = GF(3)
   T, y = polynomial_ring(R, "y")
   S = MatrixAlgebra(R, 4)

   M = S([1 2 0 2;
         1 2 1 0;
         1 2 2 1;
         2 1 2 0])

   @test minpoly(T, M) == y^2 + 2y

   R = GF(13)
   T, y = polynomial_ring(R, "y")
   S = MatrixAlgebra(R, 3)

   M = S([7 6 1;
         7 7 5;
         8 12 5])

   @test minpoly(T, M) == y^2+10*y
   S = MatrixAlgebra(R, 4)

   M = S([4 0 9 5;
         1 0 1 9;
         0 0 7 6;
         0 0 3 10])

   @test minpoly(T, M) == y^2 + 9y
   S = MatrixAlgebra(R, 6)

   M = S([2 7 0 0 0 0;
         1 0 0 0 0 0;
         0 0 2 7 0 0;
         0 0 1 0 0 0;
         0 0 0 0 4 3;
         0 0 0 0 1 0])

   @test minpoly(T, M) == (y^2+9*y+10)*(y^2+11*y+6)
   S = MatrixAlgebra(R, 6)

   M = S([2 7 0 0 0 0;
         1 0 1 0 0 0;
         0 0 2 7 0 0;
         0 0 1 0 0 0;
         0 0 0 0 4 3;
         0 0 0 0 1 0])

   @test minpoly(T, M) == (y^2+9*y+10)*(y^2+11*y+6)^2

   S = MatrixAlgebra(R, 1)
   M = S()

   @test minpoly(T, M) == y

   S = MatrixAlgebra(R, 0)
   M = S()

   @test minpoly(T, M) == 1

   R, x = polynomial_ring(ZZ, "x")
   S, y = polynomial_ring(R, "y")
   U, z = polynomial_ring(S, "z")
   T = MatrixAlgebra(S, 6)

   M = T()
   for i = 1:3
      for j = 1:3
         M[i, j] = rand(S, -1:3, -1:3, -10:10)
         M[i + 3, j + 3] = deepcopy(M[i, j])
      end
   end

   f = minpoly(U, M)

   @test degree(f) <= 3

   R, x = polynomial_ring(ZZ, "x")
   U, z = polynomial_ring(R, "z")
   T = MatrixAlgebra(R, 6)

   M = T()
   for i = 1:3
      for j = 1:3
         M[i, j] = rand(R, -1:2, -10:10)
         M[i + 3, j + 3] = deepcopy(M[i, j])
      end
   end

   p1 = minpoly(U, M)

   for i = 1:10
      similarity!(M, rand(1:6), R(rand(R, -1:2, -3:3)))
   end

   p2 = minpoly(U, M)

   @test p1 == p2

   p3 = minpoly(M)

   @test p1 == evaluate(p3, gen(U))
end

@testset "Generic.MatAlg.row_swapping" begin
   R, x = polynomial_ring(ZZ, "x")
   M = MatrixAlgebra(R, 3)

   a = M(map(R, [1 2 3; 4 5 6; 7 8 9]))

   @test swap_rows(a, 1, 3) == M(map(R, [7 8 9; 4 5 6; 1 2 3]))

   swap_rows!(a, 2, 3)

   @test a == M(map(R, [1 2 3; 7 8 9; 4 5 6]))
end

if false # see bug 160
    @testset "Generic.MatAlg.hnf_minors" begin
        R, x = polynomial_ring(QQ, "x")

        M = MatrixAlgebra(R, 3)

        A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5]))

        H = hnf_minors(A)
        @test is_upper_triangular(H)

        H, U = hnf_minors_with_transform(A)
        @test is_upper_triangular(H)
        @test is_unit(det(U))
        @test U*A == H

        # Fake up finite field of char 7, degree 2
        R, x = polynomial_ring(GF(7), "x")
        F, = residue_field(R, x^2 + 6x + 3)
        a = F(x)

        S, y = polynomial_ring(F, "y")

        N = MatrixAlgebra(S, 4)

        B = N(map(S, Any[1 0 a 0; a*y^3 0 3*a^2 0; y^4+a 0 y^2+y 5; y 1 y 2]))

        H = hnf_minors(B)
        @test is_upper_triangular(H)

        H, U = hnf_minors_with_transform(B)
        @test is_upper_triangular(H)
        @test is_unit(det(U))
        @test U*B == H
    end
end

@testset "Generic.MatAlg.hnf_kb" begin
   R, x = polynomial_ring(QQ, "x")

   M = MatrixAlgebra(R, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5]))

   H = AbstractAlgebra.hnf_kb(A)
   @test is_upper_triangular(H)

   H, U = AbstractAlgebra.hnf_kb_with_transform(A)
   @test is_upper_triangular(H)
   @test is_unit(det(U))
   @test U*A == H

   # Fake up finite field of char 7, degree 2
   R, x = polynomial_ring(GF(7), "x")
   F, = residue_field(R, x^2 + 6x + 3)
   a = F(x)

   S, y = polynomial_ring(F, "y")

   N = MatrixAlgebra(S, 3)

   B = N(map(S, Any[1 0 a; a*y^3 0 3*a^2; y^4+a 0 y^2+y]))

   H = AbstractAlgebra.hnf_kb(B)
   @test is_upper_triangular(H)

   H, U = AbstractAlgebra.hnf_kb_with_transform(B)
   @test is_upper_triangular(H)
   @test is_unit(det(U))
   @test U*B == H
end

@testset "Generic.MatAlg.hnf_cohen" begin
   R, x = polynomial_ring(QQ, "x")

   M = MatrixAlgebra(R, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5]))

   H = AbstractAlgebra.hnf_cohen(A)
   @test is_upper_triangular(H)

   H, U = AbstractAlgebra.hnf_cohen_with_transform(A)
   @test is_upper_triangular(H)
   @test is_unit(det(U))
   @test U*A == H

   # Fake up finite field of char 7, degree 2
   R, x = polynomial_ring(GF(7), "x")
   F, = residue_field(R, x^2 + 6x + 3)
   a = F(x)

   S, y = polynomial_ring(F, "y")

   N = MatrixAlgebra(S, 3)

   B = N(map(S, Any[1 0 a; a*y^3 0 3*a^2; y^4+a 0 y^2+y]))

   H = AbstractAlgebra.hnf_cohen(B)
   @test is_upper_triangular(H)

   H, U = AbstractAlgebra.hnf_cohen_with_transform(B)
   @test is_upper_triangular(H)
   @test is_unit(det(U))
   @test U*B == H
end

@testset "Generic.MatAlg.hnf" begin
   R, x = polynomial_ring(QQ, "x")

   M = MatrixAlgebra(R, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5]))

   H = hnf(A)
   @test is_upper_triangular(H)

   H, U = hnf_with_transform(A)
   @test is_upper_triangular(H)
   @test is_unit(det(U))
   @test U*A == H

   # Fake up finite field of char 7, degree 2
   R, x = polynomial_ring(GF(7), "x")
   F, = residue_field(R, x^2 + 6x + 3)
   a = F(x)

   S, y = polynomial_ring(F, "y")

   N = MatrixAlgebra(S, 3)

   B = N(map(S, Any[1 0 a; a*y^3 0 3*a^2; y^4+a 0 y^2+y]))

   H = hnf(B)
   @test is_upper_triangular(H)

   H, U = hnf_with_transform(B)
   @test is_upper_triangular(H)
   @test is_unit(det(U))
   @test U*B == H
end

@testset "Generic.MatAlg.snf_kb" begin
   R, x = polynomial_ring(QQ, "x")

   M = MatrixAlgebra(R, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5]))

   T = AbstractAlgebra.snf_kb(A)
   @test is_snf(T)

   T, U, K = AbstractAlgebra.snf_kb_with_transform(A)
   @test is_snf(T)
   @test is_unit(det(U))
   @test is_unit(det(K))
   @test U*A*K == T

   # Fake up finite field of char 7, degree 2
   R, x = polynomial_ring(GF(7), "x")
   F, = residue_field(R, x^2 + 6x + 3)
   a = F(x)

   S, y = polynomial_ring(F, "y")

   N = MatrixAlgebra(S, 3)

   B = N(map(S, Any[1 0 a; a*y^3 0 3*a^2; y^4+a 0 y^2+y]))

   T = AbstractAlgebra.snf_kb(B)
   @test is_snf(T)

   T, U, K = AbstractAlgebra.snf_kb_with_transform(B)
   @test is_snf(T)
   @test is_unit(det(U))
   @test is_unit(det(K))
   @test U*B*K == T
end

@testset "Generic.MatAlg.snf" begin
   R, x = polynomial_ring(QQ, "x")

   M = MatrixAlgebra(R, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5]))

   T = snf(A)
   @test is_snf(T)

   T, U, K = snf_with_transform(A)
   @test is_snf(T)
   @test is_unit(det(U))
   @test is_unit(det(K))
   @test U*A*K == T

   # Fake up finite field of char 7, degree 2
   R, x = polynomial_ring(GF(7), "x")
   F, = residue_field(R, x^2 + 6x + 3)
   a = F(x)

   S, y = polynomial_ring(F, "y")

   N = MatrixAlgebra(S, 3)

   B = N(map(S, Any[1 0 a; a*y^3 0 3*a^2; y^4+a 0 y^2+y]))

   T = snf(B)
   @test is_snf(T)

   T, U, K = snf_with_transform(B)
   @test is_snf(T)
   @test is_unit(det(U))
   @test is_unit(det(K))
   @test U*B*K == T
end

@testset "Generic.MatAlg.$sim_zero" for sim_zero in (similar, zero)
   test_zero = sim_zero === zero
   for R = (ZZ, GF(11))
      M = MatrixAlgebra(R, rand(0:9))
      m = R == ZZ ? rand(M, -10:10) : rand(M)
      n = sim_zero(m)
      @test !test_zero || iszero(n)
      @test parent(n) == M
      @test size(n) == (nrows(M), ncols(M))
      r = rand(0:9)
      n = sim_zero(m, r)
      @test !test_zero || iszero(n)
      @test parent(n) == MatrixAlgebra(R, r)
      @test size(n) == (r, r)
      nn = sim_zero(m, r, r)
      @test !test_zero || iszero(nn)
      @test parent(nn) == MatrixAlgebra(R, r)
      @test size(nn) == (r, r)
      @test_throws ErrorException sim_zero(m, r, r+1)
      for S = [QQ, ZZ, GF(2), GF(5)]
         n = sim_zero(m, S)
         @test !test_zero || iszero(n)
         @test parent(n) == MatrixAlgebra(S, size(n)[1])
         @test size(n) == (nrows(M), ncols(M))
         r = rand(0:9)
         n = sim_zero(m, S, r)
         @test !test_zero || iszero(n)
         @test parent(n) == MatrixAlgebra(S, r)
         @test size(n) == (r, r)
         n = sim_zero(m, S, r, r)
         @test !test_zero || iszero(n)
         @test parent(n) == MatrixAlgebra(S, r)
         @test size(n) == (r, r)
         @test_throws ErrorException sim_zero(m, S, r, r+2)
      end
   end
end

@testset "Generic.MatAlg.change_base_ring" begin
   # Tests over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   U, x = polynomial_ring(R, "x")
   S = MatrixAlgebra(R, 2)
   
   M = rand(S, -10:10)

   N = change_base_ring(U, M)

   @test isa(N, MatAlgElem)
end

@testset "Generic.MatAlg.map" begin
   # Tests over noncommutative ring
   R = MatrixAlgebra(ZZ, 2)
   U, x = polynomial_ring(R, "x")
   S = MatrixAlgebra(U, 2)

   M = rand(R, -10:10)
   N = map(U, M)
   P = map(x->x^2, M)
   Q = map(S, M)

   @test isa(N, MatAlgElem)
   @test isa(P, MatAlgElem)
   @test isa(Q, MatAlgElem)
end

@testset "Generic.MatAlg.rand" begin
   M = MatrixAlgebra(ZZ, 3)

   test_rand(M, 1:9)

   M = MatrixAlgebra(GF(7), 2)

   test_rand(M)

   sp = Random.Sampler(MersenneTwister, M)
   @test parent(rand(sp)) == M
   v = rand(sp, 3)
   @test v isa Vector{elem_type(M)}
   @test all(x -> parent(x) == M, v)
end
