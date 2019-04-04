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

function test_gen_mat_constructors()
   print("Generic.Mat.constructors...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   @test elem_type(S) == Generic.Mat{elem_type(R)}
   @test elem_type(Generic.MatSpace{elem_type(R)}) == Generic.Mat{elem_type(R)}
   @test parent_type(Generic.Mat{elem_type(R)}) == Generic.MatSpace{elem_type(R)}
   @test base_ring(S) == R
   @test nrows(S) == 3
   @test ncols(S) == 3

   @test typeof(S) <: Generic.MatSpace

   f = S(t^2 + 1)

   @test isa(f, MatElem)

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

   arr = [1 2; 3 4]
   arr2 = [1, 2, 3, 4, 5, 6]

   for T in [R, Int, BigInt, Rational{Int}, Rational{BigInt}]
      M = matrix(R, map(T, arr))
      @test isa(M, Generic.Mat{elem_type(R)})
      @test M.base_ring == R
      @test nrows(M) == 2
      @test ncols(M) == 2

      M2 = matrix(R, 2, 3, map(T, arr2))
      @test isa(M2, Generic.Mat{elem_type(R)})
      @test M2.base_ring == R
      @test nrows(M2) == 2
      @test ncols(M2) == 3
      @test_throws ErrorConstrDimMismatch matrix(R, 2, 2, map(T, arr2))
      @test_throws ErrorConstrDimMismatch matrix(R, 2, 4, map(T, arr2))
   end

   M = matrix(R, arr')
   @test isa(M, Generic.Mat{elem_type(R)})

   M = matrix(R, 2, 3, view(arr2', 1:6))
   @test isa(M, Generic.Mat{elem_type(R)})

   M3 = zero_matrix(R, 2, 3)

   @test isa(M3, Generic.Mat{elem_type(R)})
   @test M3.base_ring == R

   M4 = identity_matrix(R, 3)

   @test isa(M4, Generic.Mat{elem_type(R)})
   @test M4.base_ring == R

   x = zero_matrix(R, 2, 2)
   y = zero_matrix(ZZ, 2, 3)

   @test x in [x, y]
   @test x in [y, x]
   @test !(x in [y])

   @test x in keys(Dict(x => 1))
   @test !(y in keys(Dict(x => 1)))

   println("PASS")
end

function test_gen_mat_size()
   print("Generic.Mat.size...")

   A = matrix(QQ, [1 2 3; 4 5 6; 7 8 9])
   B = matrix(QQ, [1 2 3 4; 5 6 7 8])

   @test size(A) == (3,3)
   @test issquare(A)

   @test size(B) == (2,4)
   @test !issquare(B)

   println("PASS")
end

function test_gen_mat_manipulation()
   print("Generic.Mat.manipulation...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])

   @test iszero(zero(S))
   @test isone(one(S))

   B[1, 1] = R(3)
   @test B[1, 1] == R(3)

   B[1, 1] = 4
   @test B[1, 1] == R(4)

   B[1, 1] = BigInt(5)
   @test B[1, 1] == R(5)

   @test nrows(B) == 3
   @test ncols(B) == 3

   @test deepcopy(A) == A

   println("PASS")
end

function test_gen_mat_unary_ops()
   print("Generic.Mat.unary_ops...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([-t - 1 (-t) -R(1); -t^2 (-t) (-t); -R(-2) (-t - 2) (-t^2 - t - 1)])

   @test -A == B

   println("PASS")
end

function test_gen_mat_sub()
   print("Generic.Mat.sub...")

   S = MatrixSpace(ZZ, 3, 3)

   A = S([1 2 3; 4 5 6; 7 8 9])

   B = @inferred sub(A, 1, 1, 2, 2)

   @test typeof(B) == typeof(A)
   @test B == MatrixSpace(ZZ, 2, 2)([1 2; 4 5])

   B[1, 1] = 10
   @test A == S([1 2 3; 4 5 6; 7 8 9])

   C = @inferred sub(B, 1:2, 1:2)

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

   println("PASS")
end

function test_gen_mat_binary_ops()
   print("Generic.Mat.binary_ops...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])

   @test A + B == S([t+3 t+3 R(2); t^2 + t 2*t+1 2*t+2; R(-3) t^2 + t + 2 t^3 + 1*t^2 + t + 1])

   @test A - B == S([t-1 t-3 R(0); t^2 - t R(-1) R(-2); R(-1) (-t^2 + t + 2) (-t^3 + t^2 + t + 1)])

   @test A*B == S([t^2 + 2*t + 1 2*t^2 + 4*t + 3 t^3 + t^2 + 3*t + 1; 3*t^2 - t (t^3 + 4*t^2 + t) t^4 + 2*t^2 + 2*t; t-5 t^4 + t^3 + 2*t^2 + 3*t - 4 t^5 + 1*t^4 + t^3 + t^2 + 4*t + 2])

   println("PASS")
end

function test_gen_mat_adhoc_binary()
   print("Generic.Mat.adhoc_binary...")

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

   println("PASS")
end

function test_gen_mat_permutation()
   print("Generic.Mat.permutation...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   T = PermutationGroup(3)
   P = T([2, 3, 1])

   @test A == inv(P)*(P*A)

   println("PASS")
end

function test_gen_mat_comparison()
   print("Generic.Mat.comparison...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test A == B

   @test A != one(S)

   println("PASS")
end

function test_gen_mat_adhoc_comparison()
   print("Generic.Mat.adhoc_comparison...")

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

   println("PASS")
end

function test_gen_mat_powering()
   print("Generic.Mat.powering...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test A^5 == A^2*A^3

   @test A^0 == one(S)

   println("PASS")
end

function test_gen_mat_adhoc_exact_division()
   print("Generic.Mat.adhoc_exact_division...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test divexact(5*A, 5) == A
   @test divexact(12*A, BigInt(12)) == A
   @test divexact(12*A, Rational{BigInt}(12)) == A
   @test divexact((1 + t)*A, 1 + t) == A

   println("PASS")
end

function test_gen_mat_transpose()
   print("Generic.Mat.transpose...")

   R, t = PolynomialRing(QQ, "t")
   arr = [t + 1 t R(1); t^2 t t]
   A = matrix(R, arr)
   B = matrix(R, permutedims(arr, [2, 1]))
   @test transpose(A) == B

   arr = [t + 1 t; t^2 t]
   A = matrix(R, arr)
   B = matrix(R, permutedims(arr, [2, 1]))
   @test transpose(A) == B

   println("PASS")
end

function test_gen_mat_gram()
   print("Generic.Mat.gram...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test gram(A) == S([2*t^2 + 2*t + 2 t^3 + 2*t^2 + t 2*t^2 + t - 1; t^3 + 2*t^2 + t t^4 + 2*t^2 t^3 + 3*t; 2*t^2 + t - 1 t^3 + 3*t t^4 + 2*t^3 + 4*t^2 + 6*t + 9])

   println("PASS")
end

function test_gen_mat_tr()
   print("Generic.Mat.tr...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test tr(A) == t^2 + 3t + 2

   println("PASS")
end

function test_gen_mat_content()
   print("Generic.Mat.content...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixSpace(R, 3, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test content((1 + t)*A) == 1 + t
   println("PASS")
end

function test_gen_mat_lu()
   print("Generic.Mat.lu...")

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

   println("PASS")
end

function test_gen_mat_fflu()
   print("Generic.Mat.fflu...")

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

   println("PASS")
end

function test_gen_mat_det()
   print("Generic.Mat.det...")

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

   println("PASS")
end

function test_gen_mat_minors()
   print("Generic.Mat.minors...")

   S,z=PolynomialRing(ZZ,"z")
   n = 5
   R = MatrixSpace(S,n,n)
   for r=0:n
      M = randmat_with_rank(R, r, 0:3, 0:3)
      @test [1] == minors(M,0)
      for i=r+1:n
         for m in minors(M,i)
            @test m==0
         end
      end
      @test [det(M)] == minors(M,n)
      @test [] == minors(M,n+1)
   end

   println("PASS")
end

function test_gen_mat_rank()
   print("Generic.Mat.rank...")

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

   println("PASS")
end

function test_gen_mat_solve_lu()
   print("Generic.Mat.solve_lu...")

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

   println("PASS")
end

function test_gen_mat_solve_rational()
   print("Generic.Mat.solve_rational...")

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

   println("PASS")
end

function test_gen_mat_solve_triu()
   print("Generic.Mat.solve_triu...")

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

   println("PASS")
end

function test_gen_mat_rref()
   print("Generic.Mat.rref...")

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

   println("PASS")
end

function test_gen_mat_nullspace()
   print("Generic.Mat.nullspace...")

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

   println("PASS")
end

function test_gen_mat_inversion()
   print("Generic.Mat.inversion...")

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
      @test inv(c) isa eltype(M)
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

   println("PASS")
end

function test_gen_mat_hessenberg()
   print("Generic.Mat.hessenberg...")

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

   println("PASS")
end

function test_gen_mat_kronecker_product()
   print("Generic.Mat.kronecker_product...")

   R = ResidueRing(ZZ, 18446744073709551629)
   S = MatrixSpace(R, 2, 3)
   S2 = MatrixSpace(R, 2, 2)
   S3 = MatrixSpace(R, 3, 3)

   A = S(R.([2 3 5; 9 6 3]))
   B = S2(R.([2 3; 1 4]))
   C = S3(R.([2 3 5; 1 4 7; 9 6 3]))

   @test size(kronecker_product(A, A)) == (4,9)
   @test kronecker_product(B*A,A*C) == kronecker_product(B,A) * kronecker_product(A,C)

   println("PASS")
end

function test_gen_mat_charpoly()
   print("Generic.Mat.charpoly...")

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

   println("PASS")
end

function test_gen_mat_minpoly()
   print("Generic.Mat.minpoly...")

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

   println("PASS")
end

function test_gen_mat_row_col_swapping()
   print("Generic.Mat.row_col_swapping...")

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
   @test invert_rows(a) == matrix(R, [3 4; 1 2])
   invert_rows!(a)
   @test a == matrix(R, [3 4; 1 2])

   a = matrix(R, [1 2; 3 4])
   @test invert_cols(a) == matrix(R, [2 1; 4 3])
   invert_cols!(a)
   @test a == matrix(R, [2 1; 4 3])

   a = matrix(R, [1 2 3; 3 4 5; 5 6 7])

   @test invert_rows(a) == matrix(R, [5 6 7; 3 4 5; 1 2 3])
   invert_rows!(a)
   @test a == matrix(R, [5 6 7; 3 4 5; 1 2 3])

   a = matrix(R, [1 2 3; 3 4 5; 5 6 7])
   @test invert_cols(a) == matrix(R, [3 2 1; 5 4 3; 7 6 5])
   invert_cols!(a)
   @test a == matrix(R, [3 2 1; 5 4 3; 7 6 5])

   println("PASS")
end

function test_gen_mat_concat()
   print("Generic.Mat.concat...")

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

   println("PASS")
end

function test_gen_mat_hnf_minors()
  print("Generic.Mat.hnf_minors...")

   R, x = PolynomialRing(QQ, "x")

   M = MatrixSpace(R, 4, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))

   H = hnf_minors(A)
   @test istriu(H)

   H, U = hnf_minors_with_trafo(A)
   @test istriu(H)
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
   @test istriu(H)

   H, U = hnf_minors_with_trafo(B)
   @test istriu(H)
   @test isunit(det(U))
   @test U*B == H

   println("PASS")
end

function test_gen_mat_hnf_kb()
   print("Generic.Mat.hnf_kb...")

   R, x = PolynomialRing(QQ, "x")

   M = MatrixSpace(R, 4, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))

   H = AbstractAlgebra.hnf_kb(A)
   @test istriu(H)

   H, U = AbstractAlgebra.hnf_kb_with_trafo(A)
   @test istriu(H)
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
   @test istriu(H)

   H, U = AbstractAlgebra.hnf_kb_with_trafo(B)
   @test istriu(H)
   @test isunit(det(U))
   @test U*B == H

   println("PASS")
end

function test_gen_mat_hnf_cohen()
   print("Generic.Mat.hnf_cohen...")

   R, x = PolynomialRing(QQ, "x")

   M = MatrixSpace(R, 4, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))

   H = AbstractAlgebra.hnf_cohen(A)
   @test istriu(H)

   H, U = AbstractAlgebra.hnf_cohen_with_trafo(A)
   @test istriu(H)
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
   @test istriu(H)

   H, U = AbstractAlgebra.hnf_cohen_with_trafo(B)
   @test istriu(H)
   @test isunit(det(U))
   @test U*B == H

   println("PASS")
end

function test_gen_mat_hnf()
   print("Generic.Mat.hnf...")

   R, x = PolynomialRing(QQ, "x")

   M = MatrixSpace(R, 4, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))

   H = hnf(A)
   @test istriu(H)

   H, U = hnf_with_trafo(A)
   @test istriu(H)
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
   @test istriu(H)

   H, U = hnf_with_trafo(B)
   @test istriu(H)
   @test isunit(det(U))
   @test U*B == H

   println("PASS")
end

function test_gen_mat_snf_kb()
   print("Generic.Mat.snf_kb...")

   R, x = PolynomialRing(QQ, "x")

   M = MatrixSpace(R, 4, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))

   T = AbstractAlgebra.snf_kb(A)
   @test is_snf(T)

   T, U, K = AbstractAlgebra.snf_kb_with_trafo(A)
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

   T, U, K = AbstractAlgebra.snf_kb_with_trafo(B)
   @test is_snf(T)
   @test isunit(det(U))
   @test isunit(det(K))
   @test U*B*K == T

   println("PASS")
end

function test_gen_mat_snf()
   print("Generic.Mat.snf...")

   R, x = PolynomialRing(QQ, "x")

   M = MatrixSpace(R, 4, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))

   T = snf(A)
   @test is_snf(T)

   T, U, K = snf_with_trafo(A)
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

   T, U, K = snf_with_trafo(B)
   @test is_snf(T)
   @test isunit(det(U))
   @test isunit(det(K))
   @test U*B*K == T

   println("PASS")
end

function test_gen_mat_weak_popov()
   print("Generic.Mat.weak_popov...")

   R, x = PolynomialRing(QQ, "x")

   A = matrix(R, map(R, Any[1 2 3 x; x 2*x 3*x x^2; x x^2+1 x^3+x^2 x^4+x^2+1]))
   r = rank(A)

   P = weak_popov(A)
   @test is_weak_popov(P, r)

   P, U = weak_popov_with_trafo(A)
   @test is_weak_popov(P, r)
   @test U*A == P
   @test isunit(det(U))

   F = GF(7)

   S, y = PolynomialRing(F, "y")

   B = matrix(S, map(S, Any[ 4*y^2+3*y+5 4*y^2+3*y+4 6*y^2+1; 3*y+6 3*y+5 y+3; 6*y^2+4*y+2 6*y^2 2*y^2+y]))
   s = rank(B)

   P = weak_popov(B)
   @test is_weak_popov(P, s)

   P, U = weak_popov_with_trafo(B)
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

      P, U = weak_popov_with_trafo(A)
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

      P, U = weak_popov_with_trafo(A)
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

      P, U = weak_popov_with_trafo(A)
      @test is_weak_popov(P, r)
      @test U*A == P
      @test isunit(det(U))
   end

   println("PASS")
end

function test_gen_mat()
   test_gen_mat_constructors()
   test_gen_mat_size()
   test_gen_mat_manipulation()
   test_gen_mat_sub()
   test_gen_mat_unary_ops()
   test_gen_mat_binary_ops()
   test_gen_mat_adhoc_binary()
   test_gen_mat_permutation()
   test_gen_mat_comparison()
   test_gen_mat_adhoc_comparison()
   test_gen_mat_powering()
   test_gen_mat_adhoc_exact_division()
   test_gen_mat_transpose()
   test_gen_mat_gram()
   test_gen_mat_tr()
   test_gen_mat_content()
   test_gen_mat_lu()
   test_gen_mat_fflu()
   test_gen_mat_det()
   test_gen_mat_rank()
   test_gen_mat_solve_lu()
   test_gen_mat_solve_rational()
   test_gen_mat_solve_triu()
   test_gen_mat_rref()
   test_gen_mat_nullspace()
   test_gen_mat_inversion()
   test_gen_mat_hessenberg()
   test_gen_mat_kronecker_product()
   test_gen_mat_charpoly()
   test_gen_mat_minpoly()
   test_gen_mat_row_col_swapping()
   test_gen_mat_concat()
   test_gen_mat_hnf_minors()
   test_gen_mat_hnf_kb()
   test_gen_mat_hnf_cohen()
   test_gen_mat_hnf()
   test_gen_mat_snf_kb()
   test_gen_mat_snf()
   test_gen_mat_weak_popov()
   test_gen_mat_minors()

   println("")
end
