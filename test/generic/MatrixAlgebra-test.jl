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

function istriu(A::Generic.MatAlgElem)
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

function is_snf(A::Generic.MatAlgElem)
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

function is_weak_popov(P::Generic.MatAlgElem, rank::Int)
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

function test_gen_matalg_constructors()
   print("Generic.MatAlg.constructors...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixAlgebra(R, 3)

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

   @test_throws ErrorConstrDimMismatch S([t t^2 ; t^3 t^4])
   @test_throws ErrorConstrDimMismatch S([t t^2 t^3 ; t^4 t^5 t^6 ; t^7 t^8 t^9 ; t t^2 t^3])
   @test_throws ErrorConstrDimMismatch S([t, t^2])
   @test_throws ErrorConstrDimMismatch S([t, t^2, t^3, t^4, t^5, t^6, t^7, t^8, t^9, t^10])

   x = zero_matrix(R, 2, 2)
   y = zero_matrix(ZZ, 2, 3)

   @test x in [x, y]
   @test x in [y, x]
   @test !(x in [y])

   @test x in keys(Dict(x => 1))
   @test !(y in keys(Dict(x => 1)))

   println("PASS")
end

function test_gen_matalg_size()
   print("Generic.MatAlg.size...")

   S = MatrixAlgebra(QQ, 3)
   A = S([1 2 3; 4 5 6; 7 8 9])

   @test nrows(S) == ncols(S) == 3
   @test degree(S) == 3

   @test nrows(A) == ncols(A) == 3
   @test degree(A) == 3

   @test issquare(A)

   println("PASS")
end

function test_gen_matalg_manipulation()
   print("Generic.MatAlg.manipulation...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])

   @test isexact_type(typeof(A))

   @test iszero(zero(S))
   @test isone(one(S))

   U = MatrixAlgebra(QQ, 3)
   C = U([1 2 3; 5 6 7; 9 8 5])

   @test !isunit(A)
   @test isunit(C)

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

   @test iszero_row(C, 3)
   @test !iszero_row(C, 1)
   @test iszero_column(C, 2)
   @test !iszero_column(C, 1)

   println("PASS")
end

function test_gen_matalg_unary_ops()
   print("Generic.MatAlg.unary_ops...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([-t - 1 (-t) -R(1); -t^2 (-t) (-t); -R(-2) (-t - 2) (-t^2 - t - 1)])

   @test -A == B

   println("PASS")
end

function test_gen_matalg_binary_ops()
   print("Generic.MatAlg.binary_ops...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])

   @test A + B == S([t+3 t+3 R(2); t^2 + t 2*t+1 2*t+2; R(-3) t^2 + t + 2 t^3 + 1*t^2 + t + 1])

   @test A - B == S([t-1 t-3 R(0); t^2 - t R(-1) R(-2); R(-1) (-t^2 + t + 2) (-t^3 + t^2 + t + 1)])

   @test A*B == S([t^2 + 2*t + 1 2*t^2 + 4*t + 3 t^3 + t^2 + 3*t + 1; 3*t^2 - t (t^3 + 4*t^2 + t) t^4 + 2*t^2 + 2*t; t-5 t^4 + t^3 + 2*t^2 + 3*t - 4 t^5 + 1*t^4 + t^3 + t^2 + 4*t + 2])

   println("PASS")
end

function test_gen_matalg_adhoc_binary()
   print("Generic.MatAlg.adhoc_binary...")

   R, t = PolynomialRing(QQ, "t")
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

   println("PASS")
end

function test_gen_matalg_permutation()
   print("Generic.MatAlg.permutation...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   T = PermutationGroup(3)
   P = T([2, 3, 1])

   @test A == inv(P)*(P*A)

   println("PASS")
end

function test_gen_matalg_comparison()
   print("Generic.MatAlg.comparison...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
   B = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test A == B

   @test A != one(S)

   println("PASS")
end

function test_gen_matalg_adhoc_comparison()
   print("Generic.MatAlg.adhoc_comparison...")

   R, t = PolynomialRing(QQ, "t")
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

   println("PASS")
end

function test_gen_matalg_powering()
   print("Generic.MatAlg.powering...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test A^5 == A^2*A^3

   @test A^0 == one(S)

   println("PASS")
end

function test_gen_matalg_exact_division()
   print("Generic.MatAlg.exact_division...")

   S = MatrixAlgebra(QQ, 3)

   M = rand(S, -20:20)
   N = randmat_with_rank(S, 3, -20:20)

   @test divexact_right(M*N, N) == M
   @test divexact_left(N*M, N) == M

   println("PASS")
end

function test_gen_matalg_adhoc_exact_division()
   print("Generic.MatAlg.adhoc_exact_division...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test divexact(5*A, 5) == A
   @test divexact(12*A, BigInt(12)) == A
   @test divexact(12*A, Rational{BigInt}(12)) == A
   @test divexact((1 + t)*A, 1 + t) == A

   println("PASS")
end

function test_gen_matalg_transpose()
   print("Generic.MatAlg.transpose...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixAlgebra(R, 3)
   arr = [t + 1 t R(1); t^2 t t; t+1 t^2 R(-1)]
   A = S(arr)
   B = S(permutedims(arr, [2, 1]))
   @test transpose(A) == B

   println("PASS")
end

function test_gen_matalg_gram()
   print("Generic.MatAlg.gram...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test gram(A) == S([2*t^2 + 2*t + 2 t^3 + 2*t^2 + t 2*t^2 + t - 1; t^3 + 2*t^2 + t t^4 + 2*t^2 t^3 + 3*t; 2*t^2 + t - 1 t^3 + 3*t t^4 + 2*t^3 + 4*t^2 + 6*t + 9])

   println("PASS")
end

function test_gen_matalg_tr()
   print("Generic.MatAlg.tr...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test tr(A) == t^2 + 3t + 2

   println("PASS")
end

function test_gen_matalg_content()
   print("Generic.MatAlg.content...")

   R, t = PolynomialRing(QQ, "t")
   S = MatrixAlgebra(R, 3)

   A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

   @test content((1 + t)*A) == 1 + t
   println("PASS")
end

function test_gen_matalg_lu()
   print("Generic.MatAlg.lu...")

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
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

   R, z = PolynomialRing(ZZ, "z")
   F = FractionField(R)
   S = MatrixAlgebra(F, 3)

   A = S([F(0), F(0), F(11), 78*z^3-102*z^2+48*z+12, F(92), -16*z^2+80*z-149, -377*z^3+493*z^2-232*z-58, F(-448), 80*z^2-385*z+719])

   r, P, L, U = lu(A)

   @test r == 3
   @test P*A == L*U

   println("PASS")
end

function test_gen_matalg_fflu()
   print("Generic.MatAlg.fflu...")

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   S = MatrixAlgebra(K, 3)

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

   S = MatrixAlgebra(QQ, 3)
   A = S([0, 0, 1, 12, 1, 11, 1, 0, 1])

   r, d, P, L, U, = fflu(A)

   D = S()
   D[1, 1] = inv(U[1, 1])
   D[2, 2] = inv(U[1, 1]*U[2, 2])
   D[3, 3] = inv(U[2, 2])
   @test r == 3
   @test P*A == L*D*U

   println("PASS")
end

function test_gen_matalg_det()
   print("Generic.MatAlg.det...")

   S, x = PolynomialRing(ResidueRing(ZZ, 1009*2003), "x")

   for dim = 0:5
      R = MatrixAlgebra(S, dim)

      M = rand(R, 0:5, -100:100)

      @test det(M) == AbstractAlgebra.det_clow(M)
   end

   S, z = PolynomialRing(ZZ, "z")

   for dim = 0:5
      R = MatrixAlgebra(S, dim)

      M = rand(R, 0:3, -20:20)

      @test det(M) == AbstractAlgebra.det_clow(M)
   end

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")

   for dim = 0:7
      S = MatrixAlgebra(K, dim)

      M = rand(S, 0:2, -100:100)

      @test det(M) == AbstractAlgebra.det_clow(M)
   end

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   for dim = 0:5
      T = MatrixAlgebra(S, dim)
      M = rand(T, 0:2, 0:2, -10:10)

      @test det(M) == AbstractAlgebra.det_clow(M)
   end

   println("PASS")
end

function test_gen_matalg_rank()
   print("Generic.MatAlg.rank...")

   S = ResidueRing(ZZ, 20011*10007)
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

   S, z = PolynomialRing(ZZ, "z")
   R = MatrixAlgebra(S, 4)

   M = R([S(-2) S(0) S(5) S(3); 5*z^2+5*z-5 S(0) S(-z^2+z) 5*z^2+5*z+1; 2*z-1 S(0) z^2+3*z+2 S(-4*z); 3*z-5 S(0) S(-5*z+5) S(1)])

   @test rank(M) == 3

   R = MatrixSpace(S, 5, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, 0:3, -20:20)

      @test rank(M) == i
   end

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   S = MatrixAlgebra(K, 3)

   M = S([a a^2 + 2*a - 1 2*a^2 - 1*a; 2*a+2 2*a^2 + 2*a (-2*a^2 - 2*a); (-a) (-a^2) a^2])

   @test rank(M) == 2

   S = MatrixAlgebra(K, 5)

   for i = 0:5
      M = randmat_with_rank(S, i, 0:2, -100:100)

      @test rank(M) == i
   end

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
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

   println("PASS")
end

function test_gen_matalg_solve_lu()
   print("Generic.MatAlg.solve_lu...")

   S = QQ

   for dim = 0:5
      R = MatrixAlgebra(S, dim)
      U = MatrixAlgebra(S, dim)

      M = randmat_with_rank(R, dim, -100:100)
      b = rand(U, -100:100)

      x = Generic.solve_lu(M, b)

      @test M*x == b
   end

   S, y = PolynomialRing(ZZ, "y")
   K = FractionField(S)

   for dim = 0:5
      R = MatrixAlgebra(S, dim)
      U = MatrixAlgebra(S, dim)
      T = MatrixAlgebra(K, dim)

      M = randmat_with_rank(R, dim, 0:5, -100:100)
      b = rand(U, 0:5, -100:100);

      MK = T(elem_type(K)[ K(M[i, j]) for i in 1:nrows(M), j in 1:ncols(M) ])
      bK = T(elem_type(K)[ K(b[i, j]) for i in 1:nrows(b), j in 1:ncols(b) ])

      x = Generic.solve_lu(MK, bK)

      @test MK*x == bK
   end

   println("PASS")
end

function test_gen_matalg_rref()
   print("Generic.MatAlg.rref...")

   S = ResidueRing(ZZ, 20011*10007)
   R = MatrixAlgebra(S, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, -100:100)

      do_test = false
      r = 0
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
   R = MatrixAlgebra(S, 5)

   for i = 0:5
      M = randmat_with_rank(R, i, 0:3, -20:20)

      r, d, A = rref(M)

      @test r == i
      @test isrref(A)
   end

   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   S = MatrixAlgebra(K, 5)

   for i = 0:5
      M = randmat_with_rank(S, i, 0:2, -100:100)

      r, A = rref(M)

      @test r == i
      @test isrref(A)
   end

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   T = MatrixAlgebra(S, 5)

   for i = 0:5
      M = randmat_with_rank(T, i, 0:2, 0:2, -20:20)

      r, d, A = rref(M)

      @test r == i
      @test isrref(A)
   end

   println("PASS")
end

function test_gen_matalg_inversion()
@testset "Generic.MatAlg.inversion..." begin

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

   @testset "Matrix Algebra over ResidueRing" begin
   S = ResidueRing(ZZ, 20011*10007)

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
   S, x = PolynomialRing(ZZ, "x")

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

   @testset "Matrix Algebra over NumberField over QQ" begin
   R, x = PolynomialRing(QQ, "x")
   S, a = NumberField(x^3 + 3x + 1, "a")

   for dim = 1:5
      R = MatrixAlgebra(S, dim)

      M = randmat_with_rank(R, dim, 0:2, -100:100)

      X = inv(M)

      @test isone(M*X)
   end
   end

   @testset "Matrix Algebra over (ZZ[x])[y]" begin
   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

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
end # of @testset "Generic.MatAlg.inversion..."
end

function test_gen_matalg_hessenberg()
   print("Generic.MatAlg.hessenberg...")

   R = ResidueRing(ZZ, 18446744073709551629)

   for dim = 0:5
      S = MatrixAlgebra(R, dim)
      U, x = PolynomialRing(R, "x")

      for i = 1:10
         M = rand(S, -5:5)

         A = hessenberg(M)

         @test ishessenberg(A)
      end
   end

   println("PASS")
end

function test_gen_matalg_charpoly()
   print("Generic.MatAlg.charpoly...")

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

function test_gen_matalg_minpoly()
   print("Generic.MatAlg.minpoly...")

   R = GF(103)
   T, y = PolynomialRing(R, "y")
   S = MatrixAlgebra(R, 3)

   M = S([92 97 8;
          0 5 13;
          0 16 2])

   @test minpoly(T, M) == y^2+96*y+8

   R = GF(3)
   T, y = PolynomialRing(R, "y")
   S = MatrixAlgebra(R, 4)

   M = S([1 2 0 2;
         1 2 1 0;
         1 2 2 1;
         2 1 2 0])

   @test minpoly(T, M) == y^2 + 2y

   R = GF(13)
   T, y = PolynomialRing(R, "y")
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

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   U, z = PolynomialRing(S, "z")
   T = MatrixAlgebra(S, 6)

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
   T = MatrixAlgebra(R, 6)

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

function test_gen_matalg_row_swapping()
   print("Generic.MatAlg.row_swapping...")

   R, x = PolynomialRing(ZZ, "x")
   M = MatrixAlgebra(R, 3)

   a = M(map(R, [1 2 3; 4 5 6; 7 8 9]))

   @test swap_rows(a, 1, 3) == M(map(R, [7 8 9; 4 5 6; 1 2 3]))

   swap_rows!(a, 2, 3)

   @test a == M(map(R, [1 2 3; 7 8 9; 4 5 6]))

   println("PASS")
end

function test_gen_matalg_hnf_minors()
  print("Generic.MatAlg.hnf_minors...")

   R, x = PolynomialRing(QQ, "x")

   M = MatrixAlgebra(R, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5]))

   H = hnf_minors(A)
   @test istriu(H)

   H, U = hnf_minors_with_transform(A)
   @test istriu(H)
   @test isunit(det(U))
   @test U*A == H

   # Fake up finite field of char 7, degree 2
   R, x = PolynomialRing(GF(7), "x")
   F = ResidueField(R, x^2 + 6x + 3)
   a = F(x)

   S, y = PolynomialRing(F, "y")

   N = MatrixAlgebra(S, 4)

   B = N(map(S, Any[1 0 a 0; a*y^3 0 3*a^2 0; y^4+a 0 y^2+y 5; y 1 y 2]))

   H = hnf_minors(B)
   @test istriu(H)

   H, U = hnf_minors_with_transform(B)
   @test istriu(H)
   @test isunit(det(U))
   @test U*B == H

   println("PASS")
end

function test_gen_matalg_hnf_kb()
   print("Generic.MatAlg.hnf_kb...")

   R, x = PolynomialRing(QQ, "x")

   M = MatrixAlgebra(R, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5]))

   H = AbstractAlgebra.hnf_kb(A)
   @test istriu(H)

   H, U = AbstractAlgebra.hnf_kb_with_transform(A)
   @test istriu(H)
   @test isunit(det(U))
   @test U*A == H

   # Fake up finite field of char 7, degree 2
   R, x = PolynomialRing(GF(7), "x")
   F = ResidueField(R, x^2 + 6x + 3)
   a = F(x)

   S, y = PolynomialRing(F, "y")

   N = MatrixAlgebra(S, 3)

   B = N(map(S, Any[1 0 a; a*y^3 0 3*a^2; y^4+a 0 y^2+y]))

   H = AbstractAlgebra.hnf_kb(B)
   @test istriu(H)

   H, U = AbstractAlgebra.hnf_kb_with_transform(B)
   @test istriu(H)
   @test isunit(det(U))
   @test U*B == H

   println("PASS")
end

function test_gen_matalg_hnf_cohen()
   print("Generic.MatAlg.hnf_cohen...")

   R, x = PolynomialRing(QQ, "x")

   M = MatrixAlgebra(R, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5]))

   H = AbstractAlgebra.hnf_cohen(A)
   @test istriu(H)

   H, U = AbstractAlgebra.hnf_cohen_with_transform(A)
   @test istriu(H)
   @test isunit(det(U))
   @test U*A == H

   # Fake up finite field of char 7, degree 2
   R, x = PolynomialRing(GF(7), "x")
   F = ResidueField(R, x^2 + 6x + 3)
   a = F(x)

   S, y = PolynomialRing(F, "y")

   N = MatrixAlgebra(S, 3)

   B = N(map(S, Any[1 0 a; a*y^3 0 3*a^2; y^4+a 0 y^2+y]))

   H = AbstractAlgebra.hnf_cohen(B)
   @test istriu(H)

   H, U = AbstractAlgebra.hnf_cohen_with_transform(B)
   @test istriu(H)
   @test isunit(det(U))
   @test U*B == H

   println("PASS")
end

function test_gen_matalg_hnf()
   print("Generic.MatAlg.hnf...")

   R, x = PolynomialRing(QQ, "x")

   M = MatrixAlgebra(R, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5]))

   H = hnf(A)
   @test istriu(H)

   H, U = hnf_with_transform(A)
   @test istriu(H)
   @test isunit(det(U))
   @test U*A == H

   # Fake up finite field of char 7, degree 2
   R, x = PolynomialRing(GF(7), "x")
   F = ResidueField(R, x^2 + 6x + 3)
   a = F(x)

   S, y = PolynomialRing(F, "y")

   N = MatrixAlgebra(S, 3)

   B = N(map(S, Any[1 0 a; a*y^3 0 3*a^2; y^4+a 0 y^2+y]))

   H = hnf(B)
   @test istriu(H)

   H, U = hnf_with_transform(B)
   @test istriu(H)
   @test isunit(det(U))
   @test U*B == H

   println("PASS")
end

function test_gen_matalg_snf_kb()
   print("Generic.MatAlg.snf_kb...")

   R, x = PolynomialRing(QQ, "x")

   M = MatrixAlgebra(R, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5]))

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

   N = MatrixAlgebra(S, 3)

   B = N(map(S, Any[1 0 a; a*y^3 0 3*a^2; y^4+a 0 y^2+y]))

   T = AbstractAlgebra.snf_kb(B)
   @test is_snf(T)

   T, U, K = AbstractAlgebra.snf_kb_with_transform(B)
   @test is_snf(T)
   @test isunit(det(U))
   @test isunit(det(K))
   @test U*B*K == T

   println("PASS")
end

function test_gen_matalg_snf()
   print("Generic.MatAlg.snf...")

   R, x = PolynomialRing(QQ, "x")

   M = MatrixAlgebra(R, 3)

   A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5]))

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

   N = MatrixAlgebra(S, 3)

   B = N(map(S, Any[1 0 a; a*y^3 0 3*a^2; y^4+a 0 y^2+y]))

   T = snf(B)
   @test is_snf(T)

   T, U, K = snf_with_transform(B)
   @test is_snf(T)
   @test isunit(det(U))
   @test isunit(det(K))
   @test U*B*K == T

   println("PASS")
end

function test_gen_matalg()
   test_gen_matalg_constructors()
   test_gen_matalg_size()
   test_gen_matalg_manipulation()
   test_gen_matalg_unary_ops()
   test_gen_matalg_binary_ops()
   test_gen_matalg_adhoc_binary()
   test_gen_matalg_permutation()
   test_gen_matalg_comparison()
   test_gen_matalg_adhoc_comparison()
   test_gen_matalg_powering()
   test_gen_matalg_exact_division()
   test_gen_matalg_adhoc_exact_division()
   test_gen_matalg_transpose()
   test_gen_matalg_gram()
   test_gen_matalg_tr()
   test_gen_matalg_content()
   test_gen_matalg_lu()
   test_gen_matalg_fflu()
   test_gen_matalg_det()
   test_gen_matalg_rank()
   test_gen_matalg_solve_lu()
   test_gen_matalg_rref()
   test_gen_matalg_inversion()
   test_gen_matalg_hessenberg()
   test_gen_matalg_charpoly()
   test_gen_matalg_minpoly()
   test_gen_matalg_row_swapping()
#   test_gen_matalg_hnf_minors() # see bug 160
   test_gen_matalg_hnf_kb()
   test_gen_matalg_hnf_cohen()
   test_gen_matalg_hnf()
   test_gen_matalg_snf_kb()
   test_gen_matalg_snf()

   println("")
end
