@testset "Generic.NCPoly.constructors" begin
   R = matrix_ring(ZZ, 2)
   S1 = polynomial_ring(R, "y")
   S2 = R["y"]

   @test polynomial_ring(R, "y", cached = true)[1] === polynomial_ring(R, "y", cached = true)[1]
   @test polynomial_ring(R, "y", cached = true)[1] !== polynomial_ring(R, "y", cached = false)[1]

   for (S, y) in (S1, S2)
      @test base_ring(S) === R
      @test coefficient_ring(S) === R
      @test coefficient_ring_type(S) === typeof(R)
      @test elem_type(S) == Generic.NCPoly{elem_type(R)}
      @test elem_type(Generic.NCPolyRing{elem_type(R)}) == Generic.NCPoly{elem_type(R)}
      @test parent_type(Generic.NCPoly{elem_type(R)}) == Generic.NCPolyRing{elem_type(R)}

      @test S isa Generic.NCPolyRing

      @test isa(y, NCPolyRingElem)
   end

   S3, _ = R["x"]["y"]
   @test elem_type(S3) == Generic.NCPoly{Generic.NCPoly{Generic.MatRingElem{BigInt}}}
   @test typeof(S3) == Generic.NCPolyRing{Generic.NCPoly{Generic.MatRingElem{BigInt}}}

   S4, _ = R["x"]["y"]["z"]
   @test elem_type(S4) == Generic.NCPoly{Generic.NCPoly{Generic.NCPoly{Generic.MatRingElem{BigInt}}}}
   @test typeof(S4) == Generic.NCPolyRing{Generic.NCPoly{Generic.NCPoly{Generic.MatRingElem{BigInt}}}}

   S, y = S1
   T, z = polynomial_ring(S, "z")

   @test T isa Generic.NCPolyRing

   @test isa(z, NCPolyRingElem)

   f = one(R) + y^3 + z + 1

   @test isa(f, NCPolyRingElem)

   g = S(2)

   @test isa(g, NCPolyRingElem)

   h = S(rand(R, -10:10))

   @test isa(h, NCPolyRingElem)

   j = T(rand(R, -10:10))

   @test isa(j, NCPolyRingElem)

   k = S([rand(R, -10:10) for i in 1:3])

   @test isa(k, NCPolyRingElem)

   l = S(k)

   @test isa(l, NCPolyRingElem)

   m = S([1, 2, 3])

   @test isa(m, NCPolyRingElem)

   n = S([ZZ(1), ZZ(2), ZZ(3)])

   @test isa(n, NCPolyRingElem)
end

@testset "Generic.NCPoly.manipulation" begin
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")

   @test iszero(zero(S))

   @test isone(one(S))

   @test is_gen(gen(S))

   f = 2*y + 1

   @test leading_coefficient(f) == 2

   @test trailing_coefficient(2y^2 + 3y) == 3

   @test degree(f) == 1

   @test one(f) == one(S)
   @test zero(f) == zero(S)

   h = y^2 + 3*y + 3

   @test coeff(h, 1) == 3

   @test length(h) == 3

   @test canonical_unit(-y + 1) == -1

   @test deepcopy(h) == h

   @test is_term(rand(R, -10:10)*y^2)

   @test !is_monomial(2*y^2)

   @test is_monomial(y^2)

   @test !is_monomial(y^2 + y + 1)
end

@testset "Generic.NCPoly.rand" begin
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")

   test_rand(S, -1:10, -10:10)
end

@testset "Generic.NCPoly.binary_ops" begin
   #  Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")

   for iter = 1:100
      f = rand(S, -1:10, -10:10)
      g = rand(S, -1:10, -10:10)
      h = rand(S, -1:10, -10:10)
      @test f + g == g + f
      @test f + (g + h) == (f + g) + h
      @test f*(g + h) == f*g + f*h
      @test (f - h) + (g + h) == f + g
      @test (f + g)*(f - g) == f*f + g*f - f*g - g*g
      @test f - g == -(g - f)
   end
end

@testset "Generic.NCPoly.adhoc_binary" begin
   # Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")
   for iter = 1:500
      f = rand(S, -1:10, -10:10)
      c1 = rand(R, -10:10)
      c2 = rand(R, -10:10)

      @test c1*f - c2*f == (c1 - c2)*f
      @test c1*f + c2*f == (c1 + c2)*f

      @test f*c1 - f*c2 == f*(c1 - c2)
      @test f*c1 + f*c2 == f*(c1 + c2)
   end

   # Generic tower
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")
   T, z = polynomial_ring(S, "z")
   for iter = 1:100
      f = rand(T, -1:5, 0:5, -10:10)
      c1 = rand(R, -10:10)
      c2 = rand(R, -10:10)
      d1 = rand(S, -1:5, -10:10)
      d2 = rand(S, -1:5, -10:10)

      @test c1*f - c2*f == (c1 - c2)*f
      @test c1*f + c2*f == (c1 + c2)*f
      @test d1*f - d2*f == (d1 - d2)*f
      @test d1*f + d2*f == (d1 + d2)*f

      @test f*c1 - f*c2 == f*(c1 - c2)
      @test f*c1 + f*c2 == f*(c1 + c2)
      @test f*d1 - f*d2 == f*(d1 - d2)
      @test f*d1 + f*d2 == f*(d1 + d2)
   end

   # Fix for some ambiguity
   QQx, x = QQ["x"]
   y = gen(Generic.NCPolyRing{elem_type(QQx)}(QQx, :y, false))
   @test y * x == x * y
   @test y - 2 == -(2 - y)
   @test y - QQ(2) == -(QQ(2) - y)
end

@testset "Generic.NCPoly.comparison" begin
   # Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")
   for iter = 1:500
      f = rand(S, -1:10, -10:10)
      g = deepcopy(f)
      h = rand(S, 0:10, -10:10)

      @test f == g
      @test isequal(f, g)
      @test f != g + h
   end
end

@testset "Generic.NCPoly.adhoc_comparison" begin
   # Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")
   for iter = 1:500
      f = rand(S, 0:10, -10:10)
      c1 = rand(ZZ, -10:10)
      d1 = rand(zz, -10:10)

      @test R(c1) == c1
      @test c1 == R(c1)
      @test R(d1) == d1
      @test d1 == R(d1)

      if !iszero(f)
        @test R(c1) != c1 + f
        @test c1 != R(c1) + f
        @test R(d1) != d1 + f
        @test d1 != R(d1) + f
      end
   end

   # Generic tower
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")
   T, z = polynomial_ring(S, "z")
   for iter = 1:100
      f = rand(T, 0:10, 0:5, -10:10)
      c1 = rand(ZZ, -10:10)
      d1 = rand(S, -1:5, -10:10)

      @test S(c1) == c1
      @test c1 == S(c1)
      @test S(d1) == d1
      @test d1 == S(d1)

      @test S(c1) != c1 + f
      @test c1 != S(c1) + f
      @test S(d1) != d1 + f
      @test d1 != S(d1) + f
   end
end

@testset "Generic.NCPoly.unary_ops" begin
   #  Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")
   for iter = 1:300
      f = rand(S, -1:10, -10:10)

      @test -(-f) == f
      @test iszero(f + (-f))
   end
end

@testset "Generic.NCPoly.truncation" begin
   #  Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")
   for iter = 1:300
      f = rand(S, -1:10, -10:10)
      g = rand(S, -1:10, -10:10)
      n = rand(0:20)

      @test truncate(f*g, n) == mullow(f, g, n)
   end
end

@testset "Generic.NCPoly.reverse" begin
   #  Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")
   for iter = 1:300
      f = rand(S, -1:10, -10:10)
      len = rand(length(f):12)
      frev = reverse(f, len)

      shift = 0
      for i = 1:len
         if coeff(f, i - 1) != 0
            break
         end
         shift += 1
      end

      @test length(frev) == len - shift
      @test f == reverse(frev, len)
   end
end

@testset "Generic.NCPoly.shift" begin
   # Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")
   for iter = 1:300
      f = rand(S, -1:10, -10:10)
      s = rand(0:10)
      g = s == 0 ? S() : rand(S, 0:s - 1, -10:10)

      @test shift_right(shift_left(f, s) + g, s) == f
      @test shift_left(f, s) == y^s*f
      @test length(shift_right(f, s)) == max(0, length(f) - s)
   end
end

@testset "Generic.NCPoly.powering" begin
   # Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")

   for iter = 1:10
      f = rand(S, -1:10, -10:10)
      r2 = R(1)

      for expn = 0:10
         r1 = f^expn

         @test r1 == r2

         r2 *= f
      end
   end

   f = rand(S, 0:10, -10:10)
   @test_throws DomainError f^-1
   @test_throws DomainError f^-3
   @test_throws DomainError f^identity(-1)
   @test_throws DomainError f^-rand(2:100)
end

@testset "Generic.NCPoly.exact_division" begin
   # Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")

   for iter = 1:100
      f = rand(S, -1:10, -100:100)
      g = S()
      while rank(leading_coefficient(g)) != 2
         g = rand(S, 0:10, -100:100)
      end

      @test divexact_right(f*g, g) == f
      @test divexact_left(g*f, g) == f
   end
end

@testset "Generic.NCPoly.adhoc_exact_division" begin
   # Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")

   for iter = 1:100
      f = rand(S, -1:10, -100:100)
      g = ZZ()
      while g == 0
         g = rand(ZZ, -10:10)
      end

      @test divexact_right(f*g, g) == f
      @test divexact_left(g*f, g) == f

      h = 0
      while h == 0
         h = rand(-10:10)
      end

      @test divexact_right(f*h, h) == f
      @test divexact_left(h*f, h) == f
   end

   # Generic tower
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")
   T, z = polynomial_ring(S, "z")

   for iter = 1:100
      f = rand(T, -1:10, 0:10, -100:100)
      g = S()
      while rank(leading_coefficient(g)) != 2
         g = rand(S, 0:10, -100:100)
      end

      @test divexact_right(f*g, g) == f
      @test divexact_left(g*f, g) == f

      h = ZZ()
      while h == 0
         h = rand(ZZ, -10:10)
      end

      @test divexact_right(f*h, h) == f
      @test divexact_left(h*f, h) == f
   end
end

@testset "Generic.NCPoly.evaluation" begin
   # Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")

   for iter in 1:10
      f = rand(S, -1:4, -10:10)
      g = rand(S, -1:4, -10:10)

      d = rand(ZZ, -10:10)

      @test evaluate(f + g, d) == evaluate(f, d) + evaluate(g, d)
   end

   for iter in 1:10
      f = rand(S, -1:4, -10:10)

      d = rand(ZZ, -10:10)

      @test evaluate(f^2, d) == evaluate(f, d)^2
   end
end

@testset "Generic.NCPoly.derivative" begin
   # Exact ring
   R = matrix_ring(ZZ, 2)
   S, y = polynomial_ring(R, "y")

   for iter in 1:10
      f = rand(S, -1:4, -100:100)
      g = rand(S, -1:4, -100:100)

      @test derivative(f + g) == derivative(g) + derivative(f)

      @test derivative(g*f) == derivative(g)*f + g*derivative(f)
   end
end

@testset "Generic.NCPoly.printing" begin
   M = matrix_ring(ZZ, 3)
   _, x = M['x']
   @test string(M(-1)*x) isa String
end

@testset "Generic.NCPoly.exceptions" begin
   @test_throws MethodError polynomial_ring(Char, :x)
   @test_throws Exception dense_poly_ring_type(Char)
   @test_throws ArgumentError dense_poly_type(Char)
end

@testset "Generic.NCPoly.map_coefficients" begin
   let
      Q = matrix_ring(QQ, 2) 
      Qz, z = Q[:z]
      Qx, x = QQ[:x]
      f = map_coefficients(c -> c[1, 1], 2*z; parent = Qx)
      @test f == 2*x
      @test parent(f) === Qx
   end
end
