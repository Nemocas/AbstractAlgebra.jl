@testset "Generic.NCPoly.constructors..." begin
   R = MatrixAlgebra(ZZ, 2)
   S1 = PolynomialRing(R, "y")
   S2 = R["y"]

   for (S, y) in (S1, S2)
      @test elem_type(S) == Generic.NCPoly{elem_type(R)}
      @test elem_type(Generic.NCPolyRing{elem_type(R)}) == Generic.NCPoly{elem_type(R)}
      @test parent_type(Generic.NCPoly{elem_type(R)}) == Generic.NCPolyRing{elem_type(R)}

      @test typeof(S) <: Generic.NCPolyRing

      @test isa(y, NCPolyElem)
   end

   S3, _ = R["x"]["y"]
   @test elem_type(S3) == Generic.NCPoly{Generic.NCPoly{Generic.MatAlgElem{BigInt}}}
   @test typeof(S3) == Generic.NCPolyRing{Generic.NCPoly{Generic.MatAlgElem{BigInt}}}

   S, y = S1
   T, z = PolynomialRing(S, "z")

   @test typeof(T) <: Generic.NCPolyRing

   @test isa(z, NCPolyElem)

   f = one(R) + y^3 + z + 1

   @test isa(f, NCPolyElem)

   g = S(2)

   @test isa(g, NCPolyElem)

   h = S(rand(R, -10:10))

   @test isa(h, NCPolyElem)

   j = T(rand(R, -10:10))

   @test isa(j, NCPolyElem)

   k = S([rand(R, -10:10) for i in 1:3])

   @test isa(k, NCPolyElem)

   l = S(k)

   @test isa(l, NCPolyElem)

   m = S([1, 2, 3])

   @test isa(m, NCPolyElem)

   n = S([ZZ(1), ZZ(2), ZZ(3)])

   @test isa(n, NCPolyElem)
end

@testset "Generic.NCPoly.manipulation..." begin
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")

   @test iszero(zero(S))

   @test isone(one(S))

   @test isgen(gen(S))

   @test isunit(one(S))

   f = 2*y + 1

   @test lead(f) == 2

   @test trail(2y^2 + 3y) == 3

   @test degree(f) == 1

   @test one(f) == one(S)
   @test zero(f) == zero(S)

   h = y^2 + 3*y + 3

   @test coeff(h, 1) == 3

   @test length(h) == 3

   @test canonical_unit(-y + 1) == -1

   @test deepcopy(h) == h

   @test isterm(rand(R, -10:10)*y^2)

   @test !ismonomial(2*y^2)

   @test ismonomial(y^2)

   @test !ismonomial(y^2 + y + 1)
end

@testset "Generic.NCPoly.rand..." begin
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")

   f = rand(S, 0:10, -10:10)
   @test f isa Generic.NCPoly
   f = rand(rng, S, 0:10, -10:10)
   @test f isa Generic.NCPoly
end

@testset "Generic.NCPoly.binary_ops..." begin
   #  Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")

   for iter = 1:100
      f = rand(S, 0:10, -10:10)
      g = rand(S, 0:10, -10:10)
      h = rand(S, 0:10, -10:10)
      @test f + g == g + f
      @test f + (g + h) == (f + g) + h
      @test f*(g + h) == f*g + f*h
      @test (f - h) + (g + h) == f + g
      @test (f + g)*(f - g) == f*f + g*f - f*g - g*g
      @test f - g == -(g - f)
   end
end

@testset "Generic.NCPoly.adhoc_binary..." begin
   # Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")
   for iter = 1:500
      f = rand(S, 0:10, -10:10)
      c1 = rand(R, -10:10)
      c2 = rand(R, -10:10)

      @test c1*f - c2*f == (c1 - c2)*f
      @test c1*f + c2*f == (c1 + c2)*f

      @test f*c1 - f*c2 == f*(c1 - c2)
      @test f*c1 + f*c2 == f*(c1 + c2)
   end

   # Generic tower
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")
   for iter = 1:100
      f = rand(T, 0:5, 0:5, -10:10)
      c1 = rand(R, -10:10)
      c2 = rand(R, -10:10)
      d1 = rand(S, 0:5, -10:10)
      d2 = rand(S, 0:5, -10:10)

      @test c1*f - c2*f == (c1 - c2)*f
      @test c1*f + c2*f == (c1 + c2)*f
      @test d1*f - d2*f == (d1 - d2)*f
      @test d1*f + d2*f == (d1 + d2)*f

      @test f*c1 - f*c2 == f*(c1 - c2)
      @test f*c1 + f*c2 == f*(c1 + c2)
      @test f*d1 - f*d2 == f*(d1 - d2)
      @test f*d1 + f*d2 == f*(d1 + d2)
   end
end

@testset "Generic.NCPoly.comparison..." begin
   # Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")
   for iter = 1:500
      f = rand(S, 0:10, -10:10)
      g = deepcopy(f)
      h = S()
      while iszero(h)
         h = rand(S, 0:10, -10:10)
      end

      @test f == g
      @test isequal(f, g)
      @test f != g + h
   end
end

@testset "Generic.NCPoly.adhoc_comparison..." begin
   # Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")
   for iter = 1:500
      f = S()
      while iszero(f)
         f = rand(S, 0:10, -10:10)
      end
      c1 = rand(ZZ, -10:10)
      d1 = rand(zz, -10:10)

      @test R(c1) == c1
      @test c1 == R(c1)
      @test R(d1) == d1
      @test d1 == R(d1)

      @test R(c1) != c1 + f
      @test c1 != R(c1) + f
      @test R(d1) != d1 + f
      @test d1 != R(d1) + f
   end

   # Generic tower
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")
   for iter = 1:100
      f = T()
      while iszero(f)
         f = rand(T, 0:10, 0:5, -10:10)
      end
      c1 = rand(ZZ, -10:10)
      d1 = rand(S, 0:5, -10:10)

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

@testset "Generic.NCPoly.unary_ops..." begin
   #  Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")
   for iter = 1:300
      f = rand(S, 0:10, -10:10)

      @test -(-f) == f
      @test iszero(f + (-f))
   end
end

@testset "Generic.NCPoly.truncation..." begin
   #  Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")
   for iter = 1:300
      f = rand(S, 0:10, -10:10)
      g = rand(S, 0:10, -10:10)
      n = rand(0:20)

      @test truncate(f*g, n) == mullow(f, g, n)
   end
end

@testset "Generic.NCPoly.reverse..." begin
   #  Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")
   for iter = 1:300
      f = rand(S, 0:10, -10:10)
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

@testset "Generic.NCPoly.shift..." begin
   # Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")
   for iter = 1:300
      f = rand(S, 0:10, -10:10)
      s = rand(0:10)
      g = s == 0 ? S() : rand(S, 0:s - 1, -10:10)

      @test shift_right(shift_left(f, s) + g, s) == f
      @test shift_left(f, s) == y^s*f
      @test length(shift_right(f, s)) == max(0, length(f) - s)
   end
end

@testset "Generic.NCPoly.powering..." begin
   # Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")

   for iter = 1:10
      f = rand(S, 0:10, -10:10)
      r2 = R(1)

      for expn = 0:10
         r1 = f^expn

         @test (f == 0 && expn == 0 && r1 == 0) || r1 == r2

         r2 *= f
      end
   end

   f = rand(S, 0:10, -10:10)
   @test_throws DomainError f^-1
   @test_throws DomainError f^-3   
   @test_throws DomainError f^identity(-1)
   @test_throws DomainError f^-rand(2:100)
end

@testset "Generic.NCPoly.exact_division..." begin
   # Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")

   for iter = 1:100
      f = rand(S, 0:10, -100:100)
      g = S()
      while rank(lead(g)) != 2
         g = rand(S, 0:10, -100:100)
      end

      @test divexact_right(f*g, g) == f
      @test divexact_left(g*f, g) == f
   end
end

@testset "Generic.NCPoly.adhoc_exact_division..." begin
   # Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")

   for iter = 1:100
      f = rand(S, 0:10, -100:100)
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
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")

   for iter = 1:100
      f = rand(T, 0:10, 0:10, -100:100)
      g = S()
      while rank(lead(g)) != 2
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

@testset "Generic.NCPoly.evaluation..." begin
   # Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")

   for iter in 1:10
      f = rand(S, 0:4, -10:10)
      g = rand(S, 0:4, -10:10)

      d = rand(ZZ, -10:10)

      @test evaluate(f + g, d) == evaluate(f, d) + evaluate(g, d)
   end

   for iter in 1:10
      f = rand(S, 0:4, -10:10)

      d = rand(ZZ, -10:10)

      @test evaluate(f^2, d) == evaluate(f, d)^2
   end
end

@testset "Generic.NCPoly.derivative..." begin
   # Exact ring
   R = MatrixAlgebra(ZZ, 2)
   S, y = PolynomialRing(R, "y")

   for iter in 1:10
      f = rand(S, 0:4, -100:100)
      g = rand(S, 0:4, -100:100)

      @test derivative(f + g) == derivative(g) + derivative(f)

      @test derivative(g*f) == derivative(g)*f + g*derivative(f)
   end
end
