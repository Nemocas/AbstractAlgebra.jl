@testset "Julia.Floats.conformance_tests" begin
   ConformanceTests.test_Ring_interface(RDF)
   ConformanceTests.test_Ring_interface(RealField)
end

@testset "Julia.Floats.printing" begin
   R, x = polynomial_ring(RealField, "x")

   @test !occursin("+", string(2*x^2-3*x))
end

@testset "Julia.Floats.rand" begin
   R = RealField

   test(x) = @test 1.0 <= x <= 9.0
   test_rand(test, R, 1:9)
   test_rand(test, R, UnitRange(1.0, 9.0))
   test_rand(test, R, UnitRange(big(1.0), big(9.0)))
end

@testset "Julia.Floats.manipulation" begin
   R = RDF
   S = RealField

   @test is_unit(R(3))
   @test is_unit(S(3))

   @test is_negative(R(-3))
   @test !is_negative(R(0))
   @test !is_positive(R(-3))
   @test !is_positive(R(0))
end

@testset "Julia.Floats.exact_division" begin
   R = RDF
   S = RealField

   for iter = 1:1000
      a1 = rand(R, -10:10)
      a2 = rand(R, -10:10)
      b1 = rand(S, -10:10)
      b2 = rand(S, -10:10)
      c1 = rand(zz, -10:10)
      c2 = rand(ZZ, -10:10)

      @test a2 == 0 || isapprox(divexact(a1*a2, a2), a1)
      @test b2 == 0 || isapprox(divexact(b1*b2, b2), b1)

      @test c1 == 0 || isapprox(divexact(a1*c1, c1), a1)
      @test c2 == 0 || isapprox(divexact(b1*c2, c2), b1)

      @test c1 == 0 || isapprox(divexact(c1, R(1)*c1), R(1))
      @test c2 == 0 || isapprox(divexact(c2, S(1)*c2), S(1))
   end
end

@testset "Julia.Floats.divrem" begin
   R = RealField

   for iter = 1:1000
      r = rand(R, -100:100)
      s = zero(R)
      while iszero(s)
         s = rand(R, -100:100)
      end

      @test AbstractAlgebra.divrem(r,s) == (r/s,0)
   end
end

@testset "Julia.Floats.gcd" begin
   R = RDF
   S = RealField

   for iter = 1:1000
      r1 = rand(R, -100:100)
      r2 = rand(R, -100:100)
      s1 = rand(S, -100:100)
      s2 = rand(S, -100:100)

      @test (r1 == 0 && r2 == 0) || gcd(r1, r2) == 1
      @test (s1 == 0 && s2 == 0) || gcd(s1, s2) == 1
   end
end

@testset "Julia.Floats.is_square" begin
   R = RDF
   S = RealField

   z1 = R(0.0)
   z2 = S(0.0)

   @test is_square(z1)
   @test is_square(z2)

   fz1, z1root = is_square_with_sqrt(z1)
   @test fz1
   @test is_zero(z1root)

   fz2, z2root = is_square_with_sqrt(z2)
   @test fz2
   @test is_zero(z2root)

   r = R(-0.1)
   s = S(-0.1)
   @test !is_square(r)
   @test !is_square(s)

   f1, s1 = is_square_with_sqrt(r)

   @test !f1

   f2, s2 = is_square_with_sqrt(s)

   @test !f2
end
