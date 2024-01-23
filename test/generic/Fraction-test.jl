@testset "Generic.FracFieldElem.constructors" begin
   S, x = polynomial_ring(ZZ, "x")
   T = fraction_field(S)

   @test fraction_field(S, cached = true) === fraction_field(S, cached = true)
   @test fraction_field(S, cached = false) !== fraction_field(S, cached = true)

   @test elem_type(T) == Generic.FracFieldElem{elem_type(S)}
   @test elem_type(Generic.FracField{elem_type(S)}) == Generic.FracFieldElem{elem_type(S)}
   @test parent_type(Generic.FracFieldElem{elem_type(S)}) == Generic.FracField{elem_type(S)}

   @test isa(T, Generic.FracField)

   @test isa(T(3), Generic.FracFieldElem)

   @test isa(T(BigInt(7)), Generic.FracFieldElem)

   @test isa(T(x + 2), Generic.FracFieldElem)

   @test isa(T(3, 7), Generic.FracFieldElem)

   @test isa(T(x + 2, x + 1), Generic.FracFieldElem)

   @test isa(T(x + 2, 4), Generic.FracFieldElem)

   @test isa(T(3, x + 1), Generic.FracFieldElem)

   @test isa(T(T(x + 2)), Generic.FracFieldElem)

   @test isa(T(), Generic.FracFieldElem)

   @test isa((x + 3)//(x^2 + 2), Generic.FracFieldElem)

   @test isa((x + 3)//12, Generic.FracFieldElem)

   @test isa(12//(x + 2), Generic.FracFieldElem)

   @test isa((x + 1)//T(x + 2, x + 1), Generic.FracFieldElem)

   @test isa(T(x + 2, x + 1)//(x + 1), Generic.FracFieldElem)

   @test isa(T(x + 2, x + 1)//T(x, x + 2), Generic.FracFieldElem)

   TT = fraction_field(polynomial_ring(QQ, "x")[1])
   a = TT(1)
   b = T(2)

   @test a in [a, b]
   @test a in [b, a]
   @test !(a in [b])
   @test a in keys(Dict(a => 1))
   @test !(b in keys(Dict(a => 1)))
end

@testset "Generic.FracFieldElem.printing" begin
   S, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])

   @test string((x+y)//z) == "(x + y)//z"
   @test string(x//y//z) == "x//(y*z)"

   @test !occursin("\n", sprint(show, fraction_field(S)))
end

@testset "Generic.FracFieldElem.vars" begin
   S, (x,y,z) = polynomial_ring(QQ, ["x", "y", "z"])
   K = fraction_field(S)

   a = K(zero(S), one(S))
   res = vars(a)
   @test isempty(res) && eltype(res) == elem_type(S)

   @test parent(first(vars(1 // z))) == base_ring(a)
   @test vars((2y) // 3one(S)) == [y]
   @test vars((y + z) // (x + y)) == [y, z, x]
end

@testset "Generic.FracFieldElem.gens" begin
   # Univariate
   S, t = polynomial_ring(QQ, :t)
   K = fraction_field(S)

   @test gen(K) == K(t)

   # Multivariate
   S, (x,y) = polynomial_ring(QQ, ["x", "y"])
   K = fraction_field(S)

   @test ngens(S) == 2
   @test gens(K) == K.(gens(S))

   @test gens(K)[1] == gen(K,1)
   @test gens(K)[2] == gen(K,2)
   @test_throws ArgumentError gen(K,0)
   @test_throws ArgumentError gen(K,3)

   @test K[1] == gen(K,1)
   @test K[2] == gen(K,2)
   @test_throws ArgumentError K[0]
   @test_throws ArgumentError K[3]
end

@testset "Generic.FracFieldElem.rand" begin
   S, x = polynomial_ring(ZZ, "x")
   K = fraction_field(S)

   test_rand(K, -1:3, -3:3)
end

@testset "Generic.FracFieldElem.manipulation" begin
   R = fraction_field(ZZ)
   S, x = polynomial_ring(ZZ, "x")

   @test denominator((x + 1)//(-x^2 + 1)) == x - 1

   @test numerator((x + 1)//(-x^2 + 1)) == -1

   @test iszero(zero(R))

   @test isone(one(S))

   @test canonical_unit((x + 1)//(-x^2 + 1)) == -1//(x-1)

   @test is_unit((x + 1)//(-x^2 + 1))

   @test !is_zero_divisor((x + 1)//(-x^2 + 1))

   @test deepcopy((x + 1)//(-x^2 + 1)) == (x + 1)//(-x^2 + 1)

   @test characteristic(R) == 0
end

@testset "Generic.FracFieldElem.unary_ops" begin
   S, x = polynomial_ring(ZZ, "x")

   @test -((x + 1)//(-x^2 + 1)) == 1//(x - 1)
end

@testset "Generic.FracFieldElem.binary_ops" begin
   S, x = polynomial_ring(ZZ, "x")
   K = fraction_field(S)

   for iter = 1:100
      a = rand(K, -1:3, -3:3)
      b = rand(K, -1:3, -3:3)
      c = rand(K, -1:3, -3:3)

      @test a + b - b == a
      @test c*(a + b) == c*a + c*b
      @test c*(a - b) == c*(a - b)
      @test a - b == -(b - a)
   end
end

@testset "Generic.FracFieldElem.adhoc_binary" begin
   S, x = polynomial_ring(ZZ, "x")

   a = (-x + 1)//(2x^2 + 3)
   b = (x + 1)//(-x^2 + 1)

   @test a + 2 == (4*x^2-x+7)//(2*x^2+3)

   @test a - 2 == (-4*x^2-x-5)//(2*x^2+3)

   @test 3 + a == (6*x^2-x+10)//(2*x^2+3)

   @test 3 - a == (6*x^2+x+8)//(2*x^2+3)

   @test b*(x + 1) == (-x-1)//(x-1)

   @test (x + 1)*b == (-x-1)//(x-1)

   @test a + 2//3 == (4*x^2 - 3*x + 9)//(6*x^2 + 9)

   @test 2//3 + a == (4*x^2 - 3*x + 9)//(6*x^2 + 9)

   @test a - 2//3 == (-4*x^2 - 3*x - 3)//(6*x^2 + 9)

   @test 2//3 - a == (4*x^2 + 3*x + 3)//(6*x^2 + 9)

   @test a*(2//3) == (-2*x + 2)//(6*x^2 + 9)

   @test (2//3)*a == (-2*x + 2)//(6*x^2 + 9)
end

@testset "Generic.FracFieldElem.comparison" begin
   S, x = polynomial_ring(ZZ, "x")

   a = -((x + 1)//(-x^2 + 1))

   @test a == 1//(x - 1)

   @test isequal(a, 1//(x - 1))
end

@testset "Generic.FracFieldElem.adhoc_comparison" begin
   S, x = polynomial_ring(ZZ, "x")

   a = 1//(x - 1)

   @test 1//a == x - 1

   @test x - 1 == 1//a

   @test one(S) == 1

   @test 1 == one(S)

   @test  S(2)//3 == 2//3

   @test 2//3 == S(2)//3
end

@testset "Generic.FracFieldElem.powering" begin
   S, x = polynomial_ring(ZZ, "x")

   a = (x + 1)//(-x^2 + 1)

   @test a^-12 == x^12-12*x^11+66*x^10-220*x^9+495*x^8-792*x^7+924*x^6-792*x^5+495*x^4-220*x^3+66*x^2-12*x+1
end

@testset "Generic.FracFieldElem.inversion" begin
   S, x = polynomial_ring(ZZ, "x")

   a = (x + 1)//(-x^2 + 1)

   @test inv(a) == -x + 1
end

@testset "Generic.FracFieldElem.exact_division" begin
   S, x = polynomial_ring(ZZ, "x")
   K = fraction_field(S)

   for iter = 1:100
      a = K(); b = K(); c = K()

      while a == 0
         a = rand(K, -1:3, -3:3)
      end
      while b == 0
         b = rand(K, -1:3, -3:3)
      end
      while c == 0
         c = rand(K, -1:3, -3:3)
      end

      @test divexact(a*b, b) == a
      @test divexact((a + b)*c, c) == divexact(a*c, c) + divexact(b*c, c)
   end
end

@testset "Generic.FracFieldElem.adhoc_exact_division" begin
   S, x = polynomial_ring(ZZ, "x")

   a = (-x + 1)//(2x^2 + 3)
   b = (x + 1)//(-x^2 + 1)

   @test a//5 == (-x+1)//(10*x^2+15)

   @test a//(x + 1) == (-x+1)//(2*x^3+2*x^2+3*x+3)

   @test (x + 1)//b == -x^2+1

   @test 5//a == (-10*x^2-15)//(x-1)
end

@testset "Generic.FracFieldElem.divides" begin
   R, x = polynomial_ring(ZZ, "x")
   S = fraction_field(R)

   for i in 1:1000
     a = rand(S, -1:5, -10:10)
     b = rand(S, 0:5, -10:10)

     d, q = divides(a, b)

     if d
       @test b * q == a
     end
   end
end

@testset "Generic.FracFieldElem.evaluate" begin
   R, = residue_ring(ZZ, 5)
   S, x = polynomial_ring(R, "x")

   f = (x^2 + 2)//(x + 1)

   @test evaluate(f, 1) == R(4)
   @test evaluate(f, R(1)) == R(4)

   R, (x, y) = polynomial_ring(ZZ, ["x", "y"])

   f = (x^2 + y)//(y + 2)

   @test evaluate(f, [1, 2]) == ZZ(3)//ZZ(4)
   @test evaluate(f, [ZZ(1), ZZ(2)]) == ZZ(3)//ZZ(4)
end

@testset "Generic.FracFieldElem.derivative" begin
   R, x = polynomial_ring(QQ, "x")

   f = (x^2 + 1)//2x

   @test derivative(f) == (x^2 - 1)//2x^2

   R, (x, ) = polynomial_ring(QQ, ["x"])

   f = (x^2 + 1)//2x

   @test derivative(f, x) == (x^2 - 1)//2x^2

   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   f = (x^2 + x*y + y^2)//(x + y)

   @test derivative(f, x) == (x^2 + 2x*y)//(x^2 + 2x*y + y^2)

   @test derivative(f, 2) == (2x*y + y^2)//(x^2 + 2x*y + y^2)
end

@testset "Generic.FracFieldElem.square_root" begin
   R, x = polynomial_ring(QQ, "x")
   S = fraction_field(R)

   for i = 1:100
      a = rand(S, -1:5, -10:10)

      g = a^2

      @test is_square(g)

      @test sqrt(g)^2 == g

      f1, s1 = is_square_with_sqrt(g)

      @test f1 && s1^2 == g
   end

   f = (3x^2 + 2x + 1)^2//(x + 4)

   f1, s1 = is_square_with_sqrt(f)

   @test !f1

   f = (3x^2 + 2x + 1)//(x + 4)^2

   f2, s2 = is_square_with_sqrt(f)

   @test !f2   
end

@testset "Generic.FracFieldElem.gcd" begin
   S, x = polynomial_ring(ZZ, "x")

   a = (x + 1)//(-x^2 + 1) - x//(2x + 1)

   @test gcd(a, (x + 1)//(x - 1)) == 1//(2*x^2-x-1)
end

if @isdefined QQFieldElem
    @testset "Generic.FracFieldElem.remove_valuation" begin
        a = QQFieldElem(2, 3)

        @test remove(a, BigInt(2)) == (1, QQFieldElem(1, 3))
        @test valuation(a, BigInt(2)) == 1

        @test remove(a, BigInt(3)) == (-1, QQFieldElem(2, 1))
        @test valuation(a, BigInt(3)) == -1
    end
end

@testset "Generic.FracFieldElem.promotion" begin
   S, x = polynomial_ring(QQ, "x")
   F = fraction_field(S)
   T = elem_type(F)
   @test AbstractAlgebra.promote_rule(T, T) == T
end

@testset "Generic.FracFieldElem.factor" begin
   S, x = polynomial_ring(QQ, "x")
   F = fraction_field(S)
   a = (x + 1)//(x + 2)
   f = factor(a)
   @test length(f) == 0
   @test unit(f) == a
   @test_throws ArgumentError factor(zero(F))
end
