@testset "Generic.RationalFunctionField.constructors" begin
   T, x = RationalFunctionField(QQ, "x")

   @test RationalFunctionField(QQ, "x", cached = true)[1] === RationalFunctionField(QQ, "x", cached = true)[1]
   @test RationalFunctionField(QQ, "x", cached = false)[1] !== RationalFunctionField(QQ, "x", cached = false)[1]

   @test elem_type(T) == Generic.Rat{Rational{BigInt}}
   @test elem_type(Generic.RationalFunctionField{Rational{BigInt}}) == Generic.Rat{Rational{BigInt}}
   @test parent_type(Generic.Rat{Rational{BigInt}}) == Generic.RationalFunctionField{Rational{BigInt}}

   @test isa(T, Generic.RationalFunctionField)

   @test isa(T(3), Generic.Rat)

   @test isa(T(BigInt(7)), Generic.Rat)

   @test isa(T(x + 2), Generic.Rat)

   @test isa(T(x^2 + 2x + 1, x + 1), Generic.Rat)

   @test isa(T(), Generic.Rat)

   @test isa((x + 3)//(x^2 + 2), Generic.Rat)

   @test isa((x + 3)//12, Generic.Rat)

   @test isa(12//(x + 2), Generic.Rat)

   @test isa(((x + 1)//(x + 2))//((x + 3)//(x + 4)), Generic.Rat)

   TT, x = RationalFunctionField(QQ, "x")
   a = TT(1)
   b = T(2)

   @test a in [a, b]
   @test a in [b, a]
   @test !(a in [b])
   @test a in keys(Dict(a => 1))
   @test !(b in keys(Dict(a => 1)))
end

@testset "Generic.RationalFunctionField.printing" begin
   S, x = RationalFunctionField(QQ, "x")

   @test string((x+1)//x) == "(x + 1)//x"
   @test string(x//x//x) == "1//x"
end


@testset "Generic.RationalFunctionField.rand" begin
   K, x = RationalFunctionField(QQ, "x")

   test_rand(K, 0:3, -3:3)
end

@testset "Generic.RationalFunctionField.manipulation" begin
   R, x = RationalFunctionField(QQ, "x")

   @test denominator((x + 1)//(-x^2 + 1)) == numerator(x - 1)

   @test numerator((x + 1)//(-x^2 + 1)) == -1

   @test iszero(zero(R))

   @test isone(one(R))

   @test canonical_unit((x + 1)//(-x^2 + 1)) == -1//(x-1)

   @test isunit((x + 1)//(-x^2 + 1))

   @test deepcopy((x + 1)//(-x^2 + 1)) == (x + 1)//(-x^2 + 1)

   @test characteristic(R) == 0
end

@testset "Generic.RationalFunctionField.unary_ops" begin
   R, x = RationalFunctionField(QQ, "x")

   @test -((x + 1)//(-x^2 + 1)) == 1//(x - 1)
end

@testset "Generic.RationalFunctionField.binary_ops" begin
   K, x = RationalFunctionField(QQ, "x")

   for iter = 1:100
      a = rand(K, 0:3, -3:3)
      b = rand(K, 0:3, -3:3)
      c = rand(K, 0:3, -3:3)

      @test a + b - b == a
      @test c*(a + b) == c*a + c*b
      @test c*(a - b) == c*(a - b)
      @test a - b == -(b - a)
   end
end

@testset "Generic.RationalFunctionField.adhoc_binary" begin
   R, x = RationalFunctionField(QQ, "x")

   a = (-x + 1)//(2x^2 + 3)
   b = (x + 1)//(-x^2 + 1)

   @test a + 2 == (4*x^2-x+7)//(2*x^2+3)

   @test a - 2 == (-4*x^2-x-5)//(2*x^2+3)

   @test 3 + a == (6*x^2-x+10)//(2*x^2+3)

   @test 3 - a == (6*x^2+x+8)//(2*x^2+3)

   @test b*(x + 1) == (-x-1)//(x-1)

   @test (x + 1)*b == (-x-1)//(x-1)
end

@testset "Generic.RationalFunctionField.comparison" begin
   R, x = RationalFunctionField(QQ, "x")

   a = -((x + 1)//(-x^2 + 1))

   @test a == 1//(x - 1)

   @test isequal(a, 1//(x - 1))
end

@testset "Generic.RationalFunctionField.adhoc_comparison" begin
   R, x = RationalFunctionField(QQ, "x")

   a = 1//(x - 1)

   @test 1//a == x - 1

   @test x - 1 == 1//a

   @test one(R) == 1

   @test 1 == one(R)

   R(x + 1) == x + 1
end

@testset "Generic.RationalFunctionField.powering" begin
   R, x = RationalFunctionField(QQ, "x")

   a = (x + 1)//(-x^2 + 1)

   b = one(R)

   for iter = 0:10
      @test b == a^iter
      b *= a
   end

   b = one(R)

   for iter = 0:10
      @test 1//b == a^-iter
      b *= a
   end
end

@testset "Generic.RationalFunctionField.inversion" begin
   R, x = RationalFunctionField(QQ, "x")

   a = (x + 1)//(-x^2 + 1)

   @test inv(a) == -x + 1
   @test_throws DivideError inv(R())
end

@testset "Generic.RationalFunctionField.exact_division" begin
   K, x = RationalFunctionField(QQ, "x")

   for iter = 1:100
      a = K()
      b = K()
      c = K()

      while a == 0
         a = rand(K, 0:3, -3:3)
      end
      while b == 0
         b = rand(K, 0:3, -3:3)
      end
      while c == 0
         c = rand(K, 0:3, -3:3)
      end

      @test divexact(a*b, b) == a
      @test divexact((a + b)*c, c) == divexact(a*c, c) + divexact(b*c, c)
   end
end

@testset "Generic.RationalFunctionField.adhoc_exact_division" begin
   R, x = RationalFunctionField(QQ, "x")

   a = (-x + 1)//(2x^2 + 3)
   b = (x + 1)//(-x^2 + 1)

   @test a//5 == (-x+1)//(10*x^2+15)

   @test a//(x + 1) == (-x+1)//(2*x^3+2*x^2+3*x+3)

   @test (x + 1)//b == -x^2+1

   @test 5//a == (-10*x^2-15)//(x-1)
end

@testset "Generic.RationalFunctionField.divides" begin
   R, x = RationalFunctionField(QQ, "x")

   for i in 1:100
     a = rand(R, 1:5, -10:10)
     b = rand(R, 1:5, -10:10)

     d, q = divides(a, b)

     if d
       @test b * q == a
     end
   end
end

@testset "Generic.RationalFunctionField.evaluate" begin
   R, x = PolynomialRing(QQ, "x")

   f = (x^2 + 2)//(x + 1)

   @test evaluate(f, 1) == QQ(3, 2)
   @test evaluate(f, QQ(2)) == 2
end

@testset "Generic.RationalFunctionField.derivative" begin
   R, x = RationalFunctionField(QQ, "x")

   f = (x^2 + 1)//2x

   @test derivative(f) == (x^2 - 1)//2x^2
end

@testset "Generic.RationalFunctionField.square_root" begin
   R, x = RationalFunctionField(QQ, "x")

   for i = 1:100
      a = rand(R, 1:5, -10:10)

      @test issquare(a^2)

      @test sqrt(a^2)^2 == a^2
   end
end

@testset "Generic.RationalFunctionField.gcd" begin
   R, x = RationalFunctionField(QQ, "x")

   a = (x + 1)//(-x^2 + 1) - x//(2x + 1)

   @test gcd(a, (x + 1)//(x - 1)) == 1//(x^2-1//2*x-1//2)
end

@testset "Generic.RationalFunctionField.promotion" begin
   K, x = RationalFunctionField(QQ, "x")
   T = elem_type(K)
   @test AbstractAlgebra.promote_rule(T, T) == T
end
