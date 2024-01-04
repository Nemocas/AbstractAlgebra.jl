@testset "Generic.RationalFunctionField.constructors" begin
   # Univariate

   T, x = rational_function_field(QQ, "x")

   @test rational_function_field(QQ, "x", cached = true)[1] === rational_function_field(QQ, "x", cached = true)[1]
   @test rational_function_field(QQ, "x", cached = false)[1] !== rational_function_field(QQ, "x", cached = false)[1]

   @test elem_type(T) == Generic.RationalFunctionFieldElem{Rational{BigInt}, Generic.Poly{Rational{BigInt}}}
   @test elem_type(Generic.RationalFunctionField{Rational{BigInt}, Generic.Poly{Rational{BigInt}}}) == Generic.RationalFunctionFieldElem{Rational{BigInt}, Generic.Poly{Rational{BigInt}}}
   @test parent_type(Generic.RationalFunctionFieldElem{Rational{BigInt}, Generic.Poly{Rational{BigInt}}}) == Generic.RationalFunctionField{Rational{BigInt}, Generic.Poly{Rational{BigInt}}}

   @test isa(T, Generic.RationalFunctionField)

   @test isa(T(3), Generic.RationalFunctionFieldElem)

   @test isa(T(BigInt(7)), Generic.RationalFunctionFieldElem)

   @test isa(T(x + 2), Generic.RationalFunctionFieldElem)

   @test isa(T(numerator(x^2 + 2x + 1, false), numerator(x + 1, false)), Generic.RationalFunctionFieldElem)

   @test isa(T(), Generic.RationalFunctionFieldElem)

   @test isa((x + 3)//(x^2 + 2), Generic.RationalFunctionFieldElem)

   @test isa((x + 3)//12, Generic.RationalFunctionFieldElem)

   @test isa(12//(x + 2), Generic.RationalFunctionFieldElem)

   @test isa(((x + 1)//(x + 2))//((x + 3)//(x + 4)), Generic.RationalFunctionFieldElem)

   TT, x = rational_function_field(QQ, "x")
   a = TT(1)
   b = T(2)

   @test a in [a, b]
   @test a in [b, a]
   @test !(a in [b])
   @test a in keys(Dict(a => 1))
   @test !(b in keys(Dict(a => 1)))

   # Multivariate

   T, (x, y) = rational_function_field(QQ, ["x", "y"])

   @test rational_function_field(QQ, ["x", "y"], cached = true)[1] === rational_function_field(QQ, ["x", "y"], cached = true)[1]
   @test rational_function_field(QQ, ["x", "y"], cached = false)[1] !== rational_function_field(QQ, ["x", "y"], cached = false)[1]

   @test elem_type(T) == Generic.RationalFunctionFieldElem{Rational{BigInt}, Generic.MPoly{Rational{BigInt}}}

   @test isa(T, Generic.RationalFunctionField)

   @test isa(T(3), Generic.RationalFunctionFieldElem)

   @test isa(T(BigInt(7)), Generic.RationalFunctionFieldElem)

   @test isa(T(x + 2), Generic.RationalFunctionFieldElem)

   @test isa(T(x + y), Generic.RationalFunctionFieldElem)

   @test isa(T(numerator(x^2 + 2y + 1, false), numerator(y + 1, false)), Generic.RationalFunctionFieldElem)

   @test isa(T(), Generic.RationalFunctionFieldElem)

   @test isa((x + 3)//(y^2 + 2), Generic.RationalFunctionFieldElem)

   @test isa((x + 3)//12, Generic.RationalFunctionFieldElem)

   @test isa(12//(y + 2), Generic.RationalFunctionFieldElem)

   @test isa(((x + 1)//(y + 2))//((y + 3)//(x + 4)), Generic.RationalFunctionFieldElem)
end

@testset "Generic.RationalFunctionField.printing" begin
   # Univariate
   S, x = rational_function_field(QQ, "x")

   @test string((x+1)//x) == "(x + 1)//x"
   @test string(x//x//x) == "1//x"

   # Multivariate
   S, (x, y) = rational_function_field(QQ, ["x", "y"])

   @test string((x+1)//y) == "(x + 1)//y"
   @test string(y//x//y) == "1//x"
end

@testset "Generic.RationalFunctionField.gens" begin
   # Univariate
   K, x = rational_function_field(QQ, "x")

   @test gen(K) == x

   # Multivariate
   K, (x, y) = rational_function_field(QQ, ["x", "y"])

   @test ngens(K) == 2
   @test gens(K) == [x, y]

   @test gen(K,1) == x
   @test gen(K,2) == y
   @test_throws ArgumentError gen(K,0)
   @test_throws ArgumentError gen(K,3)

   @test K[1] == x
   @test K[2] == y
   @test_throws ArgumentError K[0]
   @test_throws ArgumentError K[3]
end

@testset "Generic.RationalFunctionField.rand" begin
   # Univariate
   K, x = rational_function_field(QQ, "x")

   test_rand(K, 0:3, -3:3)

   # Multivariate
   K, (x, y) = rational_function_field(QQ, ["x", "y"])

   test_rand(K, 0:3, 0:3, -3:3)
end

@testset "Generic.RationalFunctionField.manipulation" begin
   # Univariate
   R, x = rational_function_field(QQ, "x")

   @test denominator((x + 1)//(-x^2 + 1)) == numerator(x - 1)

   @test numerator((x + 1)//(-x^2 + 1)) == -1

   @test iszero(zero(R))

   @test isone(one(R))

   @test canonical_unit((x + 1)//(-x^2 + 1)) == -1//(x-1)

   @test is_unit((x + 1)//(-x^2 + 1))

   @test !is_zero_divisor((x + 1)//(-x^2 + 1))

   @test deepcopy((x + 1)//(-x^2 + 1)) == (x + 1)//(-x^2 + 1)

   @test characteristic(R) == 0

   # Multivariate
   R, (x, y) = rational_function_field(QQ, ["x", "y"])

   @test denominator((y + 1)//(x*y + x - y - 1)) == numerator(x - 1)

   @test numerator((y + 1)//(x*y + x - y - 1)) == 1

   @test iszero(zero(R))

   @test isone(one(R))

   @test canonical_unit((y + 1)//(x*y + x - y - 1)) == 1//(x-1)

   @test is_unit((x + 1)//(-y^2 + 1))

   @test deepcopy((x + 1)//(-y^2 + 1)) == (x + 1)//(-y^2 + 1)

   @test characteristic(R) == 0
end

@testset "Generic.RationalFunctionField.unary_ops" begin
   # Univariate
   R, x = rational_function_field(QQ, "x")

   @test -((x + 1)//(-x^2 + 1)) == 1//(x - 1)

   R, (x, y) = rational_function_field(QQ, ["x", "y"])

   @test -((y + 1)//(x*y + x - y - 1)) == -1//(x - 1)
end

@testset "Generic.RationalFunctionField.binary_ops" begin
   # Univariate
   K, x = rational_function_field(QQ, "x")

   for iter = 1:100
      a = rand(K, 0:3, -3:3)
      b = rand(K, 0:3, -3:3)
      c = rand(K, 0:3, -3:3)

      @test a + b - b == a
      @test c*(a + b) == c*a + c*b
      @test c*(a - b) == c*(a - b)
      @test a - b == -(b - a)
   end

   # Multivariate
   K, (x, y) = rational_function_field(QQ, ["x", "y"])

   for iter = 1:100
      a = rand(K, 0:3, 0:3, -3:3)
      b = rand(K, 0:3, 0:3, -3:3)
      c = rand(K, 0:3, 0:3, -3:3)

      @test a + b - b == a
      @test c*(a + b) == c*a + c*b
      @test c*(a - b) == c*(a - b)
      @test a - b == -(b - a)
   end
end

@testset "Generic.RationalFunctionField.adhoc_binary" begin
   # Univariate
   R, x = rational_function_field(QQ, "x")

   a = (-x + 1)//(2x^2 + 3)
   b = (x + 1)//(-x^2 + 1)

   @test a + 2 == (4*x^2-x+7)//(2*x^2+3)

   @test a - 2 == (-4*x^2-x-5)//(2*x^2+3)

   @test 3 + a == (6*x^2-x+10)//(2*x^2+3)

   @test 3 - a == (6*x^2+x+8)//(2*x^2+3)

   @test b*(x + 1) == (-x-1)//(x-1)

   @test (x + 1)*b == (-x-1)//(x-1)

   @test denominator(a) * a == a * denominator(a)
   @test denominator(a) + a == a + denominator(a)

   # Multivariate
   R, (x, y) = rational_function_field(QQ, ["x", "y"])

   a = (-x + 1)//(2x^2 + 3)
   b = (y + 1)//(-x*y - x + y + 1)

   @test a + 2 == (4*x^2-x+7)//(2*x^2+3)

   @test a - 2 == (-4*x^2-x-5)//(2*x^2+3)

   @test 3 + a == (6*x^2-x+10)//(2*x^2+3)

   @test 3 - a == (6*x^2+x+8)//(2*x^2+3)

   @test b*(x + 1) == (-x-1)//(x-1)

   @test (x + 1)*b == (-x-1)//(x-1)

   @test denominator(a) * a == a * denominator(a)
   @test denominator(a) + a == a + denominator(a)
end

@testset "Generic.RationalFunctionField.comparison" begin
   # Univariate
   R, x = rational_function_field(QQ, "x")

   a = -((x + 1)//(-x^2 + 1))

   @test a == 1//(x - 1)

   @test isequal(a, 1//(x - 1))

   # Multivariate
   R, (x, y) = rational_function_field(QQ, ["x", "y"])

   a = -((y + 1)//(-x*y - x + y + 1))

   @test a == 1//(x - 1)

   @test isequal(a, 1//(x - 1))
end

@testset "Generic.RationalFunctionField.adhoc_comparison" begin
   # Univariate
   R, x = rational_function_field(QQ, "x")

   a = 1//(x - 1)

   @test 1//a == x - 1

   @test x - 1 == 1//a

   @test one(R) == 1

   @test 1 == one(R)

   @test R(x + 1) == x + 1

   # Univariate
   R, (x, y) = rational_function_field(QQ, ["x", "y"])

   a = 1//(x - 1)

   @test 1//a == x - 1

   @test x - 1 == 1//a

   @test one(R) == 1

   @test 1 == one(R)

   @test R(x + 1) == x + 1
end

@testset "Generic.RationalFunctionField.powering" begin
   # Univariate
   R, x = rational_function_field(QQ, "x")

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

   # Multivariate
   R, (x, y) = rational_function_field(QQ, ["x", "y"])

   a = (x + 1)//(-y^2 + 1)

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
   # Univariate
   R, x = rational_function_field(QQ, "x")

   a = (x + 1)//(-x^2 + 1)

   @test inv(a) == -x + 1
   @test_throws NotInvertibleError inv(R())

   # Multivariate
   R, (x, y) = rational_function_field(QQ, ["x", "y"])

   a = (y + 1)//(-x*y - x + y + 1)

   @test inv(a) == -x + 1
   @test_throws NotInvertibleError inv(R())
end

@testset "Generic.RationalFunctionField.exact_division" begin
   # Univariate
   K, x = rational_function_field(QQ, "x")

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

   # Multivariate
   K, (x, y) = rational_function_field(QQ, ["x", "y"])

   for iter = 1:100
      a = K()
      b = K()
      c = K()

      while a == 0
         a = rand(K, 0:3, 0:3, -3:3)
      end
      while b == 0
         b = rand(K, 0:3, 0:3, -3:3)
      end
      while c == 0
         c = rand(K, 0:3, 0:3, -3:3)
      end

      @test divexact(a*b, b) == a
      @test divexact((a + b)*c, c) == divexact(a*c, c) + divexact(b*c, c)
   end
end

@testset "Generic.RationalFunctionField.adhoc_exact_division" begin
   # Univariate
   R, x = rational_function_field(QQ, "x")

   a = (-x + 1)//(2x^2 + 3)
   b = (x + 1)//(-x^2 + 1)

   @test a//5 == (-x+1)//(10*x^2+15)

   @test a//(x + 1) == (-x+1)//(2*x^3+2*x^2+3*x+3)

   @test (x + 1)//b == -x^2+1

   @test 5//a == (-10*x^2-15)//(x-1)

   # Multivariate
   R, (x, y) = rational_function_field(QQ, ["x", "y"])

   a = (-x + 1)//(2x^2 + 3)
   b = (y + 1)//(-x*y - x + y + 1)

   @test a//5 == (-x+1)//(10*x^2+15)

   @test a//(x + 1) == (-x+1)//(2*x^3+2*x^2+3*x+3)

   @test (x + 1)//b == -x^2+1

   @test 5//a == (-10*x^2-15)//(x-1)
end

@testset "Generic.RationalFunctionField.divides" begin
   # Univariate
   R, x = rational_function_field(QQ, "x")

   for i in 1:100
     a = rand(R, 1:5, -10:10)
     b = rand(R, 1:5, -10:10)

     d, q = divides(a, b)

     if d
       @test b * q == a
     end
   end

   # Univariate
   R, (x, y) = rational_function_field(QQ, ["x", "y"])

   for i in 1:100
     a = rand(R, 1:5, 0:3, -10:10)
     b = rand(R, 1:5, 0:3, -10:10)

     d, q = divides(a, b)

     if d
       @test b * q == a
     end
   end
end

@testset "Generic.RationalFunctionField.evaluate" begin
   # Univariate
   R, x = polynomial_ring(QQ, "x")

   f = (x^2 + 2)//(x + 1)

   @test evaluate(f, 1) == QQ(3, 2)
   @test evaluate(f, QQ(2)) == 2

   # Multivariate
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   f = (x^2 + 2)//(y + 1)

   @test evaluate(f, [1, 2]) == 1
   @test evaluate(f, [QQ(2), QQ(1)]) == 3
end

@testset "Generic.RationalFunctionField.derivative" begin
   # Univariate
   R, x = rational_function_field(QQ, "x")

   f = (x^2 + 1)//2x

   @test derivative(f) == (x^2 - 1)//2x^2

   # Multivariate
   R, (x, y) = rational_function_field(QQ, ["x", "y"])

   f = (x^2 + 1)//2y

   @test derivative(f, numerator(x)) == x//y
end

@testset "Generic.RationalFunctionField.square_root" begin
   # Univariate
   R, x = rational_function_field(QQ, "x")

   for i = 1:100
      a = rand(R, 1:5, -10:10)

      @test is_square(a^2)

      @test sqrt(a^2)^2 == a^2
   end

   # Multivariate
   R, (x, y) = rational_function_field(QQ, ["x", "y"])

   for i = 1:100
      a = rand(R, 1:5, 0:3, -10:10)

      @test is_square(a^2)

      @test sqrt(a^2)^2 == a^2
   end
end

@testset "Generic.RationalFunctionField.gcd" begin
   # Univariate
   R, x = rational_function_field(QQ, "x")

   a = (x + 1)//(-x^2 + 1) - x//(2x + 1)

   @test gcd(a, (x + 1)//(x - 1)) == 1//(x^2-1//2*x-1//2)

   # Multivariate
   R, (x, y) = rational_function_field(QQ, ["x", "y"])

   a = (x + 1)//(y + 1) - y//(2x + 1)

   @test gcd(a, (x + 1)//(y - 1)) == 1//(x*y^2 - x + 1//2*y^2 - 1//2)
end

@testset "Generic.RationalFunctionField.promotion" begin
   # Univariate
   K, x = rational_function_field(QQ, "x")
   
   T = elem_type(K)
   
   @test AbstractAlgebra.promote_rule(T, T) == T

   # Multivariate
   K, (x, y) = rational_function_field(QQ, ["x", "y"])

   T = elem_type(K)

   @test AbstractAlgebra.promote_rule(T, T) == T
end
