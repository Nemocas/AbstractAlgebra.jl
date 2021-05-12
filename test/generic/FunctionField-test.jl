R1, x1 = RationalFunctionField(QQ, "x1")
U1, z1 = R1["z1"]

P1 = [(x1 + 1)*z1 + (x1 + 2), z1 + (x1 + 1)//(x1 + 2), z1^2 + 3z1 + 1,
     (x1^2 + 1)//(x1 + 1)*z1^5 + 4z1^4 + (x1 + 2)*z1^3 + x1//(x1 + 1)*z1 + 1//(x1 + 1)]

R2, x2 = RationalFunctionField(GF(23), "x2")

U2, z2 = R2["z2"]

P2 = [(x2 + 1)*z2 + (x2 + 2), z2 + (x2 + 1)//(x2 + 2), z2^2 + 3z2 + 1,
     (x2^2 + 1)//(x2 + 1)*z2^5 + 4z2^4 + (x2 + 2)*z2^3 + x2//(x2 + 1)*z2 + 1//(x2 + 1)]

@testset "Generic.FunctionField.constructors" begin 
   @test FunctionField(P1[1], "y")[1] === FunctionField(P1[1], "y", cached=true)[1]
   @test FunctionField(P1[1], "y", cached=true)[1] !== FunctionField(P1[1], "y", cached=false)[1]

   for (S, y) in  [FunctionField(P1[3], "y"), FunctionField(P2[3], "y")]
      R = base_ring(base_ring(S))
      
      @test elem_type(S) == Generic.FunctionFieldElem{elem_type(R)}
      @test elem_type(Generic.FunctionField{elem_type(R)}) == Generic.FunctionFieldElem{elem_type(R)}
      @test parent_type(Generic.FunctionFieldElem{elem_type(R)}) == Generic.FunctionField{elem_type(R)}

      @test typeof(S) <: Generic.FunctionField

      @test isa(y, Generic.FunctionFieldElem)

      a = S()

      @test isa(a, Generic.FunctionFieldElem)

      b = S(1)

      @test isa(b, Generic.FunctionFieldElem)

      c = S(ZZ(2))

      @test isa(c, Generic.FunctionFieldElem)
   end

   S1, y1 = FunctionField(P1[3], "y1")
   S2, y2 = FunctionField(P2[3], "y2")

   d = S1(QQ(2, 3))

   @test isa(d, Generic.FunctionFieldElem)

   k = S1(x1)

   @test isa(k, Generic.FunctionFieldElem)

   m = S1(y1)

   @test isa(m, Generic.FunctionFieldElem)

   @test base_ring(S1) === R1
   @test base_field(S1) === R1

   @test parent(y1) === S1

   @test isexact_type(typeof(y1)) == true

   @test characteristic(S1) == 0
   @test characteristic(S2) == 23

   @test var(S1) == :y1

   @test y1 in [y1, y2]
   @test y1 in [y2, y1]
   @test !(y1 in [y2])

   @test y1 in keys(Dict(y1 => 1))
   @test !(y2 in keys(Dict(y1 => 1)))
end

@testset "Generic.FunctionField.printing" begin
   S, y = FunctionField(P1[4], "y")

   @test string(zero(S)) == "0"
   @test string(one(S)) == "1"
   @test string(y) == "y"
   @test string(2*y^2 - y + 1) == "2*y^2 - y + 1"
   @test string((x1 + 1)//(x1 + 2)*y^2 - 1//(x1 + 2)*y + 3) ==
                           "((x1 + 1)*y^2 - y + 3*x1 + 6)//(x1 + 2)"

   S, y = FunctionField(P2[4], "y")

   @test string(zero(S)) == "0"
   @test string(one(S)) == "1"
   @test string(y) == "y"
   @test string(x2 + y + 1) == "y + x2 + 1"
end

@testset "Generic.FunctionField.rand" begin
   for f in P1
      S, y = FunctionField(f, "y")

      # TODO: test more than just the result type
      test_rand(S, 1:10, -10:10)
   end
   
   for f in P2
      S, y = FunctionField(f, "y")

      # TODO: test more than just the result type
      test_rand(S, 1:10)
   end
end

@testset "Generic.FunctionField.manipulation" begin
   for f in union(P1, P2)
      S, y = FunctionField(f, "y")

      @test iszero(zero(S))

      @test isone(one(S))

      @test isgen(gen(S))

      @test isunit(one(S))
      @test !isunit(zero(S))
      @test isunit(gen(S))

      @test degree(S) == length(modulus(S)) - 1

      @test defining_polynomial(S) == modulus(S)

      @test isa(numerator(S), PolyElem)
      @test isa(denominator(S), PolyElem)

      T = elem_type(base_ring(base_ring(S)))

      @test isa(numerator(y), Generic.Poly{<:PolyElem{T}})
      @test isa(denominator(y), PolyElem{T})
   end

   # characteristic 0
   S, y = FunctionField(P1[4], "y")

   h = x1*y^2 + (x1 + 1)//(x1 + 2)*y + 3

   @test coeff(h, 2) == x1
   @test coeff(h, 1) == (x1 + 1)//(x1 + 2)
   @test num_coeff(h, 1) == x1 + 1

   for i = 1:10
      h = rand(S, 1:10, -10:10)

      @test deepcopy(h) == h

      @test S(numerator(h), denominator(h)) == h
   end

   # characteristic p
   S, y = FunctionField(P2[4], "y")

   h = x2*y^2 + (x2 + 1)//(x2 + 2)*y + 3

   @test coeff(h, 2) == x2
   @test coeff(h, 1) == (x2 + 1)//(x2 + 2)
   @test num_coeff(h, 1) == x2 + 1

   for i = 1:10
      h = rand(S, 1:10)

      @test deepcopy(h) == h

      @test S(numerator(h), denominator(h)) == h
   end
end

@testset "Generic.FunctionField.unary_ops" begin
   for i = 1:length(P1)
      # characteristic 0
      S, y = FunctionField(P1[i], "y")

      for iters = 1:10
         f = rand(S, 1:5, -10:10)

         @test -(-f) == f
         @test iszero(f + (-f))
      end

      # characteristic p
      S, y = FunctionField(P2[i], "y")

      for iters = 1:10
         f = rand(S, 1:5)

         @test -(-f) == f
         @test iszero(f + (-f))
      end
   end
end

@testset "Generic.FunctionField.binary_ops" begin
   for i = 1:length(P1)
      # characteristic 0
      S, y = FunctionField(P1[i], "y")

      for iters = 1:10
         f = rand(S, 1:5, -10:10)
         g = rand(S, 1:5, -10:10)
         h = rand(S, 1:5, -10:10)

         @test f + g == g + f
         @test f + (g + h) == (f + g) + h
         @test f*g == g*f
         @test f*(g + h) == f*g + f*h
         @test (f - h) + (g + h) == f + g
         @test (f + g)*(f - g) == f*f - g*g
         @test f - g == -(g - f)
      end

      # characteristic p
      S, y = FunctionField(P2[i], "y")

      for iters = 1:10
         f = rand(S, 1:5)
         g = rand(S, 1:5)
         h = rand(S, 1:5)

         @test f + g == g + f
         @test f + (g + h) == (f + g) + h
         @test f*g == g*f
         @test f*(g + h) == f*g + f*h
         @test (f - h) + (g + h) == f + g
         @test (f + g)*(f - g) == f*f - g*g
         @test f - g == -(g - f)
      end
   end
end

@testset "Generic.FunctionField.adhoc_binary" begin
   for i = 1:length(P1)
      # characteristic 0
      S, y = FunctionField(P1[i], "y")

      for iters = 1:10
         f = rand(S, 1:5, -10:10)

         c1 = rand(-10:10)
         c2 = rand(-10:10)

         @test c1*f - c2*f == (c1 - c2)*f
         @test c1*f + c2*f == (c1 + c2)*f

         @test f*c1 - f*c2 == f*(c1 - c2)
         @test f*c1 + f*c2 == f*(c1 + c2)

         @test f - c1 == f + (-c1)
         @test f + c1 == f - (-c1)
         @test c1 + f == c1 + f

         c1 = rand(-10:10)//rand(1:10)
         c2 = rand(-10:10)//rand(1:10)

         @test c1*f - c2*f == (c1 - c2)*f
         @test c1*f + c2*f == (c1 + c2)*f

         @test f*c1 - f*c2 == f*(c1 - c2)
         @test f*c1 + f*c2 == f*(c1 + c2)

         @test f - c1 == f + (-c1)
         @test f + c1 == f - (-c1)
         @test c1 + f == c1 + f

         c1 = rand(ZZ, -10:10)
         c2 = rand(ZZ, -10:10)

         @test c1*f - c2*f == (c1 - c2)*f
         @test c1*f + c2*f == (c1 + c2)*f

         @test f*c1 - f*c2 == f*(c1 - c2)
         @test f*c1 + f*c2 == f*(c1 + c2)

         @test f - c1 == f + (-c1)
         @test f + c1 == f - (-c1)
         @test c1 + f == c1 + f

         c1 = rand(QQ, -10:10)
         c2 = rand(QQ, -10:10)

         @test c1*f - c2*f == (c1 - c2)*f
         @test c1*f + c2*f == (c1 + c2)*f

         @test f*c1 - f*c2 == f*(c1 - c2)
         @test f*c1 + f*c2 == f*(c1 + c2)

         @test f - c1 == f + (-c1)
         @test f + c1 == f - (-c1)
         @test c1 + f == c1 + f

         c1 = rand(R1, 1:5, -10:10)
         c2 = rand(R1, 1:5, -10:10)

         @test c1*f - c2*f == (c1 - c2)*f
         @test c1*f + c2*f == (c1 + c2)*f

         @test f*c1 - f*c2 == f*(c1 - c2)
         @test f*c1 + f*c2 == f*(c1 + c2)

         @test f - c1 == f + (-c1)
         @test f + c1 == f - (-c1)
         @test c1 + f == c1 + f
      end

      # characteristic p
      S, y = FunctionField(P2[i], "y")

      for iters = 1:10
         f = rand(S, 1:5)

         c1 = rand(-10:10)
         c2 = rand(-10:10)

         @test c1*f - c2*f == (c1 - c2)*f
         @test c1*f + c2*f == (c1 + c2)*f

         @test f*c1 - f*c2 == f*(c1 - c2)
         @test f*c1 + f*c2 == f*(c1 + c2)

         @test f - c1 == f + (-c1)
         @test f + c1 == f - (-c1)
         @test c1 + f == c1 + f

         c1 = rand(ZZ, -10:10)
         c2 = rand(ZZ, -10:10)

         @test c1*f - c2*f == (c1 - c2)*f
         @test c1*f + c2*f == (c1 + c2)*f

         @test f*c1 - f*c2 == f*(c1 - c2)
         @test f*c1 + f*c2 == f*(c1 + c2)

         @test f - c1 == f + (-c1)
         @test f + c1 == f - (-c1)
         @test c1 + f == c1 + f

         k = base_ring(base_ring(S))
         c1 = rand(k)
         c2 = rand(k)

         @test c1*f - c2*f == (c1 - c2)*f
         @test c1*f + c2*f == (c1 + c2)*f

         @test f*c1 - f*c2 == f*(c1 - c2)
         @test f*c1 + f*c2 == f*(c1 + c2)

         @test f - c1 == f + (-c1)
         @test f + c1 == f - (-c1)
         @test c1 + f == c1 + f

         c1 = rand(R2, 1:5)
         c2 = rand(R2, 1:5)

         @test c1*f - c2*f == (c1 - c2)*f
         @test c1*f + c2*f == (c1 + c2)*f

         @test f*c1 - f*c2 == f*(c1 - c2)
         @test f*c1 + f*c2 == f*(c1 + c2)

         @test f - c1 == f + (-c1)
         @test f + c1 == f - (-c1)
         @test c1 + f == c1 + f
      end
   end
end

@testset "Generic.FunctionField.comparison" begin
   for i = 1:length(P1)
      # characteristic 0
      S, y = FunctionField(P1[i], "y")

      for iters = 1:10
         f = rand(S, 1:5, -10:10)
         g = deepcopy(f)
         h = S()
         while iszero(h)
            h = rand(S, 1:5, -10:10)
         end

         @test f == g
         @test isequal(f, g)
         @test f != g + h
      end

      # characteristic p
      S, y = FunctionField(P2[i], "y")

      for iters = 1:10
         f = rand(S, 1:5)
         g = deepcopy(f)
         h = S()
         while iszero(h)
            h = rand(S, 1:4)
         end

         @test f == g
         @test isequal(f, g)
         @test f != g + h
      end
   end
end

@testset "Generic.FunctionField.adhoc_comparison" begin
   for i = 1:length(P1)
      # characteristic 0
      S, y = FunctionField(P1[i], "y")

      for iters = 1:10
         f = S()
         while iszero(f)
            f = rand(S, 1:5, -10:10)
         end

         c1 = rand(ZZ, -10:10)
 
         @test S(c1) == c1
         @test c1 == S(c1)

         @test S(c1) != c1 + f
         @test c1 != S(c1) + f

         c1 = rand(QQ, -10:10)

         @test S(c1) == c1
         @test c1 == S(c1)

         @test S(c1) != c1 + f
         @test c1 != S(c1) + f

         c1 = rand(R1, 1:5, -10:10)

         @test S(c1) == c1
         @test c1 == S(c1)

         @test S(c1) != c1 + f
         @test c1 != S(c1) + f
      end

      # characteristic p
      S, y = FunctionField(P2[i], "y")

      for iters = 1:10
         f = S()
         while iszero(f)
            f = rand(S, 1:5)
         end

         c1 = rand(ZZ, -10:10)
 
         @test S(c1) == c1
         @test c1 == S(c1)

         @test S(c1) != c1 + f
         @test c1 != S(c1) + f

         k = base_ring(base_ring(S))
         c1 = rand(k)

         @test S(c1) == c1
         @test c1 == S(c1)

         @test S(c1) != c1 + f
         @test c1 != S(c1) + f

         c1 = rand(R2, 1:5)

         @test S(c1) == c1
         @test c1 == S(c1)

         @test S(c1) != c1 + f
         @test c1 != S(c1) + f
      end
   end
end

@testset "Generic.FunctionField.powering" begin
   for i = 1:length(P1)
      # characteristic 0
      S, y = FunctionField(P1[i], "y")

      for iters = 1:2
         f = rand(S, 1:2, -1:1)
         r2 = one(S)

         for expn = 0:2
            r1 = f^expn

            @test (f == 0 && expn == 0 && r1 == 0) || r1 == r2

            r2 *= f
         end
      end

      # characteristic p
      S, y = FunctionField(P2[i], "y")
      for iters = 1:10
         f = rand(S, 1:5)
         r2 = one(S)

         for expn = 0:5
            r1 = f^expn

            @test (f == 0 && expn == 0 && r1 == 0) || r1 == r2

            r2 *= f
         end
      end
   end
end

@testset "Generic.FunctionField.inverse" begin
   for i = 1:length(P1)
      # characteristic 0
      S, y = FunctionField(P1[i], "y")

      for iters = 1:5
         f = S()
         while iszero(f)
            f = rand(S, 1:3, -10:10)
         end

         @test isone(f*inv(f))
      end

      # characteristic p
      S, y = FunctionField(P2[i], "y")

      for iters = 1:5
         f = S()
         while iszero(f)
            f = rand(S, 1:3)
         end

         @test isone(f*inv(f))
      end
   end
end

@testset "Generic.FunctionField.exact_division" begin
   for i = 1:length(P1)
      # characteristic 0
      S, y = FunctionField(P1[i], "y")

      for iters = 1:5
         f = S()
         while iszero(f)
            f = rand(S, 1:2, -1:1)
         end
         g = rand(S, 1:2, -1:1)

         @test divexact(f*g, f) == g
      end

      # characteristic p
      S, y = FunctionField(P2[i], "y")

      for iters = 1:10
         f = S()
         while iszero(f)
            f = rand(S, 1:5)
         end
         g = rand(S, 1:5)

         @test divexact(f*g, f) == g
      end
   end
end

@testset "Generic.FunctionField.adhoc_exact_division" begin
   for i = 1:length(P1)
      # characteristic 0
      S, y = FunctionField(P1[i], "y")

      for iters = 1:10
         f = 0 # Int
         while iszero(f)
            f = rand(-10:10)
         end
         g = rand(S, 1:5, -10:10)

         @test divexact(f*g, f) == g
      end

      for iters = 1:10
         f = 0//1
         while iszero(f)
            f = rand(-10:10)//rand(1:10) # Rational{Int}
         end
         g = rand(S, 1:5, -10:10)

         @test divexact(f*g, f) == g
      end

      for iters = 1:10 # BigInt
         f = ZZ()
         while iszero(f)
            f = rand(ZZ, -10:10)
         end
         g = rand(S, 1:5, -1:10)

         @test divexact(f*g, f) == g
      end

      for iters = 1:10 # BigInt
         f = QQ()
         while iszero(f)
            f = rand(QQ, -10:10) # Rational BigInt
         end
         g = rand(S, 1:5, -10:10)

         @test divexact(f*g, f) == g
      end

      for iters = 1:10 # Rat
         f = R1()
         while iszero(f)
            f = rand(R1, 1:5, -10:10) # Rational BigInt
         end
         g = rand(S, 1:5, -10:10)

         @test divexact(f*g, f) == g
      end

      # characteristic p
      S, y = FunctionField(P2[i], "y")

      for iters = 1:10
         f = 0 # Int
         while iszero(f)
            f = rand(-10:10)
         end
         g = rand(S, 1:5)

         @test divexact(f*g, f) == g
      end

      for iters = 1:10 # BigInt
         f = ZZ()
         while iszero(f)
            f = rand(ZZ, -10:10)
         end
         g = rand(S, 1:5)

         @test divexact(f*g, f) == g
      end

      for iters = 1:10 # Rat
         f = R1()
         while iszero(f)
            f = rand(R2, 1:5) # Rational BigInt
         end
         g = rand(S, 1:5)

         @test divexact(f*g, f) == g
      end
   end
end
