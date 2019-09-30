@testset "Julia.GFElem.constructors..." begin
   R = GF(13)
   S = GF(BigInt(13))

   @test elem_type(R) == AbstractAlgebra.GFElem{Int}
   @test elem_type(S) == AbstractAlgebra.GFElem{BigInt}
   @test elem_type(AbstractAlgebra.GFField{Int}) == AbstractAlgebra.GFElem{Int}
   @test elem_type(AbstractAlgebra.GFField{BigInt}) == AbstractAlgebra.GFElem{BigInt}
   @test parent_type(AbstractAlgebra.GFElem{Int}) == AbstractAlgebra.GFField{Int}
   @test parent_type(AbstractAlgebra.GFElem{BigInt}) == AbstractAlgebra.GFField{BigInt}

   @test isa(R, AbstractAlgebra.GFField)
   @test isa(S, AbstractAlgebra.GFField)

   @test isa(R(), AbstractAlgebra.GFElem)
   @test isa(S(), AbstractAlgebra.GFElem)

   @test isa(R(11), AbstractAlgebra.GFElem)
   @test isa(S(BigInt(11)), AbstractAlgebra.GFElem)
   @test isa(S(11), AbstractAlgebra.GFElem)
   @test isa(S(BigInt(11)), AbstractAlgebra.GFElem)

   a = R(11)
   b = S(11)

   @test isa(R(a), AbstractAlgebra.GFElem)
   @test isa(S(b), AbstractAlgebra.GFElem)
end

@testset "Julia.GFElem.printing..." begin
   R = GF(13)
   S = GF(BigInt(13))

   @test string(R(3)) == "3"
   @test string(R()) == "0"
   @test string(S(3)) == "3"
   @test string(S()) == "0"
end

@testset "Julia.GFElem.manipulation..." begin
   R = GF(13)
   S = GF(BigInt(13))

   @test iszero(zero(R))
   @test iszero(zero(S))

   @test isone(one(R))
   @test isone(one(S))

   @test characteristic(R) == 13
   @test characteristic(S) == 13

   @test order(R) == 13
   @test order(S) == 13

   @test degree(R) == 1
   @test degree(S) == 1

   @test !isunit(R())
   @test !isunit(S())
   @test isunit(R(3))
   @test isunit(S(3))

   @test deepcopy(R(3)) == R(3)
   @test deepcopy(S(3)) == S(3)

   R1 = GF(13)
   S1 = GF(BigInt(13))

   @test R === R1
   @test S === S1
end

@testset "Julia.GFElem.rand..." begin
   R = GF(13)
   @test rand(R) isa AbstractAlgebra.GFElem
   @test rand(rng, R) isa AbstractAlgebra.GFElem
   @test rand(R, 2, 3) isa Matrix{<:AbstractAlgebra.GFElem}
   @test rand(rng, R, 2, 3) isa Matrix{<:AbstractAlgebra.GFElem}
end

@testset "Julia.GFElem.unary_ops..." begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:1000
      a = rand(R)
      b = rand(S)

      @test a == -(-a)
      @test b == -(-b)
   end
end

@testset "Julia.GFElem.binary_ops..." begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:1000
      a1 = rand(R)
      a2 = rand(R)
      a3 = rand(R)

      @test a1 + a2 == a2 + a1
      @test a1 - a2 == -(a2 - a1)
      @test a1 + R() == a1
      @test a1 + (a2 + a3) == (a1 + a2) + a3
      @test a1*(a2 + a3) == a1*a2 + a1*a3
      @test a1*a2 == a2*a1
      @test a1*R(1) == a1
      @test R(1)*a1 == a1

      b1 = rand(S)
      b2 = rand(S)
      b3 = rand(S)

      @test b1 + b2 == b2 + b1
      @test b1 - b2 == -(b2 - b1)
      @test b1 + S() == b1
      @test b1 + (b2 + b3) == (b1 + b2) + b3
      @test b1*(b2 + b3) == b1*b2 + b1*b3
      @test b1*b2 == b2*b1
      @test b1*S(1) == b1
      @test S(1)*b1 == b1
   end
end

@testset "Julia.GFElem.adhoc_binary..." begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:1000
      a = rand(R)

      c1 = rand(0:100)
      c2 = rand(0:100)
      d1 = rand(BigInt(0):BigInt(100))
      d2 = rand(BigInt(0):BigInt(100))

      @test a + c1 == c1 + a
      @test a + d1 == d1 + a
      @test a - c1 == -(c1 - a)
      @test a - d1 == -(d1 - a)
      @test a*c1 == c1*a
      @test a*d1 == d1*a
      @test a*c1 + a*c2 == a*(c1 + c2)
      @test a*d1 + a*d2 == a*(d1 + d2)

      b = rand(S)

      @test b + c1 == c1 + b
      @test b + d1 == d1 + b
      @test b - c1 == -(c1 - b)
      @test b - d1 == -(d1 - b)
      @test b*c1 == c1*b
      @test b*d1 == d1*b
      @test b*c1 + b*c2 == b*(c1 + c2)
      @test b*d1 + b*d2 == b*(d1 + d2)
   end
end

@testset "Julia.GFElem.powering..." begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:100
      a = R(1)
      b = S(1)

      r = rand(R)
      s = rand(S)

      for n = 0:20
         @test r == 0 || a == r^n
         @test s == 0 || b == s^n

         a *= r
         b *= s
      end
   end

   for iter = 1:100
      a = R(1)
      b = S(1)

      r = rand(R)
      s = rand(S)

      rinv = r == 0 ? R(0) : inv(r)
      sinv = s == 0 ? S(0) : inv(s)

      for n = 0:20
         @test r == 0 || a == r^(-n)
         @test s == 0 || b == s^(-n)

         a *= rinv
         b *= sinv
      end
   end
end

@testset "Julia.GFElem.comparison..." begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:1000
      a = rand(R)
      b = rand(S)

      @test a != a + 1
      @test b != b + 1

      c = rand(0:100)
      d = rand(BigInt(0):BigInt(100))

      @test R(c) == R(c)
      @test S(c) == S(c)
      @test R(d) == R(d)
      @test S(d) == S(d)
   end
end

@testset "Julia.GFElem.adhoc_comparison..." begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:1000
      c = rand(0:100)
      d = rand(BigInt(0):BigInt(100))

      @test R(c) == c
      @test c == R(c)
      @test R(d) == d
      @test d == R(d)

      @test S(c) == c
      @test c == S(c)
      @test S(d) == d
      @test d == S(d)
   end
end

@testset "Julia.GFElem.inversion..." begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:1000
      a = rand(R)
      b = rand(S)

      @test a == 0 || inv(inv(a)) == a
      @test b == 0 || inv(inv(b)) == b

      @test a == 0 || a*inv(a) == one(R)
      @test b == 0 || b*inv(b) == one(S)
   end
end

@testset "Julia.GFElem.exact_division..." begin
   R = GF(13)
   S = GF(BigInt(13))

   for iter = 1:1000
      a1 = rand(R)
      a2 = rand(R)
      b1 = rand(S)
      b2 = rand(S)

      @test a2 == 0 || divexact(a1, a2)*a2 == a1
      @test b2 == 0 || divexact(b1, b2)*b2 == b1
   end
end
