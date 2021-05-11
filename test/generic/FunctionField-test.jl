R1, x1 = RationalFunctionField(QQ, "x1")
U1, z1 = R1["z1"]

P1 = [(x1 + 1)*z1 + (x1 + 2), z1 + (x1 + 1)//(x1 + 2), z1^2 + 3z1 + 1,
     (x1^2 + 1)//(x1 + 1)*z1^5 + 4z1^4 + (x1 + 2)*z1^3 + x1//(x1 + 1)*z1 + 1//(z1 + 1)]

R2, x2 = RationalFunctionField(GF(23), "x2")

U2, z2 = R2["z2"]

P2 = [(x2 + 1)*z2 + (x2 + 2), z2 + (x2 + 1)//(x2 + 2), z2^2 + 3z2 + 1,
     (x2^2 + 1)//(x2 + 1)*z2^5 + 4z2^4 + (x2 + 2)*z2^3 + x2//(x2 + 1)*z2 + 1//(z2 + 1)]

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

