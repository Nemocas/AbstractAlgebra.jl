# TODO: instead of 'example_element' provide a function which allows sampling any number
# of elements
function test_Ring_interface(R::AbstractAlgebra.Ring, ExpectedElemType, ExpectedParentType, example_elements = [])
   @testset "Ring interface for $(R) / $(ExpectedElemType) / $(ExpectedParentType)" begin
      @test R isa ExpectedParentType
      example_elements = R.(example_elements)

      @testset "test ring conformance" begin
         @test elem_type(R) == ExpectedElemType
         @test elem_type(ExpectedParentType) == ExpectedElemType
         @test elem_type(typeof(R)) == ExpectedElemType

         @test iszero(zero(R))
         @test isone(one(R))

         @test iszero(R())
         @test iszero(R(0))
         @test isone(R(1))

         @test !isunit(R())
         @test isunit(R(1))

         @test iszero(R(characteristic(R)))
         @test iszero(characteristic(R) * one(R))
         @test iszero(one(R) * characteristic(R))
      end

      elems = [ zero(R), one(R), R(R(2)), example_elements... ]
      @testset "test element conformance" for x in elems
         @test x isa ExpectedElemType
         @test parent_type(x) == ExpectedParentType
         @test parent_type(typeof(x)) == ExpectedParentType
      end
   end
end

test_Ring_interface(zz, Int, AbstractAlgebra.Integers{Int}, [11])
test_Ring_interface(ZZ, BigInt, AbstractAlgebra.Integers{BigInt}, [11])

test_Ring_interface(qq, Rational{Int}, AbstractAlgebra.Rationals{Int}, [11])
test_Ring_interface(QQ, Rational{BigInt}, AbstractAlgebra.Rationals{BigInt}, [11])

test_Ring_interface(RDF, Float64, AbstractAlgebra.Floats{Float64},
        [11, 11//3, 1.2])
test_Ring_interface(RealField, BigFloat, AbstractAlgebra.Floats{BigFloat},
        [11, 11/3, 1.2, big(11), big(11)//big(3), big(1.2), Rational{BigInt}(11)])

test_Ring_interface(GF(3), AbstractAlgebra.GFElem{Int}, AbstractAlgebra.GFField{Int}, [2])
test_Ring_interface(GF(13), AbstractAlgebra.GFElem{Int}, AbstractAlgebra.GFField{Int}, [11])
test_Ring_interface(GF(big(13)), AbstractAlgebra.GFElem{BigInt}, AbstractAlgebra.GFField{BigInt}, [11, big(11)])


test_Ring_interface(ZZ["x"][1], Generic.Poly{BigInt}, Generic.PolyRing{BigInt})
