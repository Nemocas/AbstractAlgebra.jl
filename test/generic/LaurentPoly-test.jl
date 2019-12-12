using AbstractAlgebra.Generic: Integers, LaurentPolyWrapRing, LaurentPolyWrap

@testset "Generic.LaurentPoly" begin
   @testset "constructors" begin
      for R in (ZZ, GF(5))
         L, y = LaurentPolynomialRing(R, "y")

         @test L isa LaurentPolyWrapRing{elem_type(R)}
         @test y isa LaurentPolyWrap{elem_type(R)}

         @test parent_type(y) == typeof(L)
         @test elem_type(L) == typeof(y)

         @test parent(y) == L

         @test var(L) == :y
         @test symbols(L) == [:y]
         @test nvars(L) == 1

         @test characteristic(L) == characteristic(R)
      end
   end
end
