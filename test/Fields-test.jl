include("generic/Fraction-test.jl")

include("flint/fmpq-test.jl")
include("flint/fq-test.jl")
include("flint/fq_nmod-test.jl")

include("Padics-test.jl")

function test_fields_resultant()
   print("Fields.resultant...")

   R, x = FiniteField(7, 11, "x")
   S, y = PolynomialRing(R, "y")
   T = ResidueRing(S, y^3 + 3x*y + 1)
   U, z = PolynomialRing(T, "z")

   f = (3y^2 + y + x)*z^2 + ((x + 2)*y^2 + x + 1)*z + 4x*y + 3
   g = (7y^2 - y + 2x + 7)*z^2 + (3y^2 + 4x + 1)*z + (2x + 1)*y + 1

   s = f^12
   t = (s + g)^12
   
   r = resultant(s, t)

   @test r == (x^10+4*x^8+6*x^7+3*x^6+4*x^5+x^4+6*x^3+5*x^2+x)*y^2+(5*x^10+x^8+4*x^7+3*x^5+5*x^4+3*x^3+x^2+x+6)*y+(2*x^10+6*x^9+5*x^8+5*x^7+x^6+6*x^5+5*x^4+4*x^3+x+3)

   println("PASS")
end

function test_fields()
   test_fraction()

   test_fmpq()
   test_fq()
   test_fq_nmod()
   
   test_padics()
   
   test_fields_resultant()
end
