function test_constructors()
   print("Poly.constructors...")
 
   R, x = PolynomialRing(ZZ, "x")

   @test R <: Poly

   @test isa(x, Poly)

   S, y = PolynomialRing(R, "y")

   @test S <: Poly

   @test isa(y, Poly)

   T, z = PolynomialRing(S, "z")

   @test T <: Poly

   @test isa(z, Poly)

   f = x^2 + y^3 + z + 1

   @test isa(f, Poly)

   g = S(2)

   @test isa(g, Poly)

   h = S(x^2 + 2x + 1)

   @test isa(h, Poly)

   j = T(x + 2)

   @test isa(j, Poly)

   k = Poly(S, [x, x + 2, x^2 + 3x + 1])

   @test isa(k, Poly)

   l = S(k)

   @test isa(l, Poly)

   println("PASS")
end

function test_manipulation()
   print("Poly.manipulation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   @test iszero(zero(R))
   
   @test isone(one(S))

   @test isgen(gen(S))
   
   @test isunit(one(S))

   f = x^2 + 2x + 1

   @test lead(f) == 1

   h = x*y^2 + (x + 1)*y + 3

   @test coeff(h, 2) == x

   @test canonical_unit(-x + 1) == -1

   @test canonical_unit(-x*y + x + 1) == -1

   println("PASS")
end

function test_binary_ops()
   print("Poly.binary_ops...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 2

   @test f + g == x^3+x^2+5*x+3

   @test f*g == x^5+2*x^4+4*x^3+8*x^2+7*x+2

   k = x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y + (x^3 + 2x + 2)
   
   @test k - l == x*y^2+(-x^3-2*x+1)

   @test k*(f + g) == (x^4+x^3+5*x^2+3*x)*y^2+(x^4+2*x^3+6*x^2+8*x+3)*y+(3*x^3+3*x^2+15*x+9)

   @test (k - l) - g == x*y^2+-2*x^3-5*x-1

   println("PASS")
end

function test_adhoc_binary()
   print("Poly.adhoc_binary...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 2

   @test f + 12 == x^2+2*x+13

   @test ZZ(3)*g == 3*x^3+9*x+6

   k = x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y + (x^3 + 2x + 2)

   @test k - 4 == x*y^2+(x+1)*y-1
   
   @test ZZ(5) - l == (-x-1)*y+(-x^3-2*x+3)

   println("PASS")
end

function test_comparison()
   print("Poly.comparison...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 2

   @test f != g

   k = x*y^2 + (x + 1)*y + 3
   l = x*y^2 + (x + 1)*y + 3

   @test k == l

   println("PASS")
end

function test_adhoc_comparison()
   print("Poly.adhoc_comparison...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1

   @test f != 1 

   @test S(1) == ZZ(1) 

   println("PASS")
end

function test_unary_ops()
   print("Poly.unary_ops...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1

   @test -f == -x^2 - 2x - 1

   h = x*y^2 + (x + 1)*y + 3

   @test -h == -x*y^2 - (x + 1)*y - 3

   println("PASS")
end

function test_truncation()
   print("Poly.truncation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test truncate(f, 2) == 2*x+1

   @test mullow(f, g, 3) == 7*x^2+5*x+1

   k = x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y + (x^3 + 2x + 2)

   @test truncate(k, 1) == 3

   @test mullow(k, l, 4) == (x^2+x)*y^3+(x^4+3*x^2+4*x+1)*y^2+(x^4+x^3+2*x^2+7*x+5)*y+(3*x^3+6*x+6)

   println("PASS")
end

function test_reverse()
   print("Poly.reverse...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 3

   @test reverse(f) == 3x^2 + 2x + 1

   k = x*y^2 + (x + 1)*y + 3

   @test reverse(k, 7) == 3y^6 + (x + 1)*y^5 + x*y^4

   println("PASS")
end

function test_shift()
   print("Poly.shift...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1

   @test shift_left(f, 3) == x^5 + 2x^4 + x^3

   @test shift_right(f, 1) == x + 2

   k = x*y^2 + (x + 1)*y + 3

   @test shift_left(k, 7) == x*y^9 + (x + 1)*y^8 + 3y^7

   @test shift_right(k, 3) == 0

   println("PASS")
end

function test_powering()
   print("Poly.powering...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1

   @test f^12 == x^24+24*x^23+276*x^22+2024*x^21+10626*x^20+42504*x^19+134596*x^18+346104*x^17+735471*x^16+1307504*x^15+1961256*x^14+2496144*x^13+2704156*x^12+2496144*x^11+1961256*x^10+1307504*x^9+735471*x^8+346104*x^7+134596*x^6+42504*x^5+10626*x^4+2024*x^3+276*x^2+24*x+1 

   k = x*y^2 + (x + 1)*y + 3

   @test k^5 == (x^5)*y^10+(5*x^5+5*x^4)*y^9+(10*x^5+35*x^4+10*x^3)*y^8+(10*x^5+90*x^4+90*x^3+10*x^2)*y^7+(5*x^5+110*x^4+300*x^3+110*x^2+5*x)*y^6+(x^5+65*x^4+460*x^3+460*x^2+65*x+1)*y^5+(15*x^4+330*x^3+900*x^2+330*x+15)*y^4+(90*x^3+810*x^2+810*x+90)*y^3+(270*x^2+945*x+270)*y^2+(405*x+405)*y+243

   println("PASS")
end

function test_modular_arithmetic()
   print("Poly.modular_arithmetic...")

   R = ResidueRing(ZZ, 487326487)
   S, x = PolynomialRing(R, "x")

   f = 3*x^2 + x + 2
   g = 5*x^2 + 2*x + 1
   h = 3*x^3 + 2*x^2 + x + 7

   @test invmod(f, h) == 40508247*x^2+341251293*x+416130174 

   @test mulmod(f, g, h) == 324884334*x^2+162442132*x+162442162

   @test powmod(f, 10, h) == 485924368*x^2+380106591*x+302530457

   println("PASS")
end

function test_exact_division()
   print("Poly.exact_division...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test divexact(f*g, f) == g

   @test divexact(3*f, ZZ(3)) == f

   k = x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y + (x^3 + 2x + 2)

   @test divexact(k*l, k) == l

   @test divexact(3*k, 3) == k

   @test divexact(x*k, x) == k

   println("PASS")
end

function test_euclidean_division()
   print("Poly.euclidean_division...")

   R = ResidueRing(ZZ, 7)
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test mod(g, f) == 6*x+3

   @test divrem(g, f) == (x+5, 6*x+3)

   T = ResidueRing(S, x^3 + 3x + 1)
   U, y = PolynomialRing(T, "y")

   k = y^3 + x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y^2 + (x^3 + 2x + 2)

   @test mod(k, l) == (4*x^2+4*x+4)*y+(3*x^2+5*x+6)
   
   @test divrem(k, l) == ((5*x^2+2*x+6)*y+(2*x^2+5*x+2), (4*x^2+4*x+4)*y+(3*x^2+5*x+6))
 
   println("PASS")
end

function test_pseudodivision()
   print("Poly.pseudodivision...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test pseudorem(f, g) == x^2+2*x+1

   @test pseudodivrem(f, g) == (0, x^2+2*x+1)

   k = x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y + (x^3 + 2x + 2)

   @test pseudorem(k, l) == (x^7+3*x^5+2*x^4+x^3+5*x^2+4*x+1)

   @test pseudodivrem(k, l) == ((x^2+x)*y+(-x^4-x^2+1), (x^7+3*x^5+2*x^4+x^3+5*x^2+4*x+1))

   println("PASS")
end

function test_content_primpart_gcd()
   print("Poly.content_primpart_gcd...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1
   h = x + 1

   @test content(3*f) == 3

   @test primpart(3*f) == f

   @test gcd(f*h, g*h) == h

   @test lcm(f*h, g*h) == f*g*h

   k = x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y + (x^3 + 2x + 2)
   m = y^2 + x + 1

   @test content(k) == 1

   @test primpart(k*(x^2 + 1)) == k

   @test gcd(k*m, l*m) == m

   @test lcm(k*m, l*m) == k*l*m

   T = ResidueRing(ZZ, 7)
   U, z = PolynomialRing(T, "z")

   r = z^3 + 2z + 1
   s = z^5 + 1

   @test gcdinv(r, s) == (1, z^3+4*z^2+5*z+5)

   println("PASS")
end

function test_evaluation()
   print("Poly.evaluation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1

   @test evaluate(f, 3) == 16
   
   @test evaluate(f, ZZ(10)) == 121

   j = x*y^2 + (x + 1)*y + 3

   @test evaluate(j, 3) == 12x + 6

   @test evaluate(j, f) == x^5+4*x^4+7*x^3+7*x^2+4*x+4

   println("PASS")
end

function test_composition()
   print("Poly.composition...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test compose(f, g) == x^6+6*x^4+4*x^3+9*x^2+12*x+4

   j = x*y^2 + (x + 1)*y + 3
   k = (x + 1)*y + (x^3 + 2x + 2)

   @test compose(j, k) == (x^3+2*x^2+x)*y^2+(2*x^5+2*x^4+4*x^3+9*x^2+6*x+1)*y+(x^7+4*x^5+5*x^4+5*x^3+10*x^2+8*x+5)

   println("PASS")
end

function test_derivative()
   print("Poly.derivative...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1

   @test deriv(f) == 2x + 2

   h = x*y^2 + (x + 1)*y + 3

   @test deriv(h) == 2x*y + x + 1

   println("PASS")
end

function test_integral()
   print("Poly.integral...")

   R = ResidueRing(ZZ, 7)
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   
   @test integral(f) == 5x^3 + x^2 + x

   println("PASS")
end

function test_resultant()
   print("Poly.resultant...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = 5x^2 + 2x + 1
   g = 2x^3 + 3x + 1

   @test resultant(f, g) == 212

   j = 3x*y^2 + (x + 1)*y + 3
   k = 6(x + 1)*y + (x^3 + 2x + 2)

   @test resultant(j, k) == 3*x^7+6*x^5-6*x^3+96*x^2+192*x+96

   println("PASS")
end

function test_discriminant()
   print("Poly.discriminant...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1

   @test discriminant(f) == 0

   h = x*y^2 + (x + 1)*y + 3

   @test discriminant(h) == x^2-10*x+1

   println("PASS")
end

function test_bezout()
   print("Poly.bezout...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = 5x^2 + 2x + 1
   g = 2x^3 + 3x + 1

   @test bezout(f, g) == (212, 146*x^2-58*x+213, -365*x-1)

   j = 3x*y^2 + (x + 1)*y + 3
   k = 6(x + 1)*y + (x^3 + 2x + 2)

   @test bezout(j, k) == (3*x^7+6*x^5-6*x^3+96*x^2+192*x+96, (36*x^2+72*x+36), (-18*x^2-18*x)*y+(3*x^4-6*x-6))

   println("PASS")
end

function test_special()
   print("Poly.special...")

   R, x = PolynomialRing(ZZ, "x")

   @test chebyshev_t(R, 20) == 524288*x^20-2621440*x^18+5570560*x^16-6553600*x^14+4659200*x^12-2050048*x^10+549120*x^8-84480*x^6+6600*x^4-200*x^2+1

   @test chebyshev_u(R, 15) == 32768*x^15-114688*x^13+159744*x^11-112640*x^9+42240*x^7-8064*x^5+672*x^3-16*x
   
   @test cyclotomic(R, 120) == x^32+x^28-x^20-x^16-x^12+x^4+1

   @test swinnerton_dyer(R, 5) == x^32-448*x^30+84864*x^28-9028096*x^26+602397952*x^24-26625650688*x^22+801918722048*x^20-16665641517056*x^18+239210760462336*x^16-2349014746136576*x^14+15459151516270592*x^12-65892492886671360*x^10+172580952324702208*x^8-255690851718529024*x^6+183876928237731840*x^4-44660812492570624*x^2+2000989041197056

   @test cos_minpoly(R, 30) == x^4+x^3-4*x^2-4*x+1

   @test theta_qexp(R, 3, 30) == 72*x^29+32*x^27+72*x^26+30*x^25+24*x^24+24*x^22+48*x^21+24*x^20+24*x^19+36*x^18+48*x^17+6*x^16+48*x^14+24*x^13+8*x^12+24*x^11+24*x^10+30*x^9+12*x^8+24*x^6+24*x^5+6*x^4+8*x^3+12*x^2+6*x+1

   @test eta_qexp(R, 24, 30) == -29211840*x^29+128406630*x^28+24647168*x^27-73279080*x^26+13865712*x^25-25499225*x^24+21288960*x^23+18643272*x^22-12830688*x^21-4219488*x^20-7109760*x^19+10661420*x^18+2727432*x^17-6905934*x^16+987136*x^15+1217160*x^14+401856*x^13-577738*x^12-370944*x^11+534612*x^10-115920*x^9-113643*x^8+84480*x^7-16744*x^6-6048*x^5+4830*x^4-1472*x^3+252*x^2-24*x+1

   println("PASS")
end

function test_poly()
   test_constructors()
   test_manipulation()
   test_binary_ops()
   test_adhoc_binary()
   test_comparison()
   test_adhoc_comparison()
   test_unary_ops()
   test_truncation()
   test_reverse()
   test_shift()
   test_powering()
   test_modular_arithmetic()
   test_exact_division()
   test_euclidean_division()
   test_pseudodivision()
   test_content_primpart_gcd()
   test_evaluation()
   test_composition()
   test_derivative()
   test_integral()
   test_resultant()
   test_discriminant()
   test_bezout()
   test_special()
end
