function test_fmpz_mod_poly_modular_arithmetic()
   print("fmpz_mod_poly.modular_arithmetic...")

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

function test_fmpz_mod_poly_euclidean_division()
   print("fmpz_mod_poly.euclidean_division...")

   R = ResidueRing(ZZ, 7)
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test mod(g, f) == 6*x+3

   @test divrem(g, f) == (x+5, 6*x+3)
 
   println("PASS")
end

function test_fmpz_mod_poly_gcdinv()
   print("fmpz_mod_poly.gcdinv...")

   R = ResidueRing(ZZ, 7)
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   g = x^3 + 3x + 1

   @test gcdinv(f, g) == (1,3*x^2+2*x+2)
 
   println("PASS")
end

function test_fmpz_mod_poly_integral()
   print("fmpz_mod_poly.integral...")

   R = ResidueRing(ZZ, 7)
   S, x = PolynomialRing(R, "x")

   f = x^2 + 2x + 1
   
   @test integral(f) == 5x^3 + x^2 + x

   println("PASS")
end

function test_fmpz_mod_poly()
   test_fmpz_mod_poly_modular_arithmetic()
   test_fmpz_mod_poly_euclidean_division()
   test_fmpz_mod_poly_gcdinv()
   test_fmpz_mod_poly_integral()

   println("")
end

