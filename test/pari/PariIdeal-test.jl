function test_PariIdeal_constructors()
   print("PariIdeal.constructors...")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = parent(Ideal(O, 2, 3))

   @test isa(S, PariIdealSet)

   @test isa(Ideal(O, 2, 3), PariIdeal)

   @test isa(Ideal(O, x + 1, R(3)), PariIdeal)

   @test isa(Ideal(O, a + 1, K(3)), PariIdeal)

   @test isa(Ideal(O, O(a + 1), O(3)), PariIdeal)

   println("PASS")
end

function test_PariIdeal_lists()
   print("PariIdeal.lists...")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = prime_decomposition(O, 13)
   T = bounded_ideals(O, 100)

   @test length(T) == 100

   a = S[2]
   b = T[97][1]

   @test isa(a, PariIdeal)
   @test isa(b, PariIdeal)

   println("PASS")
end

function test_PariIdeal_basic_manipulation()
   print("PariIdeal.basic_manipulation...")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = Ideal(O, x + 1)
   T = Ideal(O, fmpq(2, 3)*x + 1)
   A = prime_decomposition(O, 3)

   @test valuation(S, A[1]) == 1
   @test numden(T) == (Ideal(O, R(55), x-26), Ideal(O, 3))

   println("PASS")
end

function test_PariIdeal_norm()
   print("PariIdeal.norm...")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = Ideal(O, x + 1)

   @test norm(S) == 3

   println("PASS")
end

function test_PariIdeal_binary_ops()
   print("PariIdeal.binary_ops...")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = Ideal(O, x + 1)
   T = Ideal(O, 3)

   @test S + T == S
   @test S*T == S^4
   @test intersect(S, T) == Ideal(O, 3)
   println("PASS")
end

function test_PariIdeal_comparison()
   print("PariIdeal.comparison...")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = Ideal(O, x + 1)
   T = Ideal(O, 3)

   @test S == S
   @test S != T
   @test Ideal(O, 3) == Ideal(O, R(3))

   println("PASS")
end

function test_PariIdeal_powering()
   print("PariIdeal.powering...")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = Ideal(O, x + 1)
   
   @test S^3 == Ideal(O, 3)

   println("PASS")
end

function test_PariIdeal_exact_division()
   print("PariIdeal.exact_division..")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = Ideal(O, x + 1)
   T = Ideal(O, 3)

   @test divexact(S*T, T) == S

   println("PASS")
end

function test_PariIdeal_inversion()
   print("PariIdeal.inversion..")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = Ideal(O, x + 1)

   @test S*inv(S) == Ideal(O, 1)

   println("PASS")
end

function test_PariIdeal_extended_gcd()
   print("PariIdeal.extended_gcd..")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = Ideal(O, x + 1)
   T = Ideal(O, 2x)

   @test gcdx(S, T) == (2*x-1,-2*x+2)

   println("PASS")
end

function test_PariIdeal_factorisation()
   print("PariIdeal.factorisation..")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = Ideal(O, 3)

   A = factor(S)
   @test factor_mul(A) == S

   println("PASS")
end

function test_PariIdeal_approximation()
   print("PariIdeal.approximation..")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = Ideal(O, 3)

   @test approx(S) == -3x

   println("PASS")
end

function test_PariIdeal_coprime_multiplier()
   print("PariIdeal.coprime_multiplier..")
 
   R, x = PolynomialRing(QQ, "x")
   K, a = NumberField(x^3 + 3x + 1, "a")
   O = MaximalOrder(K)

   S = Ideal(O, 3)
   T = Ideal(O, x + 1)

   @test coprime_multiplier(S, T) == -fmpq(1, 3)*x

   println("PASS")
end

function test_PariIdeal()
   test_PariIdeal_constructors()
   test_PariIdeal_lists()
   test_PariIdeal_basic_manipulation()
   test_PariIdeal_norm()
   test_PariIdeal_binary_ops()
   test_PariIdeal_comparison()
   test_PariIdeal_powering()
   test_PariIdeal_exact_division()
   test_PariIdeal_inversion()
   test_PariIdeal_factorisation()
   test_PariIdeal_approximation()
   test_PariIdeal_coprime_multiplier()

   println("")
end
