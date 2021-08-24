@testset "Generic.UnivPoly.constructors" begin
   for R in [ZZ, QQ]
      for iters = 1:5
         ord = rand_ordering()

         @test UniversalPolynomialRing(R; ordering=ord, cached = true) === UniversalPolynomialRing(R; ordering=ord, cached = true)
         UniversalPolynomialRing(R; ordering=ord, cached = true) !== UniversalPolynomialRing(R; ordering=ord, cached = true)

         S = UniversalPolynomialRing(R; ordering=ord)

         x = gen(S, "x")

         @test isa(x, UnivPolyElem)

         y, z = gens(S, ["y", "z"])

         @test elem_type(S) == Generic.UnivPoly{elem_type(R), Generic.MPoly{elem_type(R)}}
         @test elem_type(Generic.UnivPolyRing{elem_type(R), Generic.MPoly{elem_type(R)}}) == Generic.UnivPoly{elem_type(R), Generic.MPoly{elem_type(R)}}
         @test parent_type(Generic.UnivPoly{elem_type(R), Generic.MPoly{elem_type(R)}}) == Generic.UnivPolyRing{elem_type(R), Generic.MPoly{elem_type(R)}}

         @test base_ring(S) === R
         @test coefficient_ring(S) === R

         @test typeof(S) <: Generic.UnivPolyRing

         @test isa(x, UnivPolyElem)
         @test isa(y, UnivPolyElem)
         @test isa(z, UnivPolyElem)

         f = rand(S, 0:5, 0:10, -10:10)

         @test isa(f, UnivPolyElem)
         @test isa(S(2), UnivPolyElem)
         @test isa(S(R(2)), UnivPolyElem)
         @test isa(S(f), UnivPolyElem)

         V = [R(rand(-10:10)) for i in 1:5]
         nused = rand(0:3)
         W = [[rand(0:100) for i in 1:nused] for i = 1:5]
         
         f1 = S(V, W)

         @test isa(f1, UnivPolyElem)

         f2 = S()
         fit!(f2, 5)

         for i = 1:5
            f2 = set_exponent_vector!(f2, i, W[i])
            f2 = setcoeff!(f2, i, V[i])
         end
         f2 = sort_terms!(f2)
         f2 = combine_like_terms!(f2)

         @test f1 == f2

         C = MPolyBuildCtx(S)

         for i = 1:5
            push_term!(C, V[i], W[i])
         end
         f3 = finish(C)

         @test f1 == f3

         @test x in [x, y]
         @test x in [y, x]
         @test y in [y, z]
         @test !(x in [y])
         @test !(y in [x])
         @test x in keys(Dict(x => 1))
         @test !(y in keys(Dict(x => 1)))
         @test !(y in keys(Dict(z => 1)))
      end
   end
end

@testset "Generic.UnivPoly.parent_type" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         ord = rand_ordering()

         S = UniversalPolynomialRing(R; ordering=ord, cached=false)

         @test symbols(S) == Symbol[]
         @test nvars(S) == 0

         f1 = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         f2 = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         f3 = rand(S, 0:5, 0:10, -10:10)

         @test isdomain_type(typeof(f1))
         @test isdomain_type(typeof(f2))
         @test isdomain_type(typeof(f3))

         @test isexact_type(typeof(f1))
         @test isexact_type(typeof(f2))
         @test isexact_type(typeof(f3))

         @test parent(f1) === S
         @test parent(f2) === S
         @test parent(f3) === S

         @test base_ring(S) === R
         @test base_ring(f1) === R
         @test base_ring(f2) === R
         @test base_ring(f3) === R

         @test nvars(S) == 3

         @test symbols(S) == [:x, :y, :z]

         @test length(vars(S())) == 0
         @test vars(x^2 + 2x + 1) == [x]
         @test vars(x*y + 1) == [x, y]
         @test vars(x*y + z) == [x, y, z]

         @test ordering(S) == ord
      end
   end
end

@testset "Generic.UnivPoly.term_monomial" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         ord = rand_ordering()

         S = UniversalPolynomialRing(R; ordering=ord, cached=false)

         x = gen(S, "x")
         y, z = gens(S, ["y", "z"])
         x2 = gen(S, "x")

         @test exponent_vector(x, 1) == [1]
         @test exponent_vector(x2, 1) == [1, 0, 0]
         @test exponent_vector(y + z, 1) == [0, 1, 0]
         @test exponent_vector(y + z, 2) == [0, 0, 1]

         @test exponent(x, 1, 1) == 1
         @test exponent(x, 1, 2) == 0
         @test exponent(y, 1, 1) == 0
         @test exponent(y, 1, 2) == 1

         g = deepcopy(x)
         g = set_exponent_vector!(g, 1, [2])

         @test g == x^2

         fit!(g, 2)
         g = set_exponent_vector!(g, 2, [1, 2, 3])
         g = setcoeff!(g, 2, R(3))
         g = sort_terms!(g)

         @test g == x^2 + 3x*y^2*z^3

         g = x^2*y*z
         g = setcoeff!(g, 1, 5)

         @test g == 5x^2*y*z

         h = deepcopy(x + y + z)
         h = set_exponent_vector!(h, 1, [2])

         @test h == x^2 + y + z

         @test coeff(x, [1, 2, 3]) == 0
         @test coeff(x, [1]) == 1
         @test coeff(x, [1, 0, 0]) == 1
         @test coeff(g, [2, 1, 1]) == 5
         @test coeff(h, [2]) == 1

         k = deepcopy(x)
         k = setcoeff!(k, [1], 4)

         @test k == 4x

         l = x*y^2
         l = setcoeff!(l, [1, 2], 5)

         @test l == 5x*y^2
      end
   end
end

@testset "Generic.UnivPoly.basic_manipulation" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = UniversalPolynomialRing(R; cached=false)

         @test iszero(zero(S))
         @test isone(one(S))
         @test isunit(one(S))
         @test !isgen(S())
         @test ishomogeneous(S())
         @test !ismonomial(S())
         @test isconstant(S())
         @test !isterm(S())
         @test trailing_coefficient(S()) == 0
         @test leading_coefficient(S()) == 0
         @test constant_coefficient(S()) == 0
         @test total_degree(S()) == -1
         @test length(gens(S)) == 0
         @test characteristic(S) == 0

         x = gen(S, "x")
         y, z = gens(S, ["y", "z"])
         x2 = gen(S, "x")

         @test !iszero(x)
         @test !iszero(x2)
         @test iszero(zero(S))
         @test !isone(x)
         @test !isone(x2)
         @test isone(one(S))
         @test !isunit(x)
         @test !isunit(x2)
         @test isunit(one(S))
         @test isgen(x)
         @test isgen(x2)
         @test isgen(y)
         @test ishomogeneous(x + y + z)
         @test ishomogeneous(x2)
         @test ishomogeneous(S())
         @test ismonomial(x^2)
         @test ismonomial(x2^2)
         @test ismonomial(y)
         @test !isconstant(x)
         @test !isconstant(x2)
         @test isconstant(S(1))
         @test isterm(2x^2*y)
         @test isterm(3x2)
         @test coeff(3x^2 + 2y + z, 2) == 2
         @test coeff(2x2^2 + 3x2 + 5, 1) == 2
         @test coeff(3x^2 + 2y + z, x2^2) == 3
         @test coeff(3x2^2, x^2) == 3
         @test coeff(3x^2 + 2y, x2) == 0
         @test coeff(3x2^2, x) == 0
         @test monomial(3x^2 + 2y + z, 1) == x^2
         @test monomial(2x2^2 + 2x2 + 4, 2) == x

         v = S()
         monomial!(v, 3x^2 + 2y + z, 1)

         @test v == x^2

         v = deepcopy(x2)
         monomial!(v, 3x^2 + 2y + z, 2)

         @test v == y

         @test term(3x^2 + 2y + z, 1) == 3x^2
         @test term(2x2^2 + 2x2 + 4, 2) == 2x
         @test leading_coefficient(3x^2 + 2y + z) == 3
         @test leading_coefficient(2x2^2 + 2x2 + 4) == 2
         @test trailing_coefficient(3x^2 + 2y + z) == 1
         @test trailing_coefficient(2x2^2 + 2x2) == 2
         @test tail(3x^2 + 2y + z) == 2y + z
         @test tail(2x2^2 + 2x2 + 4) == 2x + 4
         @test constant_coefficient(3x^2 + 2y + z + 4) == 4
         @test constant_coefficient(2x2^2 + 2x2) == 0
         @test leading_monomial(3x^2 + 2y + z + 4) == x^2
         @test leading_monomial(2x2^2 + 2x2) == x^2
         @test leading_term(3x^2 + 2y + z + 4) == 3x^2
         @test leading_term(2x2^2 + 2x2) == 2x^2
         @test max_fields(3x^2 + 2y + z + 4) == ([1, 1, 2], 2)
         @test max_fields(2x2^2 + 2x2) == ([0, 0, 2], 2)
         @test degree(3x^2 + 2y + z + 4, 1) == 2
         @test degree(2x2^2 + 2x2, 1) == 2
         @test degree(3x^2*y + 2x + 3y + 1, x) == 2
         @test degree(2x2^2 + 3x2 + 1, x) == 2
         @test degree(2x2^2 + 3x2 + 1, y) == 0
         @test degrees(3x^2*y + 2x + 3y + 1) == [2, 1, 0]
         @test degrees(2x2^2 + 2x2 + 1) == [2, 0, 0]
         @test total_degree(3x^2*y + 2x + 3y + 1) == 3
         @test total_degree(2x2^2 + 3x2 + 1) == 2
         @test length(3x^2*y + 2x + 3y + 1) == 4
         @test length(2x2^2 + 3x2 + 1) == 3

         @test gen(S, :x) == x
         @test gens(S, ['x', 'y']) == (x, y)
         @test gen(S, 1) == x
         @test gen(S, 2) == y
         @test gen(S, 3) == z
         @test gens(S) == [x, y, z]

         @test var_index(x) == 1
         @test var_index(x2) == 1
         @test var_index(y) == 2

         @test vars(3x^2*y + 2x + 3y + 1) == [x, y]
         @test vars(2x2^2 + 3x2 + 1) == [x]

         @test characteristic(S) == 0
      end
   end
end

@testset "Generic.UnivPoly.multivariate_coefficients" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = UniversalPolynomialRing(R; cached=false)

         x = gen(S, "x")
         y, z = gens(S, ["y", "z"])
         x2 = gen(S, "x")

         f = 3x^3 + 2x^2 + x + 4
         g = 3x^3*y^2 + 2x^3*y*z + 2x^2*y*z + 3x + 2y + 1

         @test coeff(f, [1], [3]) == 3
         @test coeff(f, [1, 2, 3], [3, 1, 2]) == 0
         @test coeff(f, [1, 2, 3], [3, 0, 0]) == 3
         @test coeff(g, [1], [3]) == 3y^2 + 2y*z
         @test coeff(g, [2], [1]) == 2x^3*z + 2x^2*z + 2
         @test coeff(g, [1, 2, 3], [0, 0, 0]) == 1

         @test coeff(f, [x], [3]) == 3
         @test coeff(f, [x2, y, z], [3, 1, 2]) == 0
         @test coeff(f, [x, y, z], [3, 0, 0]) == 3
         @test coeff(g, [x2], [3]) == 3y^2 + 2y*z
         @test coeff(g, [y], [1]) == 2x^3*z + 2x^2*z + 2
         @test coeff(g, [x2, y, z], [0, 0, 0]) == 1
      end
   end
end

@testset "Generic.UnivPoly.unary_operations" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = UniversalPolynomialRing(R; cached=false)

         f = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         g = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         h = rand(S, 0:5, 0:10, -10:10)

         @test f + (-f) == 0
         @test g + (-g) == 0
         @test h + (-h) == 0

         @test f + (-g) == f - g
         @test f + (-h) == f - h
         @test g + (-h) == g - h
         @test g + (-f) == g - f
         @test h + (-f) == h - f
         @test h + (-g) == h - g
      end
   end
end

@testset "Generic.UnivPoly.binary_operations" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = UniversalPolynomialRing(R; cached=false)

         f = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         g = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         h = rand(S, 0:5, 0:10, -10:10)

         @test f + g == g + f
         @test f + h == h + f
         @test g + h == h + g
         
         @test f - g == -(g - f)
         @test f - h == -(h - f)
         @test g - h == -(h - g)
         
         @test f*g == g*f
         @test f*h == h*f
         @test g*h == h*g
      end
   end
end

@testset "Generic.UnivPoly.iterators" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = UniversalPolynomialRing(R; cached=false)

         x = gen(S, "x")

         h = 3x^3 + 2x^2 + x + 4

         @test collect(exponent_vectors(h)) == [[3], [2], [1], [0]]

         y, z = gens(S, ["y", "z"])

         f = 3x^3 + 2x^2 + x + 4
         g = 3x^3*y^2 + 2x^3*y*z + 2x^2*y*z + 3x + 2y + 1

         @test collect(coefficients(f)) == [R(v) for v in [3, 2, 1, 4]]
         @test collect(exponent_vectors(f)) == [[3, 0, 0], [2, 0, 0], [1, 0, 0], [0, 0, 0]]
         @test sum(terms(f)) == f
         @test collect(monomials(f)) == [x^3, x^2, x, S(1)]

         @test collect(coefficients(g)) == [R(v) for v in [3, 2, 2, 3, 2, 1]]
         @test collect(exponent_vectors(g)) == [[3, 2, 0], [3, 1, 1], [2, 1, 1], [1, 0, 0], [0, 1, 0], [0, 0, 0]]
         @test sum(terms(g)) == g
         @test collect(monomials(g)) == [x^3*y^2, x^3*y*z, x^2*y*z, x, y, S(1)]
      end
   end
end

@testset "Generic.UnivPoly.square_root" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = UniversalPolynomialRing(R; cached=false)

         x = gen(S, "x")

         f = rand(S, 0:5, 0:10, -10:10)

         @test sqrt(f^2) == f || sqrt(f^2) == -f
         @test issquare(f^2)

         y, z = gens(S, ["y", "z"])

         g = rand(S, 0:5, 0:10, -10:10)

         @test sqrt(g^2) == g || sqrt(g^2) == -g
         @test issquare(g^2)
      end
   end
end

@testset "Generic.UnivPoly.adhoc_binary_operations" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = UniversalPolynomialRing(R; cached=false)

         f = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         g = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         h = rand(S, 0:5, 0:10, -10:10)

         c1 = rand(R, -10:10)
         c2 = rand(R, -10:10)
         d1 = rand(-10:10)
         d2 = rand(-10:10)

         @test c1*f + c2*f == (c1 + c2)*f
         @test f*c1 + f*c2 == (c1 + c2)*f
         @test c1*g + c2*g == (c1 + c2)*g
         @test g*c1 + g*c2 == (c1 + c2)*g
         @test c1*h + c2*h == (c1 + c2)*h
         @test h*c1 + h*c2 == (c1 + c2)*h

         @test d1*f + d2*f == (d1 + d2)*f
         @test f*d1 + f*d2 == (d1 + d2)*f
         @test d1*g + d2*g == (d1 + d2)*g
         @test g*d1 + g*d2 == (d1 + d2)*g
         @test d1*h + d2*h == (d1 + d2)*h
         @test h*d1 + h*d2 == (d1 + d2)*h

         @test R(d1)*f + R(d2)*f == (d1 + d2)*f
         @test f*R(d1) + f*R(d2) == (d1 + d2)*f
         @test R(d1)*g + R(d2)*g == (d1 + d2)*g
         @test g*R(d1) + g*R(d2) == (d1 + d2)*g
         @test R(d1)*h + R(d2)*h == (d1 + d2)*h
         @test h*R(d1) + h*R(d2) == (d1 + d2)*h

         if !iszero(c1)
            @test divexact(c1*f, c1) == f
            @test divexact(c1*g, c1) == g
            @test divexact(c1*h, c1) == h
         end

         if !iszero(d1)
            @test divexact(d1*f, d1) == f
            @test divexact(d1*g, d1) == g
            @test divexact(d1*h, d1) == h

            @test divexact(d1*f, R(d1)) == f
            @test divexact(d1*g, R(d1)) == g
            @test divexact(d1*h, R(d1)) == h
         end
      end
   end
end

@testset "Generic.UnivPoly.comparison" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = UniversalPolynomialRing(R; cached=false)

         f = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         g = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         h = rand(S, 0:5, 0:10, -10:10)

         @test f == f
         @test g == g
         @test h == h
         @test f == deepcopy(f)
         @test g == deepcopy(g)
         @test h == deepcopy(h)

         if !iszero(g)
             @test f != g + f
             @test g + f != f
         end
         if !iszero(h)
             @test g != g + h
             @test g + h != g
             @test f != f + h
             @test f + h != f
         end

         @test y < x
         @test z < x
         @test z < y
         @test x^2*y*z < x^2*y*z^5
         @test x < x^2
         @test x < x^2*y
         @test y^2 < y^2*z^3
      end
   end
end

