@testset "Generic.UnivPoly.constructors" begin
   for R in [ZZ, QQ]
      for iters = 1:5
         ord = rand_ordering()

         @test universal_polynomial_ring(R; internal_ordering=ord, cached = true) === universal_polynomial_ring(R; internal_ordering=ord, cached = true)
         @test universal_polynomial_ring(R; internal_ordering=ord, cached = false) !== universal_polynomial_ring(R; internal_ordering=ord, cached = false)

         S = universal_polynomial_ring(R; internal_ordering=ord)

         x = gen(S, "x")
         @test x == S[1]

         @test_throws ArgumentError gen(S, 0)
         @test_throws ArgumentError gen(S, nvars(S) + 1)

         @test isa(x, UniversalPolyRingElem)

         y, z = gens(S, ["y", "z"])
         @test y == S[2]

         @test elem_type(S) == Generic.UnivPoly{elem_type(R)}
         @test elem_type(Generic.UniversalPolyRing{elem_type(R)}) == Generic.UnivPoly{elem_type(R)}
         @test parent_type(Generic.UnivPoly{elem_type(R)}) == Generic.UniversalPolyRing{elem_type(R)}

         @test base_ring(S) === R
         @test coefficient_ring(S) === R

         @test typeof(S) <: Generic.UniversalPolyRing

         @test isa(x, UniversalPolyRingElem)
         @test isa(y, UniversalPolyRingElem)
         @test isa(z, UniversalPolyRingElem)

         f = rand(S, 0:5, 0:10, -10:10)

         @test isa(f, UniversalPolyRingElem)
         @test isa(S(2), UniversalPolyRingElem)
         @test isa(S(R(2)), UniversalPolyRingElem)
         @test isa(S(f), UniversalPolyRingElem)

         V = [R(rand(-10:10)) for i in 1:5]
         nused = rand(0:3)
         W = [[rand(0:100) for i in 1:nused] for i = 1:5]
         
         f1 = S(V, W)

         @test isa(f1, UniversalPolyRingElem)

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

         S = universal_polynomial_ring(R; internal_ordering=ord, cached=false)

         @test symbols(S) == Symbol[]
         @test nvars(S) == 0

         f1 = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         f2 = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         f3 = rand(S, 0:5, 0:10, -10:10)

         @test is_domain_type(typeof(f1))
         @test is_domain_type(typeof(f2))
         @test is_domain_type(typeof(f3))

         @test is_exact_type(typeof(f1))
         @test is_exact_type(typeof(f2))
         @test is_exact_type(typeof(f3))

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

         @test internal_ordering(S) == ord
      end
   end
end

function test_elem(R::AbstractAlgebra.Generic.UniversalPolyRing{EuclideanRingResidueRingElem{BigInt}})
    return rand(R, 0:4, 0:10, -10:10)
end

@testset "Generic.UnivPoly.conformance" begin
   S = universal_polynomial_ring(residue_ring(ZZ, ZZ(6))[1])
   gen(S, "x")
   test_Ring_interface(S)
end

@testset "Generic.UnivPoly.printing" begin
   R = universal_polynomial_ring(QQ)
   x = gen(R, "x")
   @test sprint(show, x+1) == "x + 1"
   @test sprint(show, -1 + x^2 - 2x + 1) == "x^2 - 2*x"

   S, y = polynomial_ring(R, "y")
   @test sprint(show, (x+1)*y) == "(x + 1)*y"

   z = gen(R, "z")
   @test sprint(show, x*y*(z-1)) == "(x*z - x)*y"
end

@testset "Generic.UnivPoly.term_monomial" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         ord = rand_ordering()

         S = universal_polynomial_ring(R; internal_ordering=ord, cached=false)

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
         S = universal_polynomial_ring(R; cached=false)

         @test iszero(zero(S))
         @test isone(one(S))
         @test is_unit(one(S))
         @test !is_gen(S())
         @test is_homogeneous(S())
         @test !is_monomial(S())
         @test is_constant(S())
         @test !is_term(S())
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
         @test !is_unit(x)
         @test !is_unit(x2)
         @test is_unit(one(S))
         @test is_zero_divisor(S())
         @test !is_zero_divisor(x)
         @test is_gen(x)
         @test is_gen(x2)
         @test is_gen(y)
         @test is_homogeneous(x + y + z)
         @test is_homogeneous(x2)
         @test is_homogeneous(S())
         @test is_monomial(x^2)
         @test is_monomial(x2^2)
         @test is_monomial(y)
         @test !is_constant(x)
         @test !is_constant(x2)
         @test is_constant(S(1))
         @test is_term(2x^2*y)
         @test is_term(3x2)
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

         xs, ys = gens(S, :x => (1:4), :y =>(1:2, 1:3))
         @test xs isa Vector{elem_type(S)}
         @test length(xs) == 4
         @test ys isa Matrix{elem_type(S)}
         @test size(ys) == (2, 3)
      end
   end
end

@testset "Generic.UnivPoly.multivariate_coefficients" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = universal_polynomial_ring(R; cached=false)

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
         S = universal_polynomial_ring(R; cached=false)

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
         S = universal_polynomial_ring(R; cached=false)

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
         S = universal_polynomial_ring(R; cached=false)

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
         S = universal_polynomial_ring(R; cached=false)

         x = gen(S, "x")

         f = rand(S, 0:5, 0:10, -10:10)

         @test sqrt(f^2) == f || sqrt(f^2) == -f
         @test is_square(f^2)

         y, z = gens(S, ["y", "z"])

         g = rand(S, 0:5, 0:10, -10:10)

         @test sqrt(g^2) == g || sqrt(g^2) == -g
         @test is_square(g^2)
      end
   end
end

@testset "Generic.UnivPoly.adhoc_binary_operations" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = universal_polynomial_ring(R; cached=false)

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
         S = universal_polynomial_ring(R; cached=false)

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

@testset "Generic.UnivPoly.adhoc_comparison" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = universal_polynomial_ring(R; cached=false)

         f = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         g = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         h = rand(S, 0:5, 0:10, -10:10)

         c = rand(R, -10:10)
         d = rand(-10:10)

         @test S(c) == c
         @test c == S(c)
         @test S(d) == d
         @test d == S(d)
         @test S(d) == R(d)
         @test R(d) == S(d)

         if !is_constant(f)
            @test f != c
            @test c != f
            @test f != d
            @test d != f
            @test f != R(d)
            @test R(d) != f
         end

         if !is_constant(g)
            @test g != c
            @test c != g
            @test g != d
            @test d != g
            @test g != R(d)
            @test R(d) != g
         end

         if !is_constant(h)
            @test h != c
            @test c != h
            @test h != d
            @test d != h
            @test h != R(d)
            @test R(d) != h
         end
      end
   end
end

@testset "Generic.UnivPoly.powering" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = universal_polynomial_ring(R; cached=false)

         f = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         g = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         h = rand(S, 0:5, 0:10, -10:10)

         pf = S(1)
         pg = S(1)
         ph = S(1)

         for i = 0:5
            @test pf == f^i
            @test pg == g^i
            @test ph == h^i

            pf *= f
            pg *= g
            ph *= h
         end
      end
   end
end

@testset "Generic.UnivPoly.inflation_deflation" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = universal_polynomial_ring(R; cached=false)

         f = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         g = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         h = rand(S, 0:5, 0:10, -10:10)

         n = rand(0:3)
         shift = [rand(0:5) for i in 1:n]
         defl = [rand(1:5) for i in 1:n]
         
         fi = inflate(f, shift, defl)

         @test deflate(fi, shift, defl) == f        
         @test deflate(inflate(f, defl), defl) == f

         shift2, defl2 = deflation(fi)

         @test inflate(deflate(fi, shift2, defl2), shift2, defl2) == fi

         gi = inflate(g, shift, defl)

         @test deflate(gi, shift, defl) == g
         @test deflate(inflate(g, defl), defl) == g

         shift2, defl2 = deflation(gi)

         @test inflate(deflate(gi, shift2, defl2), shift2, defl2) == gi

         hi = inflate(h, shift, defl)

         @test deflate(hi, shift, defl) == h
         @test deflate(inflate(h, defl), defl) == h

         shift2, defl2 = deflation(hi)

         @test inflate(deflate(hi, shift2, defl2), shift2, defl2) == hi
      end
   end
end

@testset "Generic.UnivPoly.exact_division" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = universal_polynomial_ring(R; cached=false)

         f = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         g = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         h = rand(S, 0:5, 0:10, -10:10)

         if !iszero(f)
            @test divexact(f, f) == 1
            @test divexact(f*g, f) == g
            @test divexact(f*h, f) == h
         end

         if !iszero(g)
            @test divexact(g, g) == 1
            @test divexact(g*f, g) == f
            @test divexact(g*h, g) == h
         end

         if !iszero(h)
            @test divexact(h, h) == 1
            @test divexact(h*f, h) == f
            @test divexact(h*g, h) == g
         end

         f1, q1 = divides(f*f, f)
         @test f1 && (iszero(f) || q1 == f)

         f2, q2 = divides(f*g, f)
         @test f2 && (iszero(f) || q2 == g)

         f3, q3 = divides(f*h, f)
         @test f3 && (iszero(f) || q3 == h)

         f1, q1 = divides(g*g, g)
         @test f1 && (iszero(g) || q1 == g)

         f2, q2 = divides(g*f, g)
         @test f2 && (iszero(g) || q2 == f)

         f3, q3 = divides(g*h, g)
         @test f3 && (iszero(g) || q3 == h)

         f1, q1 = divides(h*h, h)
         @test f1 && (iszero(h) || q1 == h)

         f2, q2 = divides(h*f, h)
         @test f2 && (iszero(h) || q2 == f)

         f3, q3 = divides(h*g, h)
         @test f3 && (iszero(h) || q3 == g)
      end
   end
end

@testset "Generic.UnivPoly.euclidean_division" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = universal_polynomial_ring(R; cached=false)

         f = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         g = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         h = rand(S, 0:5, 0:10, -10:10)

         if !iszero(f)
            @test div(f, f) == 1
            @test div(f*g, f) == g
            @test div(f*h, f) == h

            q, r = divrem(f, f)
            @test q == 1 && r == 0
            q, r = divrem(f*g, f)
            @test q == g && r == 0
            q, r = divrem(f*h, f)
            @test q == h && r == 0
         end

         if !iszero(g)
            @test div(g, g) == 1
            @test div(g*f, g) == f
            @test div(g*h, g) == h

            q, r = divrem(g, g)
            @test q == 1 && r == 0
            q, r = divrem(g*f, g)
            @test q == f && r == 0
            q, r = divrem(g*h, g)
            @test q == h && r == 0
         end

         if !iszero(h)
            @test div(h, h) == 1
            @test div(h*f, h) == f
            @test div(h*g, h) == g
            
            q, r = divrem(h, h)
            @test q == 1 && r == 0
            q, r = divrem(h*f, h)
            @test q == f && r == 0
            q, r = divrem(h*g, h)
            @test q == g && r == 0
         end
      end
   end
end

@testset "Generic.UnivPoly.derivative" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = universal_polynomial_ring(R; cached=false)

         x = gen(S, "x")
         y, z = gens(S, ["y", "z"])
         x2 = gen(S, "x")

         f = 3x^3 + 2x^2 + x + 4
         g = 3x^3*y^2 + 2x^3*y*z + 2x^2*y*z + 3x + 2y + 1

         @test derivative(f, 1) == 9x^2 + 4x + 1
         @test derivative(f, 2) == 0

         @test derivative(g, 1) == 9x^2*y^2 + 6x^2*y*z + 4x*y*z + 3
         @test derivative(g, 2) == 6x^3*y + 2x^3*z + 2x^2*z + 2
         
         @test derivative(f, x) == 9x^2 + 4x + 1
         @test derivative(f, x2) == 9x^2 + 4x + 1
         @test derivative(f, y) == 0

         @test derivative(g, x) == 9x^2*y^2 + 6x^2*y*z + 4x*y*z + 3
         @test derivative(g, x2) == 9x^2*y^2 + 6x^2*y*z + 4x*y*z + 3
         @test derivative(g, y) == 6x^3*y + 2x^3*z + 2x^2*z + 2
      end
   end
end

@testset "Generic.UnivPoly.remove_valuation" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = universal_polynomial_ring(R; cached=false)

         f = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         g = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         h = rand(S, 0:5, 0:10, -10:10)

         n = rand(0:5)
         p = f*f^n
         if !iszero(f) && !is_unit(f)
            val, q = remove(p, f)
            val2 = valuation(p, f)
            @test val == val2
            @test p == q*f^val
            @test val >= n
         end

         n = rand(0:5)
         p = g*f^n
         if !iszero(f) && !is_unit(f) && !iszero(g)
            val, q = remove(p, f)
            val2 = valuation(p, f)
            @test val == val2
            @test p == q*f^val
            @test val >= n
         end

         n = rand(0:5)
         p = h*f^n
         if !iszero(f) && !is_unit(f) && !iszero(h)
            val, q = remove(p, f)
            val2 = valuation(p, f)
            @test val == val2
            @test p == q*f^val
            @test val >= n
         end

         n = rand(0:5)
         p = f*g^n
         if !iszero(f) && !iszero(g) && !is_unit(g)
            val, q = remove(p, g)
            val2 = valuation(p, g)
            @test val == val2
            @test p == q*g^val
            @test val >= n
         end

         n = rand(0:5)
         p = g*g^n
         if !iszero(g) && !is_unit(g)
            val, q = remove(p, g)
            val2 = valuation(p, g)
            @test val == val2
            @test p == q*g^val
            @test val >= n
         end

         n = rand(0:5)
         p = h*g^n
         if !iszero(g) && !is_unit(g) && !iszero(h)
            val, q = remove(p, g)
            val2 = valuation(p, g)
            @test val == val2
            @test p == q*g^val
            @test val >= n
         end

         n = rand(0:5)
         p = f*h^n
         if !iszero(f) && !iszero(h) && !is_unit(h)
            val, q = remove(p, h)
            val2 = valuation(p, h)
            @test val == val2
            @test p == q*h^val
            @test val >= n
         end

         n = rand(0:5)
         p = g*h^n
         if !iszero(g) && !is_unit(h) && !iszero(h)
            val, q = remove(p, h)
            val2 = valuation(p, h)
            @test val == val2
            @test p == q*h^val
            @test val >= n
         end

         n = rand(0:5)
         p = h*h^n
         if !iszero(h) && !is_unit(h)
            val, q = remove(p, h)
            val2 = valuation(p, h)
            @test val == val2
            @test p == q*h^val
            @test val >= n
         end
      end
   end
end

@testset "Generic.UnivPoly.evaluation" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = universal_polynomial_ring(R; cached=false)

         f = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         g = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         h = rand(S, 0:5, 0:10, -10:10)

         U, w = polynomial_ring(R, "v")

         n = rand(1:3)
         V = [rand(-10:10) for v in 1:n]

         @test evaluate(f, V) == evaluate(f, [R(v) for v in V])
         @test evaluate(f, V) == evaluate(f, [ZZ(v) for v in V])
         @test evaluate(f, V) == evaluate(f, [U(v) for v in V])
         @test evaluate(f, V) == f(V...)
         @test evaluate(f, V) == f([ZZ(v) for v in V]...)
         @test evaluate(f, V) == f([U(v) for v in V]...)

         @test evaluate(g, V) == evaluate(g, [R(v) for v in V])
         @test evaluate(g, V) == evaluate(g, [ZZ(v) for v in V])
         @test evaluate(g, V) == evaluate(g, [U(v) for v in V])
         @test evaluate(g, V) == g(V...)
         @test evaluate(g, V) == g([ZZ(v) for v in V]...)
         @test evaluate(g, V) == g([U(v) for v in V]...)

         @test evaluate(h, V) == evaluate(h, [R(v) for v in V])
         @test evaluate(h, V) == evaluate(h, [ZZ(v) for v in V])
         @test evaluate(h, V) == evaluate(h, [U(v) for v in V])
         @test evaluate(h, V) == h(V...)
         @test evaluate(h, V) == h([ZZ(v) for v in V]...)
         @test evaluate(h, V) == h([U(v) for v in V]...)

         V = [rand(-10:10) for v in 1:2]

         @test evaluate(f, [1], [V[1]]) == evaluate(f, [1], [R(V[1])])
         @test evaluate(f, [1], [V[1]]) == evaluate(f, [1], [ZZ(V[1])])
         @test evaluate(f, [1], [V[1]]) == evaluate(f, [1], [U(V[1])])
         @test evaluate(f, [1, 3], [V[1], V[2]]) == evaluate(f, [1, 3], [R(v) for v in V[1:2]])
         @test evaluate(f, [1, 3], [V[1], V[2]]) == evaluate(f, [1, 3], [ZZ(v) for v in V[1:2]])

         @test evaluate(g, [1], [V[1]]) == evaluate(g, [1], [R(V[1])])
         @test evaluate(g, [1], [V[1]]) == evaluate(g, [1], [ZZ(V[1])])
         @test evaluate(g, [1, 3], [V[1], V[2]]) == evaluate(g, [1, 3], [R(v) for v in V[1:2]])
         @test evaluate(g, [1, 3], [V[1], V[2]]) == evaluate(g, [1, 3], [ZZ(v) for v in V[1:2]])

         @test evaluate(h, [1], [V[1]]) == evaluate(h, [1], [R(V[1])])
         @test evaluate(h, [1], [V[1]]) == evaluate(h, [1], [ZZ(V[1])])
         @test evaluate(h, [1, 3], [V[1], V[2]]) == evaluate(h, [1, 3], [R(v) for v in V[1:2]])
         @test evaluate(h, [1, 3], [V[1], V[2]]) == evaluate(h, [1, 3], [ZZ(v) for v in V[1:2]])

         @test evaluate(x, [1], [y]) == evaluate(z, [3], [y])
      end
   end
end

@testset "Generic.UnivPoly.gcd" begin
   for R in [ZZ, QQ]
      for iters = 1:100
         S = universal_polynomial_ring(R; cached=false)

         f = rand(S, 0:5, 0:10, -10:10)
         x = gen(S, "x")
         g = rand(S, 0:5, 0:10, -10:10)
         y, z = gens(S, ["y", "z"])
         h = rand(S, 0:5, 0:10, -10:10)

         uf = canonical_unit(f)
         ug = canonical_unit(g)
         uh = canonical_unit(h)

         @test iszero(uf) || gcd(f, f) == divexact(f, uf)
         @test iszero(uf) || gcd(f*g, f) == divexact(f, uf)
         @test iszero(uf) || gcd(f*h, f) == divexact(f, uf)

         @test iszero(ug) || gcd(g, g) == divexact(g, ug)
         @test iszero(ug) || gcd(g*h, g) == divexact(g, ug)

         @test iszero(uh) || gcd(h, h) == divexact(h, uh)

         @test lcm(f, f) == f*uf
         @test lcm(f*g, f) == f*g*uf
         @test lcm(f*h, f) == f*h*uf

         @test lcm(g, g) == g*ug
         @test lcm(g*h, g) == g*h*ug

         @test lcm(h, h) == h*uh
      end
   end
end

@testset "Generic.UnivPoly.univariate_polynomials" begin
   for R in [ZZ, QQ]
      S = universal_polynomial_ring(R; cached=false)

      U, y = polynomial_ring(R, "y")

      x = gen(S, "x")

      @test is_univariate(S)

      y, z = gens(S, ["y", "z"])

      @test !is_univariate(S)

      f = 3x^3 + 2x^2 + x + 4
      g = 3x^3*y^2 + 2x^3*y*z + 2x^2*y*z + 3x + 2y + 1
      h = 3y^2 + 2y + 1

      @test is_univariate(f)
      @test !is_univariate(g)
      @test is_univariate(h)

      f1 = to_univariate(U, f)
      h1 = to_univariate(U, h)

      @test length(f) == length(f1)
      @test length(h) == length(h1)

      @test coefficients_of_univariate(f) == [R(4), R(1), R(2), R(3)]
      @test coefficients_of_univariate(h) == [R(1), R(2), R(3)]
   end
end

@testset "Generic.UnivPoly.map" begin
   for R in [ZZ, QQ]
      S = universal_polynomial_ring(R; cached=false)

      U, y = polynomial_ring(R, "y")

      x = gen(S, "x")
      y, z = gens(S, ["y", "z"])

      f = 3x^3 + 2x^2 + x + 4
      g = 3x^3*y^2 + 2x^3*y*z + 2x^2*y*z + 3x + 2y + 1
      h = 3y^2 + 2y + 1

      f1 = change_base_ring(U, f)
      g1 = change_base_ring(U, g)
      h1 = change_base_ring(U, h)

      @test length(f) == length(f1)
      @test length(g) == length(g1)
      @test length(h) == length(h1)

      @test base_ring(f1) === U
      @test base_ring(g1) === U
      @test base_ring(h1) === U

      f2 = map_coefficients(x->x^2, f)
      g2 = map_coefficients(x->x^2, g)
      h2 = map_coefficients(x->x^2, h)

      @test length(f) == length(f2)
      @test length(g) == length(g2)
      @test length(h) == length(h2)
   end
end

@testset "Generic.UnivPoly.unsafe_operators" begin
   for R in [ZZ, QQ]
      S = universal_polynomial_ring(R; cached=false)

      U, y = polynomial_ring(R, "y")

      x = gen(S, "x")
      y, z = gens(S, ["y", "z"])

      f = 3x^3 + 2x^2 + x + 4
      g = 3x^3*y^2 + 2x^3*y*z + 2x^2*y*z + 3x + 2y + 1
      h = 3y^2 + 2y + 1

      f1 = deepcopy(f)
      f1 = add!(f1, g, h)
      
      @test f1 == g + h

      f2 = deepcopy(f)
      f2 = add!(f2, g)

      @test f2 == f + g

      f3 = deepcopy(f)
      f3 = mul!(f3, g, h)

      @test f3 == g*h

      f4 = deepcopy(f)
      f4 = addmul!(f4, g, h)

      @test f4 == f + g*h

      f5 = deepcopy(f)
      f5 = zero!(f5)

      @test f5 == 0
   end
end

