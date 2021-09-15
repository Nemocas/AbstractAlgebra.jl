# very generic testing: just define test_elem(R) to produce elements of R,
# then invoke one of these functions, as appropriate:
# - test_NCRing_interface(R)
# - test_Ring_interface(R)
# - test_Ring_interface_recursive(R)
# - test_Field_interface(R)
# - test_Field_interface_recursive(R)
#
# The "recursive" variants perform additional tests on algebraic
# structures derived from the original ring, by calling these helpers:
# - test_EuclideanRing_interface(R)
# - test_Poly_interface(R)
# - test_MatSpace_interface(R)
# - test_MatAlgebra_interface(R)

function equality(a::T, b::T) where T <: AbstractAlgebra.NCRingElement
   if isexact_type(T)
      return a == b
   else
      return isapprox(a, b)
   end
end


function test_NCRing_interface(R::AbstractAlgebra.NCRing; reps = 50)

   T = elem_type(R)

   @testset "NCRing interface for $(R) of type $(typeof(R))" begin

      @test T <: NCRingElem || T <: RingElement

      @testset "Functions for types and parents of rings" begin
         @test elem_type(typeof(R)) == T
         @test parent_type(T) == typeof(R)
         for i in 1:reps
            a = test_elem(R)::T
            @test parent(a) == R
         end
         @test isdomain_type(T) isa Bool
         @test isexact_type(T) isa Bool

         # some rings don't support characteristic and return -1 (see issue #993)
         if characteristic(R) != -1
            @test iszero(R(characteristic(R)))
            @test iszero(characteristic(R) * one(R))
            @test iszero(one(R) * characteristic(R))
         end
      end

      @testset "Constructors" begin
         @test R() isa T
         @test R(true) isa T
         @test R(false) isa T
         @test R(0) isa T
         @test R(1) isa T
         @test R(-2) isa T
         @test R(BigInt(0)) isa T
         @test R(BigInt(1)) isa T
         @test R(BigInt(-2)) isa T
         @test R(BigInt(3)^100) isa T
         for i in 1:reps
            a = test_elem(R)::T
            @test R(a) isa T
         end
      end

      @testset "Basic functions" begin
         @test iszero(R())       # R() is supposed to construct 0 ?
         @test iszero(zero(R))
         @test isone(one(R))
         @test iszero(R(0))
         @test isone(R(1))
         @test isone(R(0)) || !isunit(R(0))
         @test isunit(R(1))
         for i in 1:reps
            a = test_elem(R)::T
            @test hash(a) isa UInt
            @test hash(a) == hash(deepcopy(a))
            @test sprint(show, "text/plain", a) isa String
         end
         @test sprint(show, "text/plain", R) isa String

         for i in 1:reps
            a = test_elem(R)::T
            b = test_elem(R)::T
            c = test_elem(R)::T
            A = deepcopy(a)
            B = deepcopy(b)
            C = deepcopy(c)
            @test equality((a + b) + c, a + (b + c))
            @test equality(a + b, b + a)
            @test equality(a - c, a + (-c))
            @test equality(a + zero(R), a)
            @test equality(a + (-a), zero(R))
            @test equality((a*b)*c, a*(b*c))
            @test equality(a*one(R), a)
            @test equality(one(R)*a, a)
            @test equality(a*(b + c), a*b + a*c)
            @test equality((a + b)*c, a*c + b*c)
            @test iszero(a*zero(R))
            @test iszero(zero(R)*a)
            @test A == a
            @test B == b
            @test C == c
         end
      end

      if !(R isa AbstractAlgebra.Ring)
         @testset "Basic functionality for noncommutative rings only" begin
            for i in 1:reps
               a = test_elem(R)::T
               b = test_elem(R)::T
               A = deepcopy(a)
               B = deepcopy(b)
               # documentation is not clear on divexact
               if isdomain_type(T)
                  @test iszero(b) || equality(divexact_left(b*a, b), a)
                  @test iszero(b) || equality(divexact_right(a*b, b), a)
               else
                  try
                     t = divexact_left(b*a, b)
                     @test equality(b*t, b*a)
                  catch
                  end
                  try
                     t = divexact_right(a*b, b)
                     @test equality(t*b, a*b)
                  catch
                  end
               end
               @test A == a
               @test B == b
            end
         end
      end

      @testset "Unsafe ring operators" begin
         for i in 1:reps
            a = test_elem(R)::T
            b = test_elem(R)::T
            c = test_elem(R)::T
            A = deepcopy(a)
            B = deepcopy(b)
            C = deepcopy(c)

            x = deepcopy(a)
            @test iszero(zero!(x))

            ab = a*b
            x = R()
            @test mul!(x, a, b) == ab
            x = deepcopy(b)
            @test mul!(x, a, b) == ab
            x = deepcopy(a)
            @test mul!(x, x, b) == ab
            x = deepcopy(b)
            @test mul!(x, a, x) == ab

            ab = a+b
            x = R()
            @test add!(x, a, b) == ab
            x = deepcopy(b)
            @test add!(x, a, b) == ab
            x = deepcopy(a)
            @test add!(x, x, b) == ab
            x = deepcopy(b)
            @test add!(x, a, x) == ab
            x = deepcopy(a)
            @test addeq!(x, b) == ab

            # isapprox as BigFloat may fuse
            # matrices don't implement addmul!
            #t = R()
            #x = deepcopy(a)
            #@test equality(addmul!(x, b, c, t), a+b*c)
            #x = deepcopy(a)
            #@test equality(addmul!(x, x, c, t), a+a*c)
            #x = deepcopy(a)
            #@test equality(addmul!(x, b, x, t), a+b*a)
            #x = deepcopy(a)
            #@test equality(addmul!(x, x, x, t), a+a*a)

            @test A == a
            @test B == b
            @test C == c
         end
      end
   end

   return nothing
end


function test_Ring_interface(R::AbstractAlgebra.Ring; reps = 50)

   T = elem_type(R)

   @testset "Ring interface for $(R) of type $(typeof(R))" begin

      @test T <: RingElement

      test_NCRing_interface(R)

      @testset "Basic functionality for commutative rings only" begin
         for i in 1:reps
            a = test_elem(R)::T
            b = test_elem(R)::T
            A = deepcopy(a)
            B = deepcopy(b)
            @test isone(AbstractAlgebra.inv(one(R)))
            @test a*b == b*a
            @test equality(a^1, a)
            @test equality(a^2, a*a)
            @test equality(a^3, a*a*a)
            @test equality(a^4, a*a*a*a)
            # documentation is not clear on divexact
            if isdomain_type(T)
               @test iszero(b) || equality(divexact(b*a, b), a)
            else
               try
                  t = divexact(b*a, b)
                  @test equality(t*b, a*b)
               catch
               end
            end
            @test A == a
            @test B == b
         end
      end
   end

   return nothing
end

function test_Field_interface(R::AbstractAlgebra.Field; reps = 50)

   T = elem_type(R)

   @testset "Field interface for $(R) of type $(typeof(R))" begin

      test_Ring_interface(R, reps = reps)

      @test iszero(R(characteristic(R)))
      @test iszero(characteristic(R) * one(R))
      @test iszero(one(R) * characteristic(R))

      for i in 1:reps
         a = test_elem(R)::T
         @test isunit(a) == !iszero(a)
      end
   end

   return nothing
end

function equality_up_to_units(a, b)
   iszero(a) && return iszero(b)
   iszero(b) && return iszero(a)
   return divides(a, b)[1] && divides(b, a)[1]
end

function test_EuclideanRing_interface(R::AbstractAlgebra.Ring; reps = 20)

   T = elem_type(R)

   isexact_type(T) || return

   @testset "Euclidean Ring interface for $(R) of type $(typeof(R))" begin

      for i in 1:reps
         f = test_elem(R)::T
         g = test_elem(R)::T
         m = test_elem(R)::T
         if iszero(m)
            m = one(R)
         end

         @test (AbstractAlgebra.div(f, m), mod(f, m)) == AbstractAlgebra.divrem(f, m)
         @test divides(mulmod(f, g, m) - mod(f*g, m), m)[1]

         fi = one(R)
         for i in 1:5
            fi *= f
            @test divides(fi - powermod(f, i, m), m)[1]
            @test divides(fi - mod(f^i, m), m)[1]
         end

         if isunit(gcd(f, m))
            a = invmod(f, m)
            @test divides(mulmod(a, f, m) - one(R), m)[1]
            @test divides(powermod(f, -1, m) - a^1, m)[1]
            @test divides(powermod(f, -2, m) - a^2, m)[1]
            @test divides(powermod(f, -3, m) - a^3, m)[1]
         end

         @test divides(f*m, m) == (true, f)
         (a, b) = divides(f*m + g, m)
         @test !a || b*m == f*m + g

         #TODO: enable once https://github.com/Nemocas/Nemo.jl/issues/1174
         #@test_throws Exception remove(f, zero(R))
         #@test_throws Exception valuation(f, zero(R))

         if !isunit(m) && !iszero(f)
            n = rand(0:3)
            f *= m^n
            (v, q) = remove(f, m)
            @test valuation(f, m) == v
            @test v >= n
            @test q*m^v == f
            @test remove(q, m) == (0, q)
            @test valuation(q, m) == 0
         end

         @test !(iszero(f) && iszero(g)) || iszero(gcd(f, g))
         @test equality_up_to_units(gcd(f, g)*lcm(f, g), f*g)

         (d, s, t) = gcdx(f, g)
         @test d == gcd(f, g)
         @test d == s*f + t*g
         @test gcdinv(f, g) == (d, s)
      end

   end

   return nothing
end


function test_elem(Rx::AbstractAlgebra.PolyRing)
   R = base_ring(Rx)
   return Rx(elem_type(R)[test_elem(R) for i in 1:rand(0:6)])
end

function test_Poly_interface(Rx::AbstractAlgebra.PolyRing; reps = 30)

   T = elem_type(Rx)

   @testset "Poly interface for $(Rx) of type $(typeof(Rx))" begin

      test_Ring_interface(Rx; reps = reps)

      x = gen(Rx)
      R = base_ring(Rx)

      @testset "Polynomial Constructors" begin
         for i in 1:reps
            a = test_elem(Rx)::T
            for b in coefficients(a)
               @assert Rx(b) isa T
            end
            @test a == Rx(collect(coefficients(a)))
         end
         @test Rx(Int[]) == zero(Rx)
         @test Rx([0, 1, 2]) == x + 2*x^2
         @test Rx([big(0), big(1), big(2)]) == x + 2*x^2
         @test Rx(map(R, [0, 1, 2])) == x + 2*x^2
      end

      if R isa AbstractAlgebra.Field
         test_EuclideanRing_interface(Rx, reps = 2 + fld(reps, 2))
      end

      @testset "Basic functionality" begin
         @test var(Rx) isa Symbol
         @test symbols(Rx) isa Vector{Symbol}
         @test length(symbols(Rx)) == 1
         @test isgen(gen(Rx))
         @test ismonic(x)
         for i in 1:reps
            a = test_elem(Rx)
            @test iszero(a) || degree(a) >= 0
            @test equality(a, leading_coefficient(a)*x^max(0, degree(a)) + tail(a))
            @test constant_coefficient(a) isa elem_type(R)
            @test trailing_coefficient(a) isa elem_type(R)
            @test isgen(x)
            @test iszero(one(Rx)) || !isgen(x^2)
            @test ismonic(a) == isone(leading_coefficient(a))
         end
      end
   end

   return nothing
end

function test_elem(S::Union{AbstractAlgebra.MatSpace,
                            AbstractAlgebra.MatAlgebra})
   R = base_ring(S)
   return S(elem_type(R)[test_elem(R) for i in 1:nrows(S), j in 1:ncols(S)])
end


function test_MatSpace_interface(S::MatSpace; reps = 20)

   ST = elem_type(S)
   R = base_ring(S)
   T = elem_type(R)

   @testset "MatSpace interface for $(S) of type $(typeof(S))" begin

      @testset "Constructors" begin
         for k in 1:reps
            a = test_elem(S)::ST
            @test nrows(a) == nrows(S)
            @test ncols(a) == ncols(S)
            @test a == S(T[a[i, j] for i in 1:nrows(a), j in 1:ncols(a)])
            @test a == S(T[a[i, j] for i in 1:nrows(a) for j in 1:ncols(a)])
            @test a == matrix(R, T[a[i, j] for i in 1:nrows(a), j in 1:ncols(a)])
            @test a == matrix(R, nrows(S), ncols(S),
                              T[a[i, j] for i in 1:nrows(a) for j in 1:ncols(a)])
         end
         @test iszero(zero_matrix(R, nrows(S), ncols(S)))
      end

      @testset "Views" begin
         M = matrix(R, 3, 3, BigInt[1, 2, 3, 2, 3, 4, 3, 4, 5])
         N1 = @view M[1:2, :]
         N2 = @view M[:, 1:2]
         @test N1*N2 == matrix(R, 2, 2, BigInt[14, 20, 20, 29])
      end

      @testset "Basic manipulation of matrices" begin
         for k in 1:reps
            a = test_elem(S)::ST
            A = deepcopy(a)
            b = zero_matrix(R, nrows(a), ncols(a))
            for i in 1:nrows(a), j in 1:ncols(a)
               b[i, j] = a[i, j]
            end
            @test b == a
            @test transpose(transpose(a)) == a
            @test a == A
         end
      end
   end

   return nothing
end

function test_MatAlgebra_interface(S::MatAlgebra; reps = 20)

   ST = elem_type(S)
   R = base_ring(S)
   T = elem_type(R)

   @test nrows(S) == ncols(S)

   @testset "MatAlgebra interface for $(S) of type $(typeof(S))" begin

      test_NCRing_interface(S, reps = reps)

      @testset "Constructors" begin
         for k in 1:reps
            a = test_elem(S)::ST
            @test nrows(a) == nrows(S)
            @test ncols(a) == ncols(S)
            @test a == S(T[a[i, j] for i in 1:nrows(a), j in 1:ncols(a)])
            @test a == S(T[a[i, j] for i in 1:nrows(a) for j in 1:ncols(a)])
         end
      end

      @testset "Basic manipulation of matrices" begin
         for k in 1:reps
            a = test_elem(S)::ST
            A = deepcopy(a)
            b = zero(S)
            for i in 1:nrows(a), j in 1:ncols(a)
               b[i, j] = a[i, j]
            end
            @test b == a
            @test transpose(transpose(a)) == a
            @test a == A
         end
      end

      @testset "Determinant" begin
         for k in 1:reps
            a = test_elem(S)::ST
            b = test_elem(S)::ST
            A = deepcopy(a)
            B = deepcopy(b)
            @test det(a*b) == det(a)*det(b)
            @test a == A
            @test b == B
         end
      end
   end

   return nothing
end

function test_Ring_interface_recursive(R::AbstractAlgebra.Ring; reps = 50)
   test_Ring_interface(R; reps = reps)
   Rx, _ = PolynomialRing(R, "x")
   test_Poly_interface(Rx, reps = 2 + fld(reps, 2))
   S = MatrixAlgebra(R, rand(0:3))
   test_MatAlgebra_interface(S, reps = 2 + fld(reps, 2))
   S = MatrixSpace(R, rand(0:3), rand(0:3))
   test_MatSpace_interface(S, reps = 2 + fld(reps, 2))
end

function test_Field_interface_recursive(R::AbstractAlgebra.Field; reps = 50)
   test_Ring_interface_recursive(R, reps = reps)
   test_Field_interface(R, reps = reps)
end
