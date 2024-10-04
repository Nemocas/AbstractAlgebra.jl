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
   if is_exact_type(T)
      return a == b
   else
      return isapprox(a, b)
   end
end

# The following functions should not expect that their input is a `NCRingElem` or similar.
# They should be usable in more general types, that don't even have a `parent/elem` correspondence
function test_mutating_op_like_zero(f::Function, f!::Function, A)
   a = deepcopy(A)
   a = f!(a)
   @test equality(a, f(A))
end

function test_mutating_op_like_neg(f::Function, f!::Function, A)
   # initialize storage var with different values to check that its value is not used
   for z in [zero(A), deepcopy(A)]
      a = deepcopy(A)
      z = f!(z, a)
      @test equality(z, f(A))
      @test a == A
   end

   a = deepcopy(A)
   a = f!(a)
   @test equality(a, f(A))
end

function test_mutating_op_like_add(f::Function, f!::Function, A, B)
   # initialize storage var with different values to check that its value is not used
   for z in [zero(A), deepcopy(A), deepcopy(B)]
      a = deepcopy(A)
      b = deepcopy(B)
      z = f!(z, a, b)
      @test equality(z, f(A, B))
      @test a == A
      @test b == B
   end

   a = deepcopy(A)
   b = deepcopy(B)
   a = f!(a, a, b)
   @test equality(a, f(A, B))
   @test b == B

   a = deepcopy(A)
   b = deepcopy(B)
   b = f!(b, a, b)
   @test equality(b, f(A, B))
   @test a == A

   a = deepcopy(A)
   b = deepcopy(B)
   a = f!(a, b, b)
   @test equality(a, f(B, B))
   @test b == B
   
   b = deepcopy(B)
   b = f!(b, b, b)
   @test equality(b, f(B, B))

   a = deepcopy(A)
   b = deepcopy(B)
   a = f!(a, b)
   @test equality(a, f(A, B))
   @test b == B

   b = deepcopy(B)
   b = f!(b, b)
   @test equality(b, f(B, B))
end

function test_mutating_op_like_addmul(f::Function, f!_::Function, Z, A, B)
   f!(z, a, b, ::Nothing) = f!_(z, a, b)
   f!(z, a, b, t) = f!_(z, a, b, t)

   # initialize storage var with different values to check that its value is not used
   # and `nothing` for the three-arg dispatch
   for t in [nothing, zero(A), deepcopy(A)]
      z = deepcopy(Z)
      a = deepcopy(A)
      b = deepcopy(B)
      z = f!(z, a, b, t)
      @test equality(z, f(Z, A, B))
      @test a == A
      @test b == B
      
      a = deepcopy(A)
      b = deepcopy(B)
      a = f!(a, a, b, t)
      @test equality(a, f(A, A, B))
      @test b == B
      
      a = deepcopy(A)
      b = deepcopy(B)
      b = f!(b, a, b, t)
      @test equality(b, f(B, A, B))
      @test a == A
      
      a = deepcopy(A)
      b = deepcopy(B)
      a = f!(a, b, b, t)
      @test equality(a, f(A, B, B))
      @test b == B
      
      b = deepcopy(B)
      b = f!(b, b, b, t)
      @test equality(b, f(B, B, B))
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
         @test is_domain_type(T) isa Bool
         @test is_exact_type(T) isa Bool

         @test base_ring_type(R) == typeof(base_ring(R))
         @test base_ring_type(zero(R)) == typeof(base_ring(zero(R)))
         @test base_ring_type(typeof(R)) == typeof(base_ring(R))
         @test base_ring_type(T) == typeof(base_ring(zero(R)))

         # some rings don't support characteristic and raise an exception (see issue #993)
         try ch = characteristic(R)
            @test iszero(R(characteristic(R)))
            @test iszero(characteristic(R) * one(R))
            @test iszero(one(R) * characteristic(R))
         catch
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
         @test isone(R(0)) || !is_unit(R(0))
         @test is_unit(R(1))
         for i in 1:reps
            a = test_elem(R)::T
            @test hash(a) isa UInt
            A = deepcopy(a)
            @test !ismutable(a) || a !== A
            @test equality(a, A)
            @test hash(a) == hash(A)
            @test parent(a) === parent(A)
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
            @test equality(+(a), a)
            @test equality(*(a), a)
            @test equality(a^1, a)
            @test equality(a^2, a*a)
            @test equality(a^3, a*a*a)
            @test equality(a^4, a*a*a*a)
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
               if is_domain_type(T)
                  @test iszero(b) || equality(divexact_left(b*a, b), a)
                  @test iszero(b) || equality(divexact_left(b*a, b, check = true), a)
                  @test iszero(b) || equality(divexact_left(b*a, b, check = false), a)
                  @test iszero(b) || equality(b \ (b*a), a)

                  @test iszero(b) || equality(divexact_right(a*b, b), a)
                  @test iszero(b) || equality(divexact_right(a*b, b, check = true), a)
                  @test iszero(b) || equality(divexact_right(a*b, b, check = false), a)
                  @test iszero(b) || equality((a*b) / b, a)
               else
                  try
                     t = divexact_left(b*a, b)
                     @test equality(b*t, b*a)
                     t = divexact_left(b*a, b, check = true)
                     @test equality(b*t, b*a)
                     t = divexact_left(b*a, b, check = false)
                     @test equality(b*t, b*a)
                     t = b \ (b*a)
                     @test equality(b*t, b*a)
                  catch
                  end
                  try
                     t = divexact_right(a*b, b)
                     @test equality(t*b, a*b)
                     t = divexact_right(a*b, b, check = true)
                     @test equality(t*b, a*b)
                     t = divexact_right(a*b, b, check = false)
                     @test equality(t*b, a*b)
                     t = (a*b) / b
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
            
            test_mutating_op_like_zero(zero, zero!, a)
            test_mutating_op_like_zero(one, one!, a)

            test_mutating_op_like_neg(-, neg!, a)

            test_mutating_op_like_add(+, add!, a, b)
            test_mutating_op_like_add(-, sub!, a, b)
            test_mutating_op_like_add(*, mul!, a, b)

            test_mutating_op_like_addmul((a, b, c) -> a + b*c, addmul!, a, b, c)
            test_mutating_op_like_addmul((a, b, c) -> a - b*c, submul!, a, b, c)
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
         @test isone(AbstractAlgebra.inv(one(R)))
         test_mutating_op_like_neg(AbstractAlgebra.inv, inv!, one(R))
         test_mutating_op_like_neg(AbstractAlgebra.inv, inv!, -one(R))
         for i in 1:reps
            a = test_elem(R)::T
            b = test_elem(R)::T
            A = deepcopy(a)
            B = deepcopy(b)
            @test a*b == b*a
            # documentation is not clear on divexact
            if is_domain_type(T)
               @test iszero(b) || equality(divexact(b*a, b), a)
               @test iszero(b) || equality(divexact(b*a, b, check = true), a)
               @test iszero(b) || equality(divexact(b*a, b, check = false), a)
               if T isa RingElem
                  @test iszero(b) || equality((b*a) / b, a)
               end
               iszero(b) || test_mutating_op_like_add(divexact, divexact!, b*a, b)
            else
               try
                  t = divexact(b*a, b)
                  @test equality(t*b, a*b)
                  t = divexact(b*a, b, check = true)
                  @test equality(t*b, a*b)
                  t = divexact(b*a, b, check = false)
                  @test equality(t*b, a*b)
                  if T isa RingElem
                     t = (b*a) / b
                     @test equality(t*b, a*b)
                  end
               catch
               end
            end
            try
               (f, h) = is_zero_divisor_with_annihilator(a)
               @test parent(h) == R
               @test f == is_zero_divisor(a)
               if f
                  @test !is_zero(h)
                  @test is_zero(a*h)
               end
            catch
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
         A = deepcopy(a)
         @test is_unit(a) == !iszero(a)
         if !is_zero(a)
            @test is_one(a * inv(a))
            @test is_one(inv(a) * a)
            test_mutating_op_like_neg(inv, inv!, a)
         end
         @test A == a
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

   is_exact_type(T) || return

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

         if is_unit(gcd(f, m))
            a = invmod(f, m)
            @test divides(mulmod(a, f, m) - one(R), m)[1]
            @test divides(powermod(f, -1, m) - a^1, m)[1]
            @test divides(powermod(f, -2, m) - a^2, m)[1]
            @test divides(powermod(f, -3, m) - a^3, m)[1]
         end

         @test divides(f*m, m) == (true, f)
         (a, b) = divides(f*m + g, m)
         @test !a || b*m == f*m + g

         @test_throws Exception remove(f, zero(R))
         @test_throws Exception valuation(f, zero(R))

         if !is_unit(m) && !iszero(f)
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

         g1 = gcd(f, gcd(g, m)) 
         g2 = gcd(gcd(f, g), m)
         g3 = gcd(f, g, m)
         g4 = gcd([f, g, m])
         @test equality_up_to_units(g1, g2)
         @test equality_up_to_units(g2, g3)
         @test equality_up_to_units(g3, g4)

         l1 = lcm(f, lcm(g, m))
         l2 = lcm(lcm(f, g), m)
         l3 = lcm(f, g, m)
         l4 = lcm([f, g, m])
         @test equality_up_to_units(l1, l2)
         @test equality_up_to_units(l2, l3)
         @test equality_up_to_units(l3, l4)

         (d, s, t) = gcdx(f, g)
         @test d == gcd(f, g)
         @test d == s*f + t*g
         @test gcdinv(f, g) == (d, s)

         test_mutating_op_like_add(AbstractAlgebra.div, div!, f, m)
         test_mutating_op_like_add(mod, mod!, f, m)
         test_mutating_op_like_add(gcd, gcd!, f, m)
         test_mutating_op_like_add(lcm, lcm!, f, m)
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

            B = MPolyBuildCtx(Rx)   # TODO rename to BuildCtx
            for (c, e) in zip(AbstractAlgebra.coefficients(a), AbstractAlgebra.exponent_vectors(a))
               push_term!(B, c, e)
            end
            @test finish(B) == a
         end
         @test Rx(Int[]) == zero(Rx)
         @test Rx([0, 1, 2]) == x + 2*x^2
         @test Rx([big(0), big(1), big(2)]) == x + 2*x^2
         @test Rx(map(R, [0, 1, 2])) == x + 2*x^2
      end

      if R isa AbstractAlgebra.Field
         test_EuclideanRing_interface(Rx, reps = 2 + fld(reps, 2))
         @testset "Half-GCD" begin
            for i in 1:reps
               a = test_elem(Rx)
               b = test_elem(Rx)
               for j in 1:8
                  q = test_elem(Rx)
                  a, b = q*a + b, a
               end
               g, s, t = gcdx(a, b)
               @test g == gcd(a, b)
               @test g == s*a + t*b
               @test (g, s) == gcdinv(a, b)
               if degree(a) < degree(b)
                  a, b = b, a
               end
               degree(a) > degree(b) >= 0 || continue
               (A, B, m11, m12, m21, m22, s) = hgcd(a, b)
               @test degree(A) >= cld(degree(a), 2) > degree(B)
               @test m11*A + m12*B == a
               @test m21*A + m22*B == b
               @test m11*m22 - m21*m12 == s
               @test s^2 == 1
            end
         end
      end

      @testset "Basic functionality" begin
         @test var(Rx) isa Symbol
         @test symbols(Rx) isa Vector{Symbol}
         @test length(symbols(Rx)) == 1
         @test is_gen(gen(Rx))
         @test is_monic(x)
         for i in 1:reps
            a = test_elem(Rx)
            @test iszero(a) || degree(a) >= 0
            @test equality(a, leading_coefficient(a)*x^max(0, degree(a)) + tail(a))
            @test constant_coefficient(a) isa elem_type(R)
            @test trailing_coefficient(a) isa elem_type(R)
            @test is_gen(x)
            @test iszero(one(Rx)) || !is_gen(x^2)
            @test is_monic(a) == isone(leading_coefficient(a))
         end
      end
   end

   return nothing
end

function test_elem(S::Union{AbstractAlgebra.MatSpace,
                            AbstractAlgebra.MatRing})
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

function test_MatAlgebra_interface(S::MatRing; reps = 20)

   ST = elem_type(S)
   R = base_ring(S)
   T = elem_type(R)

   @test nrows(S) == ncols(S)

   @testset "MatRing interface for $(S) of type $(typeof(S))" begin

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
   Rx, _ = polynomial_ring(R, "x")
   test_Poly_interface(Rx, reps = 2 + fld(reps, 2))
   S = matrix_ring(R, rand(0:3))
   test_MatAlgebra_interface(S, reps = 2 + fld(reps, 2))
   S = matrix_space(R, rand(0:3), rand(0:3))
   test_MatSpace_interface(S, reps = 2 + fld(reps, 2))
end

function test_Field_interface_recursive(R::AbstractAlgebra.Field; reps = 50)
   test_Ring_interface_recursive(R, reps = reps)
   test_Field_interface(R, reps = reps)
end
