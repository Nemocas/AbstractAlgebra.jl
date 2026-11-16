@testset "Generic.RationalFunctionFieldHom.canonical" begin
  F, t = rational_function_field(QQ, :t)
  E, u = rational_function_field(QQ, :u)
  a = (2u + 3)//(u + 4)
  f = hom(F, E, a)
  @test domain(f) === F
  @test codomain(f) === E
  @test f(t) == a
  @test f(0) == zero(E)
  @test f(1) == one(E)
  @test f(2//3) == E(2//3)
  @test image(f, t) == a
  @test image_fn(f)(t) == a
  x = (t^3 + 2t + 1)//(t^2 - 2)
  @test f(x) == (a^3 + 2a + 1)//(a^2 - 2)
  @test parent(f(x)) === E
  @test Generic.morphism_type(F, E) === typeof(f)
  @test Generic.morphism_type(typeof(F), typeof(E)) === typeof(f)
  @test sprint(show, f) isa String
  @test sprint(show, "text/plain", f) isa String

  for sigma in (nothing, identity, identity_map(QQ))
    g = hom(F, E, sigma, a)
    @test typeof(g) === typeof(f)
    @test Generic.morphism_type(F, E, sigma) === typeof(g)
    @test Generic._coefficient_map(g) === nothing
    @test g(x) == f(x)
  end

  @test_throws ArgumentError hom(F, E, 1)
  @test_throws ArgumentError hom(F, E, identity_map(GF(7)), u)
  @test_throws ArgumentError hom(F, QQ, QQ(1))
  @test_throws ArgumentError hom(F, GF(7), 1)

  # All constructor validation can be disabled independently of coercion.
  wrong_domain = map_from_func(GF(7), E, c -> E(c))
  @test_throws ArgumentError hom(F, E, wrong_domain, u)
  @test Generic._coefficient_map(hom(F, E, wrong_domain, u; check=false)) === wrong_domain
  ring_map = map_from_func(QQ, ZZ, c -> ZZ(c))
  @test_throws ArgumentError hom(F, E, ring_map, u)
  @test Generic._coefficient_map(hom(F, E, ring_map, u; check=false)) === ring_map
  k = GF(7)
  @test codomain(hom(F, k, one(k); check=false)) === k

  F1, t1 = rational_function_field(QQ, :t; cached=false)
  @test_throws ArgumentError f(t1)
  E1, u1 = rational_function_field(QQ, :u; cached=false)
  @test_throws ArgumentError hom(F, E, u1)
  M, (v, w) = rational_function_field(QQ, [:v, :w])
  @test_throws MethodError hom(M, E, u)
  @test_throws MethodError hom(F, ZZ, 1)
end

@testset "Generic.RationalFunctionFieldHom.coefficients" begin
  K, s = rational_function_field(QQ, :s)
  F, t = rational_function_field(K, :t)
  L, v = rational_function_field(QQ, :v)
  E, u = rational_function_field(L, :u)
  sigma = hom(K, E, u)
  f = hom(F, E, sigma, E(v))
  x = (F(s)*t + 1)//(t - F(s))
  expected = (u*E(v) + 1)//(E(v) - u)
  @test f(F(s)) == u
  @test f(t) == E(v)
  @test f(x) == expected
  @test f(x + t) == f(x) + f(t)
  @test f(x*t) == f(x)*f(t)
  @test f(inv(x)) == inv(f(x))
  @test Generic.morphism_type(F, E, sigma) === typeof(f)

  sigma_func = c -> sigma(c)
  g = hom(F, E, sigma_func, E(v))
  @test g(x) == expected
  @test Generic.morphism_type(F, E, sigma_func) === typeof(g)

  # Coefficient maps into the target's coefficient field also work.
  tau = hom(K, L, v + 1)
  h = hom(F, E, tau, u^2)
  @test h(F(s)) == E(v + 1)
  @test h(t) == u^2
  @test h(x) == (E(v + 1)*u^2 + 1)//(u^2 - E(v + 1))

  # A constant outer image can be transcendental over a smaller coefficient image.
  Qx, z = rational_function_field(QQ, :z)
  h = hom(Qx, E, E(v))
  @test h((z + 1)//(z^2 + 2)) == E((v + 1)//(v^2 + 2))
end

@testset "Generic.RationalFunctionFieldHom.arbitrary_codomain" begin
  F, t = rational_function_field(QQ, :t)
  K, x = rational_function_field(QQ, :x)
  R, z = polynomial_ring(K, :z)
  E, y = function_field(z^2 - x, :y)
  f = hom(F, E, y)
  @test codomain(f) === E
  @test f(t) == y
  @test f(t^2) == E(x)
  @test f((t + 1)//(t^2 + 2)) == (y + 1)//(E(x) + 2)
  @test f(t^30 + t + 1) == y^30 + y + 1
  @test iszero(f(zero(F)))
  @test parent(f(zero(F))) === E
  @test_throws ArgumentError inv(f)

  # Multivariate targets are fields too; only the source is restricted.
  M, (u, v) = rational_function_field(QQ, [:u, :v])
  g = hom(F, M, u + v)
  @test g((t + 1)//(t - 1)) == (u + v + 1)//(u + v - 1)

  k = GF(3)
  Fp, s = rational_function_field(k, :s)
  Ep, w = rational_function_field(k, :w)
  h = hom(Fp, Ep, w^3)
  @test h((s + 1)//s) == (w^3 + 1)//w^3
  @test_throws ArgumentError inv(h)
  @test_throws ArgumentError hom(Fp, k, one(k))
  @test codomain(hom(Fp, k, one(k); check=false)) === k
end

@testset "Generic.RationalFunctionFieldHom.composition" begin
  F, t = rational_function_field(QQ, :t)
  E, u = rational_function_field(QQ, :u)
  H, v = rational_function_field(QQ, :v)
  f = hom(F, E, u^2)
  g = hom(E, H, v + 1)
  h = f*g
  @test domain(h) === F
  @test codomain(h) === H
  @test h(t) == (v + 1)^2
  @test Generic._coefficient_map(h) === nothing
  @test h((t + 1)//(t - 1)) == g(f((t + 1)//(t - 1)))
  @test identity_map(F)*f === f
  @test f*identity_map(E) === f
  @test_throws ErrorException f*hom(F, H, v)

  K, s = rational_function_field(QQ, :s)
  F, t = rational_function_field(K, :t)
  L, v = rational_function_field(QQ, :v)
  E, u = rational_function_field(L, :u)
  M, w = rational_function_field(QQ, :w)
  H, z = rational_function_field(M, :z)
  f = hom(F, E, hom(K, E, u; check=false), E(v); check=false)
  g = hom(E, H, hom(L, H, z; check=false), H(w); check=false)
  h = f*g
  @test h(F(s)) == H(w)
  @test h(t) == z
  a = (F(s)*t + 1)//(t - F(s))
  @test h(a) == g(f(a))
end

@testset "Generic.RationalFunctionFieldHom.inverse" begin
  for k in (QQ, GF(2), GF(7))
    F, t = rational_function_field(k, :t)
    images = [t + 1, inv(t), (t + 1)//t]
    for a in images
      f = hom(F, F, a)
      g = inv(f)
      @test domain(g) === F
      @test codomain(g) === F
      x = (t^2 + t + 1)//(t + 1)
      y = (t^3 + 1)//(t^2 + t + 1)
      @test g(f(x)) == x
      @test f(g(y)) == y
      @test preimage(f, f(x)) == x
      @test f \ f(x) == x
    end
  end

  F, t = rational_function_field(QQ, :t)
  f = hom(F, F, (2t + 3)//(t + 4))
  @test inv(f)(t) == (4t - 3)//(2 - t)
  # The underlying fraction interface also permits unreduced representatives.
  Q = Generic.underlying_fraction_field(F)
  R = base_ring(Q)
  z = gen(R)
  unreduced = F(Q((z + 1)*(2z + 3), (z + 1)*(z + 4); reduce=false))
  @test inv(hom(F, F, unreduced))(t) == inv(f)(t)
  @test_throws ArgumentError hom(F, F, F(Q(z + 1, z + 1; reduce=false)))
  @test_throws ArgumentError inv(hom(F, F, t^2))
  @test_throws ArgumentError inv(hom(F, F, F(2); check=false))
  E, u = rational_function_field(QQ, :u)
  @test_throws ArgumentError preimage(f, u)
  @test_throws ArgumentError inv(hom(F, E, u))
  E, u = rational_function_field(QQ, :t; cached=false)
  @test_throws ArgumentError inv(hom(F, E, u))

  # The inverse coefficient map must act on the Möbius coefficients too.
  K, s = rational_function_field(QQ, :s)
  F, t = rational_function_field(K, :t)
  sigma = hom(K, K, s + 1)
  a = (F(s)*t + F(s + 1))//(t + F(s + 2))
  f = hom(F, F, sigma, a)
  g = inv(f)
  x = (F(s)*t + 1)//(t + F(s + 1))
  y = (F(s)*t + 2)//(t + F(s + 1))
  @test g(f(x)) == x
  @test f(g(y)) == y
  @test g(F(s)) == F(s - 1)
  @test inv(f*g)(x) == x
  canonical = hom(F, F, t + 1)
  @test inv(canonical*f)((canonical*f)(x)) == x
  @test inv(f*canonical)((f*canonical)(x)) == x

  sigma_with_inverse = map_from_func(K, K, c -> sigma(c), c -> inv(sigma)(c))
  @test inv(hom(F, F, sigma_with_inverse, a))(f(x)) == x
  sigma_without_inverse = map_from_func(K, K, c -> sigma(c))
  @test_throws ArgumentError inv(hom(F, F, sigma_without_inverse, a))

  # Degree one does not compensate for a nonsurjective coefficient map.
  @test_throws ArgumentError inv(hom(F, F, hom(K, K, s^2), t))
end
