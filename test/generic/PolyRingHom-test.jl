@testset "PolyRingHom" begin
  GF2 = GF(2)
  Qx, x = QQ["x"]
  f = hom(Qx, Qx, x)
  @test sprint(show, f) isa String
  @test sprint(show, "text/plain", f) isa String
  @test f(x) == x
  @test f(x^2) == x^2
  @test f(1) == 1
  f = hom(Qx, Qx, 1)
  @test AbstractAlgebra.Generic.morphism_type(Qx, Qx) === typeof(f)
  @test f(x) == 1
  @test parent(f(x)) === Qx
  @test f(x + 1) == 2

  f = hom(Qx, Qx, c -> zero(Qx), x)
  @test is_zero(f(x))
  @test is_zero(f(x + 1))
  @test sprint(show, f) isa String
  @test sprint(show, "text/plain", f) isa String
  ff = hom(Qx, Qx, c -> c, x)
  @test (f * ff)(x) == 0
  @test (ff * f)(x) == 0
  ff = hom(Qx, Qx, c -> one(GF2), x)
  @test_throws ErrorException f * ff
  ff = hom(Qx, Qx, map_from_func(x -> Qx(x), QQ, Qx), x)
  @test (f * ff)(x + 1) == 0
  @test (ff * f)(x + 1) == 0
  ff = hom(Qx, Qx, map_from_func(x -> x, QQ, QQ), x)
  @test (f * ff)(x + 1) == 0
  @test (ff * f)(x + 1) == 0

  f = identity_map(Qx)
  ff = hom(Qx, Qx, map_from_func(x -> Qx(x), QQ, Qx), x)
  @test (f * ff)(x + 1) == x + 1
  @test (ff * f)(x + 1) == x + 1
  ff = hom(Qx, Qx, map_from_func(x -> x, QQ, QQ), x)
  @test (f * ff)(x + 1) == x + 1
  @test (ff * f)(x + 1) == x + 1

  Zy, y = ZZ["y"]
  @test_throws Exception hom(Qx, Qx, y)
  @test_throws Exception hom(Qx, Qx, identity, y)

  f = hom(Qx, Zy, y + 1)
  @test f(x^2) == y^2 + 2y + 1

  g = hom(Zy, Qx, x^2)
  @test (g * f)(y) == (y + 1)^2
  @test (f * g)(x) == x^2 + 1

  R = Qx
  Ry, y = R["y"]
  g = hom(Qx, Qx, x + 1)
  f = hom(Ry, Ry, g, y^2)
  @test f(x * y) == (x + 1) * y^2
  (f * f)(x * y) == (x + 2) * y^4
  @test sprint(show, f) isa String
  @test sprint(show, "text/plain", f) isa String

  ff = hom(Ry, Ry, identity, y^2)
  @test AbstractAlgebra.Generic.morphism_type(Ry, Ry, identity) === typeof(ff)
  @test ff(x * y) == x*y^2
  @test (f * ff)(x * y) == (x + 1) * y^4
  @test (ff * f)(x * y) == (x + 1) * y^4

  fff = hom(Ry, Ry, hom(Qx, Ry, y), y^3)
  @test fff(x * y) == y^4
  @test (f * fff)(x * y) == (y + 1)*y^6
  @test (fff * f)(x * y) == y^8

  Qx, x = QQ["x"]
  S = matrix_ring(ZZ, 2)
  m = S([1 1; 0 1])
  h = hom(Qx, S, S([1 1; 0 1]))
  @test h(x^2) == S([1 2; 0 1])
  f = hom(Qx, Qx, x + 1)
  @test (f * h)(x^2) == (m + 1)^2
end
