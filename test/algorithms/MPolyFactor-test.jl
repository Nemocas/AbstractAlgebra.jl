@testset "MPolyFactor" begin
   # not much to test because we would need working univariate factorization
   @testset "hlift_have_lcs" begin

      R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])

      fac = [y*x^2+z, (z+1)*x^3+x*y+z, (y*z+1)*x^2+1]

      ok, f = AbstractAlgebra.MPolyFactor.hlift_have_lcs(prod(fac),
                                          [x^2+1, 2*x^3+x+1, 2*x^2+1],
                                          [y, z+1, y*z+1],
                                          1, [2, 3], [QQ(1), QQ(1)])
      @test ok
      @test f == fac
   end

   @testset "hlift_with_lcc" begin

      R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])

      fac = [y*x^2+z, (z+1)*x^3+x*y+z, (y*z+1)*x^2+1]

      ok, f = AbstractAlgebra.MPolyFactor.hlift_with_lcc(prod(fac),
                                          [x^2+1, 2*x^3+x+1, 2*x^2+1],
                                          [y, z+1, y*z+1],
                                          1, [2, 3], [QQ(1), QQ(1)])
      @test ok
      @test f == fac
   end

   @testset "hlift_bivar_combine" begin

      R, (x, y) = polynomial_ring(QQ, ["x", "y"])

      p = y*(y*x+1)*((y+1)*x+y)*((y+2)*x+y)

      ok, content, fac = AbstractAlgebra.MPolyFactor.hlift_bivar_combine(p,
                                            1, 2, QQ(1), [x+1, 2*x+1, 3*x+1])
      @test ok
      @test degrees(content) == [0, 1]
      @test p == content*prod(fac)
   end

   @testset "make_bases_coprime!" begin

      R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
      P = Pair{elem_type(R), Int}

      # a, b: squarefree factorizations. Afterwards the bases across a and b
      # must be equal or coprime, those within each pairwise coprime, and the
      # products unchanged.
      function check(a, b)
         A = prod(p^e for (p, e) in a)
         B = prod(p^e for (p, e) in b)
         AbstractAlgebra.MPolyFactor.make_bases_coprime!(a, b)
         @test A == prod(p^e for (p, e) in a)
         @test B == prod(p^e for (p, e) in b)
         for (p, _) in a, (q, _) in b
            @test p == q || is_constant(gcd(p, q))
         end
         for c in (a, b), i in 1:length(c), j in i+1:length(c)
            @test is_constant(gcd(c[i].first, c[j].first))
         end
      end

      check(P[x*(x + y) => 1, y => 1], P[x => 2])
      check(P[y => 1, x*(x + y) => 1], P[x => 2])
      check(P[y => 1, x*(x + y) => 1], P[x => 2, y + 1 => 3])
      check(P[x*(x + y) => 1, y*(y + 1) => 2], P[(x + y)*(y + 1) => 3, x*y => 1])
   end

   R, (x, y) = ZZ[:x,:y]
   f = x^2
   fa = AbstractAlgebra.MPolyFactor.mfactor_squarefree_char_zero(x^2)
   @test is_unit(unit(fa))
   @test length(fa) == 1 && f == unit(fa) * prod(p^e for (p, e) in fa)
   fa = AbstractAlgebra.MPolyFactor.mfactor_squarefree_char_zero(x^2)
   @test is_unit(unit(fa))
   @test length(fa) == 1 && f == unit(fa) * prod(p^e for (p, e) in fa)
end
