@testset "DensePoly.ZZ" begin
   R = ZZ
   Rx, x = R["x"]
   for blen in 0:16, clen in 0:16, alen in 0:(blen + clen)
      b = [rand(R, -10000:10000) for i in 1:blen]
      c = [rand(R, -10000:10000) for i in 1:clen]
      for dr in (-1, 0, 1)
         cutoff = rand(1:100)
         a = [rand(R, -10000:10000) for i in 1:alen]
         d = deepcopy(a)
         AbstractAlgebra.DensePoly.macc_classical(dr, a, 1, alen, b, 1, blen, c, 1, clen, R)
         AbstractAlgebra.DensePoly.macc(dr, d, 1, alen, b, 1, blen, c, 1, clen, R, cutoff)
         @test a == d
         if dr == 0
            @test Rx(a) == mullow(Rx(b), Rx(c), alen)
            AbstractAlgebra.DensePoly.mullow_fast!(d, alen, b, blen, c, clen, R, cutoff)
            @test a == d
         end
      end
   end
end

@testset "DensePoly.ZZ[y]" begin
   R, y = ZZ["y"]
   Rx, x = R["x"]
   for blen in 0:12, clen in 0:12, alen in 0:(blen + clen)
      b = [rand(R, 0:rand(1:9), -10000:10000) for i in 1:blen]
      c = [rand(R, 0:rand(1:9), -10000:10000) for i in 1:clen]
      for dr in (-1, 0, 1)
         cutoff = rand(1:100)
         a = [rand(R, 0:rand(1:9), -10000:10000) for i in 1:alen]
         d = deepcopy(a)
         AbstractAlgebra.DensePoly.macc_classical(dr, a, 1, alen, b, 1, blen, c, 1, clen, R)
         AbstractAlgebra.DensePoly.macc(dr, d, 1, alen, b, 1, blen, c, 1, clen, R, cutoff)
         @test a == d
         @test dr != 0 || Rx(a) == mullow(Rx(b), Rx(c), alen)
      end
   end
end

@testset "DensePoly.GF" begin
   R = GF(5)
   Rx, x = R["x"]
   for blen in 0:20, clen in 0:20, alen in 0:(blen + clen)
      b = [R(rand(-10000:10000)) for i in 1:blen]
      c = [R(rand(-10000:10000)) for i in 1:clen]
      for dr in (-1, 0, 1)
         cutoff = rand(1:100)
         a = [R(rand(-10000:10000)) for i in 1:alen]
         d = deepcopy(a)
         AbstractAlgebra.DensePoly.macc_classical(dr, a, 1, alen, b, 1, blen, c, 1, clen, R)
         AbstractAlgebra.DensePoly.macc(dr, d, 1, alen, b, 1, blen, c, 1, clen, R, cutoff)
         @test a == d
         @test dr != 0 || Rx(a) == mullow(Rx(b), Rx(c), alen)
      end
   end
end
