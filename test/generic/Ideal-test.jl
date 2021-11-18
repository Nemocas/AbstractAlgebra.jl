function spoly(f::T, g::T) where T <: MPolyElem
   fc = leading_coefficient(f)
   gc = leading_coefficient(g)
   mf = exponent_vector(f, 1)
   mg = exponent_vector(g, 1)
   c = lcm(fc, gc)
   llcm = max.(mf, mg)
   infl = [1 for i in 1:nvars(parent(f))]
   shiftf = llcm .- mf
   shiftg = llcm .- mg
   s = divexact(c, fc)*inflate(f, shiftf, infl) - divexact(c, gc)*inflate(g, shiftg, infl)
end

function gpoly(f::T, g::T) where T <: MPolyElem
   fc = leading_coefficient(f)
   gc = leading_coefficient(g)
   mf = exponent_vector(f, 1)
   mg = exponent_vector(g, 1)
   _, s, t = gcdx(fc, gc)
   llcm = max.(mf, mg)
   infl = [1 for i in 1:nvars(parent(f))]
   shiftf = llcm .- mf
   shiftg = llcm .- mg
   g = s*inflate(f, shiftf, infl) + t*inflate(g, shiftg, infl)
end

function testit(R, V)
   I = Ideal(R, V)
   G = I.gens
   for v in V
      r = normal_form(v, G)
      if r != 0
         println("v = ", v)
         println("r = ", r)
         println("I = ", G)
         println("wrong basis")
         return false
      end
   end
   for i = 1:length(G)
      for j = i + 1:length(G)
         s = spoly(G[i], G[j])
         r = normal_form(s, G)
         if r != 0
            println("s = ", s)
            println("i = ", i, ", j = ", j)
            println("r = ", r)
            println("I = ", G)
            println("spolys not zero")
            return false
         end
      end
   end
   for i = 1:length(G)
      for j = i + 1:length(G)
         g = gpoly(G[i], G[j])
         r = normal_form(g, G)
         if r != 0
            println("V = ", V)
            println("g = ", g)
            println("i = ", i, ", j = ", j)
            println("r = ", r)
            println("I = ", G)
            println("gpolys not zero")
            return false
         end
      end
   end
   return true
end

@testset "Generic.Ideal.ideal_reduction(multivariate)" begin
   R, (x, y) = PolynomialRing(ZZ, ["x", "y"]; ordering=:degrevlex)

   for V in [
      [-10*x^2*y - 2, 4*x^2*y^2],
      [R(8), -6*x*y^2 - 5],
      [9*x*y^2 - x*y, 18*y^2 - 2*y],
      [10*x^2*y + 2, 10*x^2*y^2 - 7],
      [-5*y + 6, 4*x^2*y^2 - 8*y],
      [R(6), -R(9), 7*y^2 + 3*x],
      [-8*y^2, 2*x*y^2 - 5],
      [-R(10), 8*x*y + 7*y^2, -10*y],
      [-6*x*y^2 + 4, -6*x^2*y - x^2],
      [-7*x^2*y^2 + 2*x^2*y, -x^2*y^2 - x],
      [R(0), 6*y^2 - 7, -8*x*y - 2*x],
      [4*x^2*y + 3*x*y, -4*x^2*y^2 + 2],
      [-5*x*y^2 - 2*y, 7*x^2*y - 5],
      [4*x^2*y^2 - 5*x, 9*x*y + 2*y^2],
      [-2*x^3*y - 7*x, -2*x*y^2 - 5*y^3, -2*y^3 - 8*x + y],
      [4*x^2*y - 9*x*y, 4*x^2*y + 3*y^2],
      [3*x^2*y - 3*y^2, 9*x^2*y + 7*x*y],
      [3*x^2*y - 6*x*y^2, -10*x*y + y^2],
      [8*x - 10*y, 3*x^2*y],
      [3*x*y^2 - 5*x, 9*x^2*y^2 + 5],
      [-6*x*y + 3*y^2, -6*x^2*y^3 + 5*x*y^3 + 7*x^2, -9*x^3 + 10],
      [3*x^2*y^2 + 5*x*y, 3*x^2*y^2 + 10*x],
      [-9*x^2*y^2 + 10*x*y, 3*x^2*y^2 - 5*x^2, 6*y^2],
      [R(0), -4*x - 3*y - 9, -9*x^3 + 3*x*y, -2*x^3*y^3],
      [7*x^3*y^3 + 9*y^2 - 5*x, -5*x^3*y^2 + 4*y, 4*x^2*y - 4*y^3],
      [-3*x^2*y - 5*x*y, 2*x^2, 9*x^2*y^2 - 7*y],
      [-6*x^2*y^2 + 10*x, 9*x^2*y^3 - 7*x*y^2, -5*x^2 - 2*x*y],
      [6*x^3*y - 4*x*y^2, 3*x^2*y^3, -x^3*y + 5*y^2 - 9*y],
      [7*x^3*y^3 - 10*x^2*y - 3*x*y^2, 9*x*y^3 + 3*x^3 - 9*y^2, 2*x*y^3 + 4*x^2*y],
      [5*x^3*y^3 - 5*x^3 + 5*x^2, -9*x^3 - 6*x, -3*x^3*y^3 - 6*x^3*y - 10*y^3],
      [-2*y^3 + 10*x - 2, -8*x*y^3 + 2*x, 6*x^3*y - 5*y^2 - 3*x],
      [3*x^2*y^2 + 10*x*y, 9*y^2 - 7*x + 7, 7*x^2*y^3 + 2*y^3]
   ]
      @test testit(R, V)
   end

   # random examples
   for i = 1:100
      n = rand(0:3)
      V = elem_type(R)[]
      for j = 1:n
         push!(V, rand(R, 0:3, 0:3, -10:10))
      end
      @test testit(R, V)
   end

   R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"]; ordering=:degrevlex)

   for i = 1:100
      n = rand(0:2)
      V = elem_type(R)[]
      for j = 1:n
         push!(V, rand(R, 0:2, 0:2, -10:10))
      end
      @test testit(R, V)
   end
end

