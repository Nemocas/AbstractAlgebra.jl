function mix_ideal(I::Ideal{T}) where T <: RingElement
   G = gens(I)
   R = base_ring(I)
   if length(G) == 0
      return I
   end
   if length(G) == 1
      n = rand(0:5)
      H = T[rand(-10:10)*G[1] for i in 1:n]
      if rand(0:1) == 1
         return Ideal(R, vcat([-G[1]], H))
      else
         return Ideal(R, vcat([G[1]], H))
      end
   end
   n = rand(0:length(G))
   H = T[sum(rand(-10:10)*G[i] for i = 1:length(G)) for j = 1:n]
   G = vcat(G, H)
   for i = 1:length(G)
      for j = 1:length(G)
         if i != j
            G[i] += rand(-10:10)*G[j]
         end
      end
   end
   return Ideal(R, G)   
end

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
   if Ideal(R, G) != I
      println("I = ", G)
      println("I not reduced")
      return false
   end
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

@testset "Generic.Ideal.constructors" begin
   I = Ideal(ZZ, 3, 5)
   S = parent(I)

   @test typeof(IdealSet(ZZ)) == Generic.IdealSet{BigInt}

   @test typeof(S) == Generic.IdealSet{BigInt}

   @test base_ring(S) == ZZ

   @test parent(I) == S

   @test elem_type(S) == Generic.Ideal{BigInt}

   @test parent_type(I) == Generic.IdealSet{BigInt}

   J = Ideal(ZZ)

   @test parent(J) == S
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

@testset "Generic.Ideal.ideal_reduction(univariate)" begin
   R, x = PolynomialRing(ZZ, "x")

   for i = 1:300
      n = rand(0:5)
      V = elem_type(R)[rand(R, 0:10, -10:10) for i in 1:n]
      I = Ideal(R, V...)

      for v in V
         @test normal_form(v, I) == 0
      end

      G = gens(I)

      for i = 2:length(G)
         @test length(G[i]) > length(G[i - 1])
         @test divides(leading_coefficient(G[i - 1]), leading_coefficient(G[i]))[1]
      end

      @test Ideal(R, gens(I)) == I
   end
end

@testset "Generic.Ideal.ideal_reduction(integer)" begin
   for i = 1:300
      n = rand(0:10)
      V = elem_type(ZZ)[rand(ZZ, -10:10) for i in 1:n]
      I = Ideal(ZZ, V...)
      G = gens(I)

      @test length(G) == 1 || (length(V) == 0 && length(G) == 0) || (iszero(V) || length(G) == 0)

      if !isempty(G)
         for v in V
            @test divides(v, G[1])[1]
         end
      end

      @test Ideal(ZZ, gens(I)) == I
   end
end

@testset "Generic.Ideal.ideal_reduction(Fp[x])" begin
   Fp = GF(31)
   R, x = PolynomialRing(Fp, "x")

   for i = 1:300
      n = rand(0:10)
      V = elem_type(R)[rand(R, 0:5) for i in 1:n]
      I = Ideal(R, V...)
      G = gens(I)

      @test length(G) == 1 || (length(V) == 0 && length(G) == 0) || (iszero(V) || length(G) == 0)

      if !isempty(G)
         for v in V
            @test divides(v, G[1])[1]
         end
      end

      @test Ideal(R, gens(I)) == I
   end
end

@testset "Generic.Ideal.comparison" begin
   # multivariate
   R, (x, y) = PolynomialRing(ZZ, ["x", "y"]; ordering=:degrevlex)

   # random examples
   for i = 1:100
      n = rand(0:3)
      V = elem_type(R)[]
      for j = 1:n
         push!(V, rand(R, 0:3, 0:3, -10:10))
      end
      
      I = Ideal(R, V)

      @test I == mix_ideal(I)
   end

   # univariate
   R, x = PolynomialRing(ZZ, "x")

   for i = 1:300
      n = rand(0:5)
      V = elem_type(R)[rand(R, 0:10, -10:10) for i in 1:n]
      I = Ideal(R, V)

      @test I == mix_ideal(I)
   end

   # Fp[x]
   Fp = GF(31)
   R, x = PolynomialRing(Fp, "x")

   for i = 1:300
      n = rand(0:10)
      V = elem_type(R)[rand(R, 0:5) for i in 1:n]
      I = Ideal(R, V)

      @test I == mix_ideal(I)
   end

   # integer
   for i = 1:300
      n = rand(0:10)
      V = elem_type(ZZ)[rand(ZZ, -10:10) for i in 1:n]
      I = Ideal(ZZ, V)

      @test I == mix_ideal(I)
   end
end

@testset "Generic.Ideal.containment" begin
   # multivariate
   R, (x, y) = PolynomialRing(ZZ, ["x", "y"]; ordering=:degrevlex)

   # random examples
   for i = 1:100
      n = rand(0:3)
      m = rand(0:3)
      V = elem_type(R)[]
      W = elem_type(R)[]
      for j = 1:n
         push!(V, rand(R, 0:3, 0:3, -10:10))
      end
      for j = 1:m
         push!(W, rand(R, 0:3, 0:3, -10:10))
      end

      I = Ideal(R, V)
      J = Ideal(R, vcat(V, W))

      @test contains(J, I)
   end

   # univariate
   R, x = PolynomialRing(ZZ, "x")

   for i = 1:300
      n = rand(0:5)
      m = rand(0:5)
      V = elem_type(R)[rand(R, 0:10, -10:10) for i in 1:n]
      W = elem_type(R)[rand(R, 0:10, -10:10) for i in 1:m]

      I = Ideal(R, V)
      J = Ideal(R, vcat(V, W))

      @test contains(J, I)
   end

   # Fp[x]
   Fp = GF(31)
   R, x = PolynomialRing(Fp, "x")

   for i = 1:300
      n = rand(0:10)
      m = rand(0:10)
      V = elem_type(R)[rand(R, 0:5) for i in 1:n]
      W = elem_type(R)[rand(R, 0:5) for i in 1:m]

      I = Ideal(R, V)
      J = Ideal(R, vcat(V, W))

      @test contains(J, I)
   end

   # integer
   for i = 1:300
      n = rand(0:10)
      m = rand(0:10)
      V = elem_type(ZZ)[rand(ZZ, -10:10) for i in 1:n]
      W = elem_type(ZZ)[rand(ZZ, -10:10) for i in 1:m]

      I = Ideal(ZZ, V)
      J = Ideal(ZZ, vcat(V, W))

      @test contains(J, I)
   end

   I = Ideal(ZZ, 2)

   @test contains(I, Ideal(ZZ, BigInt[]))
   @test !contains(Ideal(ZZ, BigInt[]), I)
end

@testset "Generic.Ideal.addition" begin
   # multivariate
   R, (x, y) = PolynomialRing(ZZ, ["x", "y"]; ordering=:degrevlex)

   # random examples
   for i = 1:100
      n = rand(0:3)
      m = rand(0:3)
      V = elem_type(R)[]
      W = elem_type(R)[]
      for j = 1:n
         push!(V, rand(R, 0:3, 0:3, -10:10))
      end
      for j = 1:m
         push!(W, rand(R, 0:3, 0:3, -10:10))
      end

      I = Ideal(R, V)
      J = Ideal(R, W)

      @test contains(I + J, I)
      @test contains(I + J, J)
   end

   # univariate
   R, x = PolynomialRing(ZZ, "x")

   for i = 1:300
      n = rand(0:5)
      m = rand(0:5)
      V = elem_type(R)[rand(R, 0:10, -10:10) for i in 1:n]
      W = elem_type(R)[rand(R, 0:10, -10:10) for i in 1:m]

      I = Ideal(R, V)
      J = Ideal(R, W)

      @test contains(I + J, I)
      @test contains(I + J, J)
   end

   # Fp[x]
   Fp = GF(31)
   R, x = PolynomialRing(Fp, "x")

   for i = 1:300
      n = rand(0:10)
      m = rand(0:10)
      V = elem_type(R)[rand(R, 0:5) for i in 1:n]
      W = elem_type(R)[rand(R, 0:5) for i in 1:m]

      I = Ideal(R, V)
      J = Ideal(R, W)

      @test contains(I + J, I)
      @test contains(I + J, J)
   end

   # integer
   for i = 1:300
      n = rand(0:10)
      m = rand(0:10)
      V = elem_type(ZZ)[rand(ZZ, -10:10) for i in 1:n]
      W = elem_type(ZZ)[rand(ZZ, -10:10) for i in 1:m]

      I = Ideal(ZZ, V)
      J = Ideal(ZZ, W)

      @test contains(I + J, I)
      @test contains(I + J, J)
   end
end

@testset "Generic.Ideal.multiplication" begin
   # multivariate
   R, (x, y) = PolynomialRing(ZZ, ["x", "y"]; ordering=:degrevlex)

   # random examples
   for i = 1:50
      n = rand(0:3)
      m = rand(0:3)
      V = elem_type(R)[]
      W = elem_type(R)[]
      X = elem_type(R)[]
      for j = 1:n
         push!(V, rand(R, 0:2, 0:3, -10:10))
      end
      for j = 1:m
         push!(W, rand(R, 0:3, 0:2, -10:10))
      end
      for j = 1:m
         push!(X, rand(R, 0:3, 0:3, -10:10))
      end

      I = Ideal(R, V)
      J = Ideal(R, W)
      K = Ideal(R, X)

      @test I*(J + K) == I*J + I*K
   end

   # univariate
   R, x = PolynomialRing(ZZ, "x")

   for i = 1:300
      n = rand(0:5)
      m = rand(0:5)
      k = rand(0:5)
      V = elem_type(R)[rand(R, 0:5, -10:10) for i in 1:n]
      W = elem_type(R)[rand(R, 0:5, -10:10) for i in 1:m]
      X = elem_type(R)[rand(R, 0:5, -10:10) for i in 1:k]

      I = Ideal(R, V)
      J = Ideal(R, W)
      K = Ideal(R, X)

      @test I*(J + K) == I*J + I*K
   end

   # Fp[x]
   Fp = GF(31)
   R, x = PolynomialRing(Fp, "x")

   for i = 1:300
      n = rand(0:10)
      m = rand(0:10)
      k = rand(0:10)
      V = elem_type(R)[rand(R, 0:5) for i in 1:n]
      W = elem_type(R)[rand(R, 0:5) for i in 1:m]
      X = elem_type(R)[rand(R, 0:5) for i in 1:k]

      I = Ideal(R, V)
      J = Ideal(R, W)
      K = Ideal(R, X)

      @test I*(J + K) == I*J + I*K
   end

   # integer
   for i = 1:300
      n = rand(0:10)
      m = rand(0:10)
      k = rand(0:10)
      V = elem_type(ZZ)[rand(ZZ, -10:10) for i in 1:n]
      W = elem_type(ZZ)[rand(ZZ, -10:10) for i in 1:m]
      X = elem_type(ZZ)[rand(ZZ, -10:10) for i in 1:k]

      I = Ideal(ZZ, V)
      J = Ideal(ZZ, W)
      K = Ideal(ZZ, X)

      @test I*(J + K) == I*J + I*K
   end
end

@testset "Generic.Ideal.adhoc_multiplication" begin
   # multivariate
   R, (x, y) = PolynomialRing(ZZ, ["x", "y"]; ordering=:degrevlex)

   # random examples
   for i = 1:100
      n = rand(0:3)
      V = elem_type(R)[]
      for j = 1:n
         push!(V, rand(R, 0:3, 0:3, -10:10))
      end
      c = rand(R, 0:3, 0:3, -10:10)
      d = rand(R, 0:3, 0:3, -10:10)

      I = Ideal(R, V)

      @test I*c == Ideal(R, gens(I*c))
      @test c*I == Ideal(R, gens(I*c))
      @test c*I*d == d*I*c
   end

   # univariate
   R, x = PolynomialRing(ZZ, "x")

   for i = 1:300
      n = rand(0:5)
      V = elem_type(R)[rand(R, 0:10, -10:10) for i in 1:n]
      I = Ideal(R, V)
      c = rand(R, 0:10, -10:10)
      d = rand(R, 0:10, -10:10)

      @test I*c == Ideal(R, gens(I*c))
      @test c*I == Ideal(R, gens(I*c))
      @test c*I*d == d*I*c

      m = rand(ZZ, -10:10)

      @test m*I == I*m
   end

   # Fp[x]
   Fp = GF(31)
   R, x = PolynomialRing(Fp, "x")

   for i = 1:300
      n = rand(0:10)
      V = elem_type(R)[rand(R, 0:5) for i in 1:n]
      I = Ideal(R, V)
      c = rand(R, 0:5)
      d = rand(R, 0:5)

      @test I*c == Ideal(R, gens(I*c))
      @test c*I == Ideal(R, gens(I*c))
      @test c*I*d == d*I*c
   end

   # integer
   for i = 1:300
      n = rand(0:10)
      V = elem_type(ZZ)[rand(ZZ, -10:10) for i in 1:n]
      I = Ideal(ZZ, V)
      c = rand(ZZ, -10:10)
      d = rand(ZZ, -10:10)

      @test I*c == Ideal(ZZ, gens(I*c))
      @test c*I == Ideal(ZZ, gens(I*c))
      @test c*I*d == d*I*c
   end
end

@testset "Generic.Ideal.intersection" begin
   # multivariate
   R, (x, y) = PolynomialRing(ZZ, ["x", "y"]; ordering=:degrevlex)

   # random examples
   for i = 1:50
      n = rand(0:3)
      m = rand(0:3)
      V = elem_type(R)[]
      W = elem_type(R)[]
      for j = 1:n
         push!(V, rand(R, 0:2, 0:3, -10:10))
      end
      for j = 1:m
         push!(W, rand(R, 0:3, 0:2, -10:10))
      end

      I = Ideal(R, V)
      J = Ideal(R, W)
      
      K = intersection(I, J)

      @test contains(I, K)
      @test contains(J, K)
   end

   # univariate
   R, x = PolynomialRing(ZZ, "x")

   for i = 1:300
      n = rand(0:5)
      m = rand(0:5)
      V = elem_type(R)[rand(R, 0:5, -10:10) for i in 1:n]
      W = elem_type(R)[rand(R, 0:5, -10:10) for i in 1:m]

      I = Ideal(R, V)
      J = Ideal(R, W)

      K = intersection(I, J)

      @test contains(I, K)
      @test contains(J, K)
   end

   # Fp[x]
   Fp = GF(31)
   R, x = PolynomialRing(Fp, "x")

   for i = 1:300
      n = rand(0:10)
      m = rand(0:10)
      V = elem_type(R)[rand(R, 0:5) for i in 1:n]
      W = elem_type(R)[rand(R, 0:5) for i in 1:m]
      
      I = Ideal(R, V)
      J = Ideal(R, W)

      K = intersection(I, J)

      @test contains(I, K)
      @test contains(J, K)
   end

   # integer
   for i = 1:300
      n = rand(0:10)
      m = rand(0:10)
      V = elem_type(ZZ)[rand(ZZ, -10:10) for i in 1:n]
      W = elem_type(ZZ)[rand(ZZ, -10:10) for i in 1:m]
      
      I = Ideal(ZZ, V)
      J = Ideal(ZZ, W)
      
      K = intersection(I, J)

      @test contains(I, K)
      @test contains(J, K)
   end
end
