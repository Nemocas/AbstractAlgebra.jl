@testset "Generic.OreOperator" begin
  R,x = polynomial_ring(QQ,:x)
  @testset "Differential operators" verbose=true begin
    δ = derivation(R)
    AD,D = ore_extension(R,:D,δ)

    @test order(D) == 1
    @test length(D) == 2

    a0 = zero(AD)
    @test order(a0) == -1
    @test length(a0) == 0
    @test a0.coeffs == elem_type(R)[]
    @test iszero(a0)
    @test !is_unit(a0)

    @test leading_coefficient(a0) |> iszero
    @test leading_term(a0)  |> iszero
    @test leading_monomial(a0) |> iszero

    a1 = one(AD)
    @test order(a1) == 0
    @test length(a1) == 1
    @test is_unit(a1)

    @test leading_coefficient(a1) |> isone
    @test leading_term(a1) |> isone
    @test leading_monomial(a1) |> isone

    @test a0 * a1 == a0
    @test a1 * a1 == a1

    a2 = D^2
    @test order(a2) == 2
    @test coefficients(a2) == R.([0,0,1])

    a3 = 3*a2 + D
    @test leading_coefficient(a3) == R(3)
    @test leading_monomial(a3) == D^2
    @test leading_term(a3) == leading_coefficient(a3)*leading_monomial(a3)

#    @testset "Division with remainder" verbose=true begin
#      @test rdivrem(D,D) == (one(AD),zero(AD))
#      @test rdivrem(D^2+D,D) == (D+1,zero(AD))
#    end

    @testset "Arithmetic" verbose=true begin
      c4 = [x, x^3, x^2, x^5]
      a4 = (c4[1]*D + c4[2])*(c4[3]*D + c4[4])
      @test coefficients(a4) == [c4[1]*δ(c4[4]) + c4[2]*c4[4], c4[1]*δ(c4[3]) + c4[2]*c4[3] + c4[1]*c4[4], c4[1]*c4[3]]
      @test length(a4) == 3
      @test order(a4) == 2
    end
  end

  @testset "Euler operators" verbose=true begin
    σ = identity_map(R)
    δ = derivation(R,σ,x)
    AE,θ = ore_extension(R,:θ,δ)

    a1 = θ*x
    @test order(a1) == 1
    @test coefficients(a1) == [x,x]

    a2 = θ*x^2
    @test order(a2) == 1
    @test coefficients(a2) == [2*x^2, x^2]
  end

  @testset "Forward differences" verbose=true begin
    σ = hom(R,R,x+1)
    δ = derivation(R,σ)
    AF,F = ore_extension(R,:F,δ)

    a1 = F*x
    @test order(a1) == 1
    @test coefficients(a1) == [one(R),x+1]

    a2 = F*x^2
    @test order(a2) == 1
    @test coefficients(a2) == [2*x + 1, x^2 + 2*x + 1]
  end

  @testset "Recurrence operators" verbose=true begin
    σ = hom(R,R,x+1)
    δ = trivial_derivation(R,σ)
    AS,S = ore_extension(R,:S,δ)

#    @testset "Division with remainder" verbose=true begin
#      @test rdivrem(S,S) == (one(AS),zero(AS))
#      @test rdivrem(S^2+S,S) == (S+1,zero(AS))
#    end

    @testset "Arithmetic" verbose=true begin
      c1 = [x, x^3, x^2, x^5]
      a1 = (c1[1]*S + c1[2])*(c1[3]*S + c1[4])
      @test coefficients(a1) == [c1[2]*c1[4],c1[1]*σ(c1[4]) + c1[2]*c1[3],c1[1]*σ(c1[3])]
      @test length(a1) == 3
      @test order(a1) == 2
    end
  end
end

#@testset "Univariate Ore algebra over $name" verbose=true for (name,(R,x)) in [
#    ("univariate polynomial fraction field", begin
#      R,x = polynomial_ring(QQ,:x)
#      K = fraction_field(R)
#      K,K(x)
#     end),
#    ("univariate polynomial function field", rational_function_field(QQ)),
#  ]
#
#  @testset "Differential operators" verbose=true begin
#    δ = derivation(R)
#    σ = sigma_endomorphism(δ)
#    AD,D = ore_extension(R,:D,δ)
#
#    U = AD([3x+1,2x-3,-(2x+1),3x+5])
#    V = AD([-(3x+5),x+2])
#    q,r = rdivrem(U,V)
#    @test q*V + r == U
#
#    U = (x-1)*D^5 + 5*D^4
#    V = (x-1)*D^3 + (6 - 3x)*D^2 + (27x^2 + 9x -42)*D + 117 - 54x
#    g = grcd(U,V)
#    @test g == D - (2x^2-x-7)//((x+2)*(x+1)*(x-1))
#
#    U = (2x+3)*x*D^2 + 2*(4x^2 + 3x - 3)*D + 2*(4x^2 - 3)
#    V = (x+1)*(x-1)*D^2 + 2*(2x^2 - x - 2)*D + 2*(2x^2 - 2x - 1)
#    f = lclm(U,V)
#    @test iszero(rrem(f,U))
#    @test iszero(rrem(f,V))
#
#    @testset "Hermite normal form" verbose=true begin
#      H1 = AD[D-1 D; 0 D^2]
#      @test hnf(H1) == H1
#
#      H2 = AD[D-(3x^2+4)//(x*(3x+2)) 6*(x-2)//(x*(3x+2)); -(3-x)//(3x+2) D-11//(3x+2)]
#      @test hnf(H2) == AD[1 (3x+2)//(x-3)*D-11/(x-3); 0 D^2+(6-x^2)//(x*(x-3))*D+3*(x-2)//(x*(x-3))]
#    end
#
#    @testset "Uncoupling" verbose=true begin
#      H3 = R[(3x^2+4)//(x*(3x+2)) -6*(x-2)//(x*(3x+2)); (3-x)//(3x+2) 11//(3x+2)]
#      P3 = R[1 0; (3x^2+4)//(x*(3x+2)) -6*(x-2)//(x*(3x+2))]
#      @test uncouple(H3,δ) == P3
#      @test gauge_transform(H3,δ) == (σ.(P3)*H3 + δ.(P3))*inv(P3)
#    end
#  end
#  @testset "Recurrence operators" verbose=true begin
#    σ = hom(R,R,x+1)
#    δ = trivial_derivation(R,σ)
#    AS,S = ore_extension(R,:S,δ)
#
#    U = AS([3x+1,2x-3,-(2x+1),3x+5])
#    V = AS([-(3x+5),x+2])
#    q,r = rdivrem(U,V)
#    @test q*V + r == U
#
#    a1 = (x*S + x) * (2*x*S + x)
#    @test a1 == x*2*(x+1)*S^2 + (x*(x+1) + x*2*x)*S + x*x
#
#    U = S^7 + 5S^6 + 9S^5 +5S^4 - 5S^3 - 9S^2 - 5S - 1
#    V = (x+5)*S^5 + 6S^4 - 6S^3 - 2*(x-1)*S^2 - (x+19)*S + 2*(x+6)
#    g = grcd(U,V)
#    @test g == S^2 + 10*(x+4)//((x+2)*(2x+7))*S - (2x+9)*(x+6)//((2x+7)*(x+2))
#
#    U = (x+2)*(2x-1)*S^2 - 8*(x^2+3x-1)*S + 4*(x+4)*(2x+1)
#    V = (x+4)*(x+3)*S^2 - 2*(x+5)*(2x+5)*S + 4*(x+4)*(x+5)
#    f = lclm(U,V)
#    @test iszero(rrem(f,U))
#    @test iszero(rrem(f,V))
#
#    @testset "Division with remainder" verbose=true begin
#      @test rdivrem(S,S) == (one(AS),zero(AS))
#      @test rdivrem(S^2+S,S) == (S+1,zero(AS))
#    end
#
#    @testset "Arithmetic" verbose=true begin
#      c1 = [x, x^3, x^2, x^5]
#      a1 = (c1[1]*S + c1[2])*(c1[3]*S + c1[4])
#      @test coefficients(a1) == [c1[2]*c1[4],c1[1]*σ(c1[4]) + c1[2]*c1[3],c1[1]*σ(c1[3])]
#      @test length(a1) == 3
#      @test order(a1) == 2
#    end
#  end
#end
#
