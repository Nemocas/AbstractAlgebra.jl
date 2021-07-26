import AbstractAlgebra: NCRingElem, NCRingElement

const reps = 50

function equality(a::T, b::T) where T <: NCRingElement
   if isexact_type(T)
      return a == b
   else
      return isapprox(a, b)
   end
end

function test_NCRing_interface(
   R::AbstractAlgebra.NCRing,
   example_elements::Vector{T}
) where T <: NCRingElement

   @testset "NCRing interface for $(R)" begin

      @testset "Functions for types and parents of rings" begin
         @test elem_type(R) == T
         @test elem_type(typeof(R)) == T
         @test parent_type(T) == typeof(R)
         for x in example_elements
            @test parent(x) == R
         end
         @test isdomain_type(T) isa Bool
         @test isexact_type(T) isa Bool
         @test iszero(R(characteristic(R)))
         @test iszero(characteristic(R) * one(R))
         @test iszero(one(R) * characteristic(R))
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
         for x in example_elements
            @test R(x) isa T
         end
      end

      elems = deepcopy(example_elements)
      append!(elems, T[R(0), R(1), R(-2), R(BigInt(0)), R(BigInt(1)),
                       R(BigInt(-2))#=, R(BigInt(3)^100)=#])

      @testset "Basic functions" begin
         @test iszero(R())       # R() is supposed to construct 0 ?
         @test iszero(zero(R))
         @test isone(one(R))
         @test iszero(R(0))
         @test isone(R(1))
         @test isone(R(0)) || !isunit(R(0))
         @test isunit(R(1))
         for x in example_elements
            @test hash(x) isa UInt
            @test sprint(show, "text/plain", x) isa String
         end
         @test sprint(show, "text/plain", R) isa String

         for i in 1:reps
            a = rand(elems)
            b = rand(elems)
            c = rand(elems)
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

      if !(R isa AbstractAlgebra.Ring) && !iszero(one(R)) && isdomain_type(T)
         @testset "Basic functionality for noncommutative rings only" begin
            for i in 1:reps
               a = rand(elems)
               b = rand(elems)
               A = deepcopy(a)
               B = deepcopy(b)
               @test iszero(b) || equality(divexact_left(b*a, b), a)
               @test iszero(b) || equality(divexact_right(a*b, b), a)
               @test A == a
               @test B == b
            end
         end
      end

      @testset "Unsafe ring operators" begin
         for i in 1:reps
            a = rand(elems)
            b = rand(elems)
            c = rand(elems)
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

            # hmmm, doc says three args but describes four args 
            # BigFloat may fuse, so isapprox
            t = R()
            x = deepcopy(a)
            @test equality(addmul!(x, b, c, t), a+b*c)
            x = deepcopy(a)
            @test equality(addmul!(x, x, c, t), a+a*c)
            x = deepcopy(a)
            @test equality(addmul!(x, b, x, t), a+b*a)
            x = deepcopy(a)
            @test equality(addmul!(x, x, x, t), a+a*a)

            @test A == a
            @test B == b
            @test C == c
         end         
      end
   end
end


function test_Ring_interface(
   R::AbstractAlgebra.Ring,
   example_elements::Vector{T}
) where T <: NCRingElement

   test_NCRing_interface(R, example_elements)


   elems = deepcopy(example_elements)
   append!(elems, T[R(0), R(1), R(-2), R(BigInt(0)), R(BigInt(1)),
                    R(BigInt(-2))#=, R(BigInt(3)^100)=#])

   if !iszero(one(R)) && isdomain_type(T)
      @testset "Basic functionality for commutative rings only" begin
         for i in 1:reps
            a = rand(elems)
            b = rand(elems)
            A = deepcopy(a)
            B = deepcopy(b)
            @test isone(inv(one(R)))
            @test a*b == b*a
            @test iszero(b) || divexact(b*a, b) == a
            @test A == a
            @test B == b
         end
      end
   end
end


#InexactError on R(BigInt(3) ^ 100)
#failure on iszero(b) || divexact(b * a, b) == a
#test_Ring_interface(zz, Int[11, 2^30, 2^40])

test_Ring_interface(ZZ, BigInt[11])

test_Ring_interface(QQ, Rational{BigInt}[11//5])

#no method matching (::AbstractAlgebra.Floats{Float64})(::Bool)
#test_Ring_interface(RDF, Float64[11, 11//3, 1.2])

#no method matching (::AbstractAlgebra.Floats{BigFloat})(::Bool)
#test_Ring_interface(RealField, BigFloat[11, 11/3, 1.2, big(11), big(11)//big(3), big(1.2), Rational{BigInt}(11)])

R = GF(3)
test_Ring_interface(R, AbstractAlgebra.GFElem{Int}[R(2)])

#R = ResidueRing(ZZ, 1)
#falure on isone(one(R))
#test_Ring_interface(R, AbstractAlgebra.Generic.Res{BigInt}[R(2)])

R = GF(13)
test_Ring_interface(R, AbstractAlgebra.GFElem{Int}[R(11)])

R = GF(big(13))
test_Ring_interface(R, AbstractAlgebra.GFElem{BigInt}[R(11), R(big(11))])

R, x = ZZ["x"]
test_Ring_interface(R, Generic.Poly{BigInt}[x, x^2-big(2)^100, x^5])

