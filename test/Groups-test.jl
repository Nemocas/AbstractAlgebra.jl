include("generic/YoungTabs-test.jl")
include("generic/Perm-test.jl")
include("generic/PermGroupAPI-test.jl")

@testset "Convenience methods for comm and conj" begin
   # square and invertible
   m = matrix(ZZ, [1 -1 ; 0 1])
   @test conj(m, m) == m
   @test is_one(comm(m, m))

   # square but not invertible
   m2 = matrix(ZZ, [1 -1 ; 0 0])
   @test_throws ErrorException conj(m, m2)
   @test_throws ErrorException comm(m, m2)

   # not square
   m3 = matrix(ZZ, [1 2 3; 4 5 6])
   @test_throws DomainError conj(m, m3)
   @test_throws ErrorException comm(m, m3)
end
