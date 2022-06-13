include("Rings-conformance-tests.jl")
# TODO: add more

# add methods for test_elem here

function test_elem(R::AbstractAlgebra.Integers{BigInt})
   n = big(2)^rand(1:100)
   return rand(ZZ, -n:n)
end

function test_elem(R::AbstractAlgebra.Rationals{BigInt})
   n = big(2)^rand(1:100)
   m = big(2)^rand(1:100)
   return rand(ZZ, -n:n)//rand(ZZ, 1:m)
end

function test_elem(R::AbstractAlgebra.GFField{Int64})
   return rand(R)
end
