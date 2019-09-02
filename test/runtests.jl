using AbstractAlgebra

using SparseArrays, LinearAlgebra
using AbstractAlgebra: mul! # disambiguate from LinearAlgebra.mul!

if VERSION < v"0.7.0-DEV.2004"
   using Base.Test
else
   using Test
end

include("../test/AbstractAlgebra-test.jl")

test_all()
