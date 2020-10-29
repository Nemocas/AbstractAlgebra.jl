using AbstractAlgebra

using SparseArrays, LinearAlgebra
using AbstractAlgebra: mul! # disambiguate from LinearAlgebra.mul!

using Random: Random, MersenneTwister, randsubseq, AbstractRNG
using RandomExtensions: make

const rng = MersenneTwister()

if VERSION < v"0.7.0-DEV.2004"
   using Base.Test
else
   using Test
end

include("AbstractAlgebra-test.jl")
