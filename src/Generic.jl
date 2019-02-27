module Generic

import LinearAlgebra: det, norm,
                      nullspace, rank, transpose!, hessenberg

if VERSION < v"1.0.0"
   import Base: atan2, base, contains, nextpow2, prevpow2
end

import LinearAlgebra: lu, lu!, tr

using Markdown, Random, InteractiveUtils

import Base: Array, abs, asin, asinh, atan, atanh, bin, checkbounds,
             conj, convert, cmp, cos, cosh, dec, deepcopy,
             deepcopy_internal, div, divrem,
             exp, exponent, gcd, gcdx, getindex, hash, hcat, hex, intersect, inv,
             invmod, isapprox, isequal, isfinite, isless, isqrt, isreal, iszero, lcm,
             ldexp, length, log, mod, ndigits,
             oct, one, parent, parse, precision,
             rand, Rational, rem, reverse, 
             setindex!, show, similar, sign, sin, sinh, size, string,
             tan, tanh, trailing_zeros, transpose, truncate,
             typed_hvcat, typed_hcat, vcat, xor, zero, zeros, +, -, *, ==, ^,
             &, |, <<, >>, ~, <=, >=, <, >, //, /, !=

if VERSION >= v"0.7.0-DEV.1144"
import Base: isone
end

import Base: floor, ceil, hypot, log, log1p, exp, expm1, sin, cos, sinpi,
             cospi, tan, cot, sinh, cosh, tanh, coth, atan, asin, acos, atanh,
             asinh, acosh, sinpi, cospi

import AbstractAlgebra: Integers, Rationals, NCRing, NCRingElem, Ring, RingElem,
       RingElement, Field, FieldElement, Map, promote_rule

using AbstractAlgebra

include("generic/GenericTypes.jl")

include("generic/PermGroups.jl")

include("generic/YoungTabs.jl")

include("generic/Residue.jl")

include("generic/ResidueField.jl")

include("generic/Poly.jl")

include("generic/NCPoly.jl")

include("generic/MPoly.jl")

include("generic/SparsePoly.jl")

include("generic/RelSeries.jl")

include("generic/AbsSeries.jl")

include("generic/LaurentSeries.jl")

include("generic/PuiseuxSeries.jl")

include("generic/Matrix.jl")

include("generic/MatrixAlgebra.jl")

include("generic/Fraction.jl")

include("generic/FreeModule.jl")

include("generic/Submodule.jl")

include("generic/QuotientModule.jl")

include("generic/ModuleHomomorphism.jl")

include("generic/Map.jl")

include("generic/MapWithInverse.jl")

include("generic/MapCache.jl")

end # generic
