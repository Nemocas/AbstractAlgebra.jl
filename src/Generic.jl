module Generic

import LinearAlgebra: det, issymmetric, norm,
                      nullspace, rank, transpose!, hessenberg

import LinearAlgebra: istriu, lu, lu!, tr

using Markdown, Random, InteractiveUtils

using Random: SamplerTrivial, GLOBAL_RNG
using RandomExtensions: RandomExtensions, make, Make, Make2, Make3, Make4

import Base: Array, abs, asin, asinh, atan, atanh, axes, bin, checkbounds, cmp, conj,
             convert, copy, cos, cosh, dec, deepcopy, deepcopy_internal,
             exponent, gcd, gcdx, getindex, hash, hcat, hex, intersect,
             invmod, isapprox, isempty, isequal, isfinite, isless, isone, isqrt,
             isreal, iszero, lcm, ldexp, length, Matrix, mod, ndigits, oct, one,
             parent, parse, powermod,
             precision, rand, Rational, rem, reverse, setindex!,
             show, similar, sign, sin, sinh, size, string, tan, tanh,
             trailing_zeros, transpose, truncate, typed_hvcat, typed_hcat,
             vcat, xor, zero, zeros, +, -, *, ==, ^, &, |, <<, >>, ~, <=, >=,
             <, >, //, /, !=

import Base: floor, ceil, hypot, log1p, expm1, sin, cos, sinpi, cospi,
             tan, cot, sinh, cosh, tanh, coth, atan, asin, acos, atanh, asinh,
             acosh, sinpi, cospi

# The type and helper function for the dictionaries for hashing
import ..AbstractAlgebra: CacheDictType, get_cached!

import ..AbstractAlgebra: Integers, Rationals, NCRing, NCRingElem, Ring,
                          RingElem, RingElement, Field, FieldElement,
                          isexact_type, isdomain_type, Map, promote_rule

import ..AbstractAlgebra: base_ring, canonical_unit, change_base_ring,
                          check_parent, denominator, div, divrem, elem_type,
                          exp, expressify, factor, gen, gens, inv, isgen,
                          ismonomial, ismonomial_recursive, isunit,
                          leading_coefficient, log, numerator, parent_type,
                          show_via_expressify, sqrt, trailing_coefficient

using ..AbstractAlgebra

include("generic/GenericTypes.jl")

include("generic/PermGroups.jl")

include("generic/YoungTabs.jl")

include("generic/Residue.jl")

include("generic/ResidueField.jl")

include("generic/Poly.jl")

include("generic/NCPoly.jl")

include("generic/MPoly.jl")

include("generic/SparsePoly.jl")

include("generic/LaurentPoly.jl")

include("generic/RelSeries.jl")

include("generic/AbsSeries.jl")

include("generic/AbsMSeries.jl")

include("generic/LaurentSeries.jl")

include("generic/PuiseuxSeries.jl")

include("generic/Matrix.jl")

include("generic/MatrixAlgebra.jl")

include("generic/Fraction.jl")

include("generic/RationalFunctionField.jl")

include("generic/FunctionField.jl")

include("generic/FreeModule.jl")

include("generic/Submodule.jl")

include("generic/QuotientModule.jl")

include("generic/DirectSum.jl")

include("generic/ModuleHomomorphism.jl")

include("generic/Module.jl")

include("generic/InvariantFactorDecomposition.jl")

include("generic/Map.jl")

include("generic/MapWithInverse.jl")

include("generic/MapCache.jl")

###############################################################################
#
#   Temporary miscellaneous files being moved from Hecke.jl
#
###############################################################################

include("generic/Misc/Poly.jl")
include("generic/Misc/Rings.jl")
include("generic/Misc/Localization.jl")

end # generic
