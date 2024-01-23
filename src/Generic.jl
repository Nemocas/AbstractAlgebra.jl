module Generic

import LinearAlgebra: det
import LinearAlgebra: hessenberg
import LinearAlgebra: issymmetric
import LinearAlgebra: istriu
import LinearAlgebra: lu
import LinearAlgebra: lu!
import LinearAlgebra: norm
import LinearAlgebra: nullspace
import LinearAlgebra: rank
import LinearAlgebra: tr

using Random, InteractiveUtils

using Random: SamplerTrivial, GLOBAL_RNG
using RandomExtensions: RandomExtensions, make, Make, Make2, Make3, Make4

import Base: abs
import Base: acos
import Base: acosh
import Base: Array
import Base: asin
import Base: asinh
import Base: atan
import Base: atanh
import Base: axes
import Base: bin
import Base: ceil
import Base: checkbounds
import Base: cmp
import Base: conj
import Base: convert
import Base: copy
import Base: cos
import Base: cosh
import Base: cospi
import Base: cot
import Base: coth
import Base: dec
import Base: deepcopy
import Base: deepcopy_internal
import Base: expm1
import Base: exponent
import Base: floor
import Base: gcd
import Base: gcdx
import Base: getindex
import Base: hash
import Base: hcat
import Base: hex
import Base: hypot
import Base: intersect
import Base: invmod
import Base: isapprox
import Base: isempty
import Base: isequal
import Base: isfinite
import Base: isless
import Base: isone
import Base: isqrt
import Base: isreal
import Base: iszero
import Base: lcm
import Base: ldexp
import Base: length
import Base: log1p
import Base: Matrix
import Base: mod
import Base: oct
import Base: one
import Base: parent
import Base: parse
import Base: powermod
import Base: precision
import Base: rand
import Base: Rational
import Base: rem
import Base: reverse
import Base: setindex!
import Base: show
import Base: sign
import Base: similar
import Base: sin
import Base: sinh
import Base: sinpi
import Base: size
import Base: string
import Base: tan
import Base: tanh
import Base: trailing_zeros
import Base: transpose
import Base: truncate
import Base: typed_hcat
import Base: typed_hvcat
import Base: vcat
import Base: xor
import Base: zero
import Base: zeros

import Base: +
import Base: -
import Base: *
import Base: ==
import Base: ^
import Base: &
import Base: |
import Base: <<
import Base: >>
import Base: ~
import Base: <=
import Base: >=
import Base: <
import Base: >
import Base: //
import Base: /
import Base: !=

# The type and helper function for the dictionaries for hashing
import ..AbstractAlgebra: @attributes
import ..AbstractAlgebra: @enable_all_show_via_expressify
import ..AbstractAlgebra: add!
import ..AbstractAlgebra: addeq!
import ..AbstractAlgebra: addmul!
import ..AbstractAlgebra: base_ring
import ..AbstractAlgebra: base_ring_type
import ..AbstractAlgebra: CacheDictType
import ..AbstractAlgebra: can_solve_with_solution_fflu
import ..AbstractAlgebra: can_solve_with_solution_lu
import ..AbstractAlgebra: canonical_unit
import ..AbstractAlgebra: change_base_ring
import ..AbstractAlgebra: characteristic
import ..AbstractAlgebra: check_parent
import ..AbstractAlgebra: codomain
import ..AbstractAlgebra: coeff
import ..AbstractAlgebra: coefficient_ring
import ..AbstractAlgebra: coefficients
import ..AbstractAlgebra: coefficients_of_univariate
import ..AbstractAlgebra: compose
import ..AbstractAlgebra: constant_coefficient
import ..AbstractAlgebra: content
import ..AbstractAlgebra: CycleDec
import ..AbstractAlgebra: data
import ..AbstractAlgebra: Dedent
import ..AbstractAlgebra: deflate
import ..AbstractAlgebra: deflation
import ..AbstractAlgebra: degree
import ..AbstractAlgebra: degrees
import ..AbstractAlgebra: degrees_range
import ..AbstractAlgebra: denominator
import ..AbstractAlgebra: derivative
import ..AbstractAlgebra: div
import ..AbstractAlgebra: divexact
import ..AbstractAlgebra: divides
import ..AbstractAlgebra: divrem
import ..AbstractAlgebra: domain
import ..AbstractAlgebra: elem_type
import ..AbstractAlgebra: evaluate
import ..AbstractAlgebra: exp
import ..AbstractAlgebra: exponent_vectors
import ..AbstractAlgebra: expressify
import ..AbstractAlgebra: factor
import ..AbstractAlgebra: factor_squarefree
import ..AbstractAlgebra: Field
import ..AbstractAlgebra: FieldElement
import ..AbstractAlgebra: gen
import ..AbstractAlgebra: gens
import ..AbstractAlgebra: get_cached!
import ..AbstractAlgebra: GFElem
import ..AbstractAlgebra: identity_matrix
import ..AbstractAlgebra: image_fn
import ..AbstractAlgebra: inflate
import ..AbstractAlgebra: Indent
import ..AbstractAlgebra: Integers
import ..AbstractAlgebra: integral
import ..AbstractAlgebra: inv
import ..AbstractAlgebra: is_constant
import ..AbstractAlgebra: is_domain_type
import ..AbstractAlgebra: is_exact_type
import ..AbstractAlgebra: is_finite
import ..AbstractAlgebra: is_gen
import ..AbstractAlgebra: is_monomial
import ..AbstractAlgebra: is_square
import ..AbstractAlgebra: is_square_with_sqrt
import ..AbstractAlgebra: is_term
import ..AbstractAlgebra: is_unit
import ..AbstractAlgebra: is_univariate
import ..AbstractAlgebra: is_zero_divisor
import ..AbstractAlgebra: is_zero_divisor_with_annihilator
import ..AbstractAlgebra: isreduced_form
import ..AbstractAlgebra: leading_coefficient
import ..AbstractAlgebra: leading_exponent_vector
import ..AbstractAlgebra: leading_monomial
import ..AbstractAlgebra: leading_term
import ..AbstractAlgebra: Lowercase
import ..AbstractAlgebra: log
import ..AbstractAlgebra: Map
import ..AbstractAlgebra: map_coefficients
import ..AbstractAlgebra: max_precision
import ..AbstractAlgebra: minpoly
import ..AbstractAlgebra: modulus
import ..AbstractAlgebra: monomials
import ..AbstractAlgebra: mul_classical
import ..AbstractAlgebra: mul_karatsuba
import ..AbstractAlgebra: mul!
import ..AbstractAlgebra: mullow
import ..AbstractAlgebra: number_of_columns
import ..AbstractAlgebra: number_of_digits
import ..AbstractAlgebra: number_of_generators
import ..AbstractAlgebra: NCRing
import ..AbstractAlgebra: NCRingElem
import ..AbstractAlgebra: number_of_rows
import ..AbstractAlgebra: numerator
import ..AbstractAlgebra: number_of_variables
import ..AbstractAlgebra: O
import ..AbstractAlgebra: order
import ..AbstractAlgebra: parent_type
import ..AbstractAlgebra: Perm
import ..AbstractAlgebra: pol_length
import ..AbstractAlgebra: preimage
import ..AbstractAlgebra: pretty
import ..AbstractAlgebra: primpart
import ..AbstractAlgebra: promote_rule
import ..AbstractAlgebra: pseudodivrem
import ..AbstractAlgebra: pseudorem
import ..AbstractAlgebra: Rationals
import ..AbstractAlgebra: reduced_form
import ..AbstractAlgebra: remove
import ..AbstractAlgebra: renormalize!
import ..AbstractAlgebra: Ring
import ..AbstractAlgebra: RingElem
import ..AbstractAlgebra: RingElement
import ..AbstractAlgebra: set_coefficient!
import ..AbstractAlgebra: set_length!
import ..AbstractAlgebra: set_precision!
import ..AbstractAlgebra: set_valuation!
import ..AbstractAlgebra: shift_left
import ..AbstractAlgebra: shift_right
import ..AbstractAlgebra: snf
import ..AbstractAlgebra: sqrt
import ..AbstractAlgebra: symbols
import ..AbstractAlgebra: tail
import ..AbstractAlgebra: term_degree
import ..AbstractAlgebra: terms
import ..AbstractAlgebra: terms_degrees
import ..AbstractAlgebra: to_univariate
import ..AbstractAlgebra: trailing_coefficient
import ..AbstractAlgebra: LowercaseOff
import ..AbstractAlgebra: use_karamul
import ..AbstractAlgebra: valuation
import ..AbstractAlgebra: var
import ..AbstractAlgebra: var_index
import ..AbstractAlgebra: vars
import ..AbstractAlgebra: zero!


using ..AbstractAlgebra

include("generic/GenericTypes.jl")

include("generic/PermGroups.jl")

include("generic/YoungTabs.jl")

include("generic/Residue.jl")

include("generic/ResidueField.jl")

include("generic/Poly.jl")

include("generic/NCPoly.jl")

include("generic/MPoly.jl")

include("generic/UnivPoly.jl")

include("generic/SparsePoly.jl")

include("generic/LaurentPoly.jl")

include("generic/LaurentMPoly.jl")

include("generic/RelSeries.jl")

include("generic/AbsSeries.jl")

include("generic/AbsMSeries.jl")

include("generic/LaurentSeries.jl")

include("generic/PuiseuxSeries.jl")

include("generic/Matrix.jl")

include("generic/MatrixAlgebra.jl")

include("generic/FreeAssAlgebra.jl")

include("generic/Fraction.jl")

include("generic/TotalFraction.jl")

include("generic/FactoredFraction.jl")

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

include("generic/Ideal.jl")

include("generic/AhoCorasick.jl")

include("generic/FreeAssAlgebraGroebner.jl")

###############################################################################
#
#   Temporary miscellaneous files being moved from Hecke.jl
#
###############################################################################

include("generic/Misc/Poly.jl")
include("generic/Misc/Rings.jl")
include("generic/Misc/Localization.jl")

# TODO/FIXME: deprecate aliases, remove in the future
import ..AbstractAlgebra: @alias
Base.@deprecate_binding ResF EuclideanRingResidueFieldElem
Base.@deprecate_binding ResField EuclideanRingResidueField
Base.@deprecate_binding Res EuclideanRingResidueRingElem
Base.@deprecate_binding ResRing EuclideanRingResidueRing

Base.@deprecate_binding Rat RationalFunctionFieldElem

Base.@deprecate_binding AbsSeriesRing AbsPowerSeriesRing
Base.@deprecate_binding AbsSeriesElem AbsPowerSeriesRingElem
Base.@deprecate_binding RelSeriesRing RelPowerSeriesRing
Base.@deprecate_binding RelSeriesElem RelPowerSeriesRingElem

# Deprecated in 0.34.*
@alias Frac FracFieldElem
@alias FactoredFrac FactoredFracFieldElem

# Deprecated in 0.35.*
Base.@deprecate_binding ResidueField EuclideanRingResidueField
Base.@deprecate_binding ResidueFieldElem EuclideanRingResidueFieldElem
Base.@deprecate_binding ResidueRing EuclideanRingResidueRing
Base.@deprecate_binding ResidueRingElem EuclideanRingResidueRingElem

end # generic
