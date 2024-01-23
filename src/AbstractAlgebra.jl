@doc raw"""
 AbstractAlgebra is a pure Julia package for computational abstract algebra.

 For more information see https://github.com/Nemocas/AbstractAlgebra.jl
"""
module AbstractAlgebra

using Random: SamplerTrivial, GLOBAL_RNG
using RandomExtensions: RandomExtensions, make, Make, Make2, Make3, Make4

using InteractiveUtils

using Preferences

using Test # for "interface-conformance" functions

# A list of all symbols external packages should not import from AbstractAlgebra
const import_exclude = [:import_exclude, :QQ, :ZZ,
                  :RealField, :GF,
                  :AbstractAlgebra,
                  :inv, :log, :exp, :sqrt, :div, :divrem,
                  :numerator, :denominator,
                  :promote_rule,
                  :Set, :Module, :Group,
                 ]

# If you want to add methods to functions in LinearAlgebra they should be
# imported here and in Generic.jl, and exported below.
# They should not be imported/exported anywhere else.

import LinearAlgebra: det
import LinearAlgebra: hessenberg
import LinearAlgebra: ishermitian
import LinearAlgebra: issymmetric
import LinearAlgebra: isdiag
import LinearAlgebra: istril
import LinearAlgebra: istriu
import LinearAlgebra: lu
import LinearAlgebra: lu!
import LinearAlgebra: norm
import LinearAlgebra: nullspace
import LinearAlgebra: rank
import LinearAlgebra: tr

################################################################################
#
#  Import/export philosophy
#
#  For certain julia Base types and Base function, e.g. BigInt and div or exp, we
#  need a different behavior. These functions are not exported.
#
#  Take for example exp. Since exp is not imported from Base, there are two exp
#  functions, AbstractAlgebra.exp and Base.exp. Inside AbstractAlgebra, exp
#  will always refer to AbstractAlgebra.exp. When calling the function, one
#  should just use "exp" without namespace qualifcation.
#
#  On the other hand, if an AbstractAlgebra type wants to add a method to exp,
#  it must add a method to "Base.exp".
#
#  The rational for this is as follows: If we do "using AbstractAlgebra" in the
#  REPL, then "exp" will refer to the Base.exp. So if we want to make exp(a)
#  work in the REPL for an AbstractAlgebra type, we have to overload Base.exp.
#
################################################################################

# This is the list of functions for which we locally have a different behavior.
const Base_import_exclude = [:exp, :log, :sqrt, :inv, :div, :divrem, :numerator,
                             :denominator]

################################################################################
#
#  Functions that we do not import from Base
#
################################################################################

function exp(a::T) where T
   return Base.exp(a)
end

function log(a::T) where T
   return Base.log(a)
end

function sqrt(a::T; check::Bool=true) where T
  return Base.sqrt(a; check=check)
end

function divrem(a::T, b::T) where T
  return Base.divrem(a, b)
end

function div(a::T, b::T) where T
  return Base.div(a, b)
end

function inv(a::T) where T
  return Base.inv(a)
end

function numerator(a::T, canonicalise::Bool=true) where T
  return Base.numerator(a, canonicalise)
end

function denominator(a::T, canonicalise::Bool=true) where T
  return Base.denominator(a, canonicalise)
end

# If you want to add methods to functions in Base they should be imported here
# and in Generic.jl.
# They should not be imported/exported anywhere else.

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
import Base: fill
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
import Base: mod
import Base: ndigits
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
import Base: sincos
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
import Base: typed_vcat
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

using Random: Random, AbstractRNG, SamplerTrivial
using RandomExtensions: RandomExtensions, make, Make2

export AbsPowerSeriesRingElem
export add!
export addeq!
export AdditiveGroupElem
export crt
export crt_with_lcm
export elem_type
export ErrorConstrDimMismatch
export factor
export factor_squarefree
export Field
export FieldElem
export FieldElement
export FinField
export FinFieldElem
export FracElem
export FracField
export FreeAssAlgebra
export FreeAssAlgElem
export FunctionalMap
export Group
export GroupElem
export hgcd
export Ideal
export IdealSet
export IdentityMap
export InfiniteOrderError
export is_irreducible
export is_squarefree
export is_perfect
export ItemQuantity
export LaurentMPolyRing
export LaurentMPolyRingElem
export LaurentPolyRing
export LaurentPolyRingElem
export Map
export MatAlgebra
export MatAlgElem
export MatElem
export MatSpace
export ModuleElem
export MPolyRing
export MPolyRingElem
export mul!
export NCPolyRingElem
export NCRing
export NCRingElem
export NCRingElement
export NotImplementedError
export NotInvertibleError
export NumField
export NumFieldElem
export ordinal_number_string
export parent_type
export pluralize
export PolyRing
export PolyRingElem
export qq
export QQ
export RDF
export RealField
export RelPowerSeriesRingElem
export ResElem
export ResidueRing
export Ring
export RingElem
export RingElement
export SeriesElem
export SeriesRing
export SetElem
export SetMap
export SimpleNumField
export SimpleNumFieldElem
export sub!
export UniversalPolyRing
export UniversalPolyRingElem
export VarName
export zero!
export zeros
export zz
export ZZ

include("Attributes.jl")
include("AliasMacro.jl")
include("PrintHelper.jl")

# alternative names for some functions from Base
export is_empty
export is_equal
export is_even
export is_finite
export is_inf
export is_integer
export is_less
export is_odd
export is_one
export is_real
export is_subset
export is_valid
export is_zero
export number_of_digits

@alias is_empty isempty
@alias is_even iseven
@alias is_equal isequal
@alias is_finite isfinite
@alias is_inf isinf
@alias is_integer isinteger
@alias is_less isless
@alias is_odd isodd
@alias is_one isone
@alias is_real isreal
@alias is_subset issubset
@alias is_valid isvalid
@alias is_zero iszero
@alias number_of_digits ndigits

function order end

# alternative names for some functions from LinearAlgebra
# we don't use the `@alias` macro here because we provide custom
# docstrings for these aliases
const is_diagonal = isdiag
const is_hermitian = ishermitian
const is_symmetric = issymmetric
const is_lower_triangular = istril
const is_upper_triangular = istriu

# alternative names for some of our own functions
function number_of_columns end
function number_of_generators end
function number_of_rows end
function number_of_variables end
export number_of_columns
export number_of_generators
export number_of_rows
export number_of_variables
@alias ncols number_of_columns
@alias ngens number_of_generators
@alias nrows number_of_rows
@alias nvars number_of_variables

###############################################################################
# Macros for fancy printing. to use, enable attribute storage for your struct,
# i.e.m change
#
#   mutable struct bla..
#   ...
#   end
#
# to
#
#   @attributes mutable struct bla ..
#   ...
#   end
#
# Then, in the `show` method, start with
#   @show_name(io, obj)
# If the user assigned a name to the object (in the REPL mainly) by doing
# A = bla...
# then, in the compact printing only the name "A" is printed
# also adding
# @show_special(io, obj)
# allows, if present to call a different printing function for this instance
# See FreeModule for an example
#
###############################################################################

function set_name!(G::Any, name::String)
   set_attribute!(G, :name => name)
end

function set_name!(G)
   s = get_attribute(G, :name)
   s === nothing || return
   sy = find_name(G)
   sy === nothing && return
   set_name!(G, string(sy))
end

extra_name(G) = nothing

macro show_name(io, O)
  return :( begin
    local i = $(esc(io))
    local o = $(esc(O))
    s = get_attribute(o, :name)
    if s === nothing
      sy = find_name(o)
      if sy === nothing
        sy = extra_name(o)
      end
      if sy !== nothing
        s = string(sy)
        set_name!(o, s)
      end
    end
    if get(i, :compact, false) && s !== nothing
      if AbstractAlgebra.PrettyPrinting._supports_io_custom(i)
        print(i, LowercaseOff())
      end
      print(i, s)
      return
    end
  end )
end

const CurrentModule = Ref(Main)

function set_current_module(m)
  CurrentModule[] = m
end

function get_current_module()
  return CurrentModule[]
end

function find_name(A, M = Main; all::Bool = false)
  # in Documenter, the examples are not run in the REPL.
  # in the REPL: A = ... adds A to the global name space (Main....)
  # in Documenter (doctests) all examples are run in their own module
  # which is stored in CurrentModule, hence we need to search there as well
  #furthermore, they are not exported, hence the "all" option
  if M === Main && AbstractAlgebra.get_current_module() != Main
    a = find_name(A, AbstractAlgebra.get_current_module(), all = true)
    if a !== nothing
      return a
    end
  end
  for a = names(M, all = all)
    a === :ans && continue
    if isdefined(M, a) && getfield(M, a) === A
        return a
    end
  end
end

macro show_special(io, O)
  return :( begin
    local i = $(esc(io))
    local o = $(esc(O))
    s = get_attribute(o, :show)
    if s !== nothing
      s(i, o)
      return
    end
  end )
end

macro show_special_elem(io, e)
  return :( begin
    local i = $(esc(io))
    local a = $(esc(e))
    local o = parent(a)
    s = get_attribute(o, :show_elem)
    if s !== nothing
      s(i, a)
      return
    end
  end )
end

###############################################################################
# generic fall back if no immediate coercion is possible
# can/ should be called for more generic general coercion mechanisms

#tries to turn b into an element of a
# applications (in outside AbstractAlgebra so far)
#  - number fields (different cyclotomics, ie. coerce zeta_n into
#    cyclo(m*n)
#  - finite fields (although they roll their own)
#  - unram. local fields
#  - modules, abelian groups
#
# intended usage
# (a::Ring)(b::elem_type(a))
#   parent(b) == a && return a
#   return force_coerce(a, b)
#
function force_coerce(a, b, throw_error::Type{Val{T}} = Val{true}) where {T}
  if throw_error === Val{true}
    error("coercion not possible")
  end
  return nothing
end

#to allow +(a::T, b::T) where a, b have different parents, but
# a common over structure
# designed(?) to be minimally invasive in AA and Nemo, but filled with
# content in Hecke/Oscar
function force_op(op::Function, throw_error::Type{Val{T}}, a...) where {T}
  if throw_error === Val{true}
    error("no common overstructure for the arguments found")
  end
  return false
end

function force_op(op::Function, a...)
  return force_op(op, Val{true}, a...)
end

###############################################################################
#
#  Weak key id dictionaries
#
###############################################################################

include("WeakKeyIdDict.jl")

###############################################################################
#
#  Weak value dictionaries
#
###############################################################################

include("WeakValueDict.jl")

###############################################################################
#
#  Type for the Hash dictionary
#
###############################################################################

const CacheDictType = WeakValueDict

function get_cached!(default::Base.Callable, dict,
                                             key,
                                             use_cache::Bool)
   return use_cache ? Base.get!(default, dict, key) : default()
end

###############################################################################
#
#  Types
#
################################################################################

include("AbstractTypes.jl")

const PolynomialElem{T} = Union{PolyRingElem{T}, NCPolyRingElem{T}}
const MatrixElem{T} = Union{MatElem{T}, MatAlgElem{T}}

###############################################################################
#
#   Julia types
#
###############################################################################

include("julia/JuliaTypes.jl")

###############################################################################
#
#   Fundamental interface for AbstractAlgebra
#
###############################################################################

include("fundamental_interface.jl")

################################################################################
#
#   Printing
#
################################################################################

include("PrettyPrinting.jl")

import .PrettyPrinting: @enable_all_show_via_expressify
import .PrettyPrinting: allow_unicode
import .PrettyPrinting: canonicalize
import .PrettyPrinting: expr_to_latex_string
import .PrettyPrinting: expr_to_string
import .PrettyPrinting: expressify
import .PrettyPrinting: get_html_as_latex
import .PrettyPrinting: get_syntactic_sign_abs
import .PrettyPrinting: is_syntactic_one
import .PrettyPrinting: is_syntactic_zero
import .PrettyPrinting: is_unicode_allowed
import .PrettyPrinting: obj_to_latex_string
import .PrettyPrinting: obj_to_string
import .PrettyPrinting: obj_to_string_wrt_times
import .PrettyPrinting: print_integer_string
import .PrettyPrinting: print_obj
import .PrettyPrinting: printer
import .PrettyPrinting: set_html_as_latex
import .PrettyPrinting: show_obj
import .PrettyPrinting: show_via_expressify
import .PrettyPrinting: with_unicode
import .PrettyPrinting: pretty
import .PrettyPrinting: LowercaseOff
import .PrettyPrinting: Lowercase
import .PrettyPrinting: Indent
import .PrettyPrinting: Dedent

export @enable_all_show_via_expressify

###############################################################################
#
#   Generic algorithms defined on abstract types
#
###############################################################################

include("algorithms/LaurentPoly.jl")
include("algorithms/FinField.jl")
include("algorithms/GenericFunctions.jl")

include("CommonTypes.jl") # types needed by AbstractAlgebra and Generic
include("Poly.jl")
include("NCPoly.jl")
include("Matrix.jl")
include("Matrix-Strassen.jl")
include("MatrixAlgebra.jl")
include("AbsSeries.jl")
include("RelSeries.jl")
include("LaurentPoly.jl")
include("FreeModule.jl")
include("Submodule.jl")
include("QuotientModule.jl")
include("Module.jl")
include("InvariantFactorDecomposition.jl")
include("DirectSum.jl")
include("Map.jl")
include("MapCache.jl")
include("MapWithInverse.jl")
include("ModuleHomomorphism.jl")
include("Ideal.jl")
include("YoungTabs.jl")
include("PermGroups.jl")
include("LaurentSeries.jl")
include("PuiseuxSeries.jl")
include("SparsePoly.jl")
include("AbsMSeries.jl")
include("RationalFunctionField.jl")
include("Residue.jl")
include("ResidueField.jl")
include("Fraction.jl")
include("TotalFraction.jl")
include("MPoly.jl")
include("UnivPoly.jl")
include("FreeAssAlgebra.jl")
include("LaurentMPoly.jl")
include("MatrixNormalForms.jl")

###############################################################################
#
#   Generic submodule
#
###############################################################################

include("Generic.jl")

# Do not import div, divrem, exp, inv, log, sqrt, numerator and denominator
# as we have our own
import .Generic: @perm_str
import .Generic: abs_series_type
import .Generic: base_field
import .Generic: basis
import .Generic: character
import .Generic: collength
import .Generic: combine_like_terms!
import .Generic: cycles
import .Generic: defining_polynomial
import .Generic: degrees
import .Generic: dense_matrix_type
import .Generic: dim
import .Generic: disable_cache!
import .Generic: downscale
import .Generic: EuclideanRingResidueField
import .Generic: EuclideanRingResidueRing
import .Generic: enable_cache!
import .Generic: exp_gcd
import .Generic: exponent
import .Generic: exponent_vector
import .Generic: exponent_word
import .Generic: finish
import .Generic: fit!
import .Generic: function_field
import .Generic: gcd
import .Generic: gcdx
import .Generic: groebner_basis
import .Generic: has_bottom_neighbor
import .Generic: has_left_neighbor
import .Generic: hash
import .Generic: hooklength
import .Generic: image_fn
import .Generic: image_map
import .Generic: interreduce!
import .Generic: inv!
import .Generic: inverse_fn
import .Generic: inverse_image_fn
import .Generic: inverse_mat
import .Generic: invmod
import .Generic: is_compatible
import .Generic: is_divisible_by
import .Generic: is_homogeneous
import .Generic: is_power
import .Generic: is_rimhook
import .Generic: is_submodule
import .Generic: is_unit
import .Generic: isone
import .Generic: laurent_ring
import .Generic: laurent_series
import .Generic: lcm
import .Generic: leading_coefficient
import .Generic: leading_exponent_vector
import .Generic: leading_exponent_word
import .Generic: leading_monomial
import .Generic: leading_term
import .Generic: leglength
import .Generic: length
import .Generic: main_variable
import .Generic: main_variable_extract
import .Generic: main_variable_insert
import .Generic: map1
import .Generic: map2
import .Generic: matrix_repr
import .Generic: max_fields
import .Generic: mod
import .Generic: monomial
import .Generic: monomial_iszero
import .Generic: monomial_set!
import .Generic: monomial!
import .Generic: monomials
import .Generic: MPolyBuildCtx
import .Generic: mullow_karatsuba
import .Generic: norm
import .Generic: normal_form
import .Generic: normalise
import .Generic: num_coeff
import .Generic: one
import .Generic: order
import .Generic: ordering
import .Generic: parity
import .Generic: partitionseq
import .Generic: perm
import .Generic: permtype
import .Generic: polcoeff
import .Generic: poly
import .Generic: poly_ring
import .Generic: precision
import .Generic: preimage_map
import .Generic: prime
import .Generic: push_term!
import .Generic: reduce!
import .Generic: rel_series_type
import .Generic: rels
import .Generic: rescale!
import .Generic: retraction_map
import .Generic: reverse
import .Generic: rising_factorial
import .Generic: rising_factorial2
import .Generic: rowlength
import .Generic: section_map
import .Generic: set_exponent_vector!
import .Generic: set_exponent_word!
import .Generic: set_limit!
import .Generic: setcoeff!
import .Generic: setpermstyle
import .Generic: size
import .Generic: sort_terms!
import .Generic: summands
import .Generic: supermodule
import .Generic: term
import .Generic: terms
import .Generic: to_univariate
import .Generic: total_degree
import .Generic: trailing_coefficient
import .Generic: truncate
import .Generic: unit
import .Generic: upscale
import .Generic: weights
import .Generic: zero

# Moved from Hecke into Misc
import .Generic: Loc
import .Generic: localization
import .Generic: LocElem
import .Generic: roots
import .Generic: sturm_sequence

###############################################################################
#
#   Linear solving submodule
#
###############################################################################

include("Solve.jl")

# Do not export inv, div, divrem, exp, log, sqrt, numerator and denominator as we define our own
export _check_dim
export _checkbounds
export @alias
export @attr
export @attributes
export @free_associative_algebra
export @laurent_polynomial_ring
export @perm_str
export @polynomial_ring
export @power_series_ring
export @rational_function_field
export abs_series
export abs_series_type
export AbsPowerSeriesRing
export add_column
export add_column!
export add_row
export add_row!
export addmul_delayed_reduction!
export addmul!
export AllParts
export AllPerms
export allow_unicode
export base_field
export base_ring
export base_ring_type
export basis
export block_diagonal_matrix
export cached
export can_solve
export can_solve_left_reduced_triu
export can_solve_with_kernel
export can_solve_with_solution
export can_solve_with_solution_interpolation
export canonical_unit
export change_base_ring
export change_coefficient_ring
export character
export characteristic
export charpoly
export charpoly_danilevsky_ff!
export charpoly_danilevsky!
export charpoly_hessenberg!
export chebyshev_t
export chebyshev_u
export check_composable
export check_parent
export codomain
export coeff
export coefficient_ring
export coefficients
export coefficients_of_univariate
export collength
export combine_like_terms!
export comm
export compose
export constant_coefficient
export content
export cycles
export data
export defining_polynomial
export deflate
export deflation
export degree
export degrees
export denest
export dense_matrix_type
export dense_poly_ring_type
export dense_poly_type
export derivative
export det
export det_popov
export diagonal_matrix
export dim
export direct_sum
export disable_cache!
export discriminant
export divexact
export divexact_left
export divexact_low
export divexact_right
export divhigh
export divides
export domain
export downscale
export echelon_form
export echelon_form_with_transformation
export elem_type
export enable_cache!
export EuclideanRingResidueField
export EuclideanRingResidueRing
export evaluate
export exp_gcd
export exponent
export exponent_vector
export exponent_vectors
export exponent_word
export exponent_words
export extended_weak_popov
export extended_weak_popov_with_transform
export exterior_power
export Fac
export FactoredFractionField
export fflu
export fflu!
export find_pivot_popov
export finish
export fit!
export fraction_field
export free_associative_algebra
export free_module
export FreeModule
export function_field
export gcd
export gcd_with_cofactors
export gcdinv
export gcdx
export gen
export gens
export get_attribute
export get_attribute!
export gram
export has_attribute
export has_bottom_neighbor
export has_gens
export has_left_neighbor
export hash
export hermite_form
export hermite_form_with_transformation
export hessenberg
export hessenberg!
export hnf
export hnf_cohen
export hnf_cohen_with_transform
export hnf_kb
export hnf_kb_with_transform
export hnf_kb!
export hnf_minors
export hnf_minors_with_transform
export hnf_via_popov
export hnf_via_popov_with_transform
export hnf_with_transform
export hooklength
export ideal
export identity_map
export identity_matrix
export image
export image_fn
export image_map
export inflate
export integral
export interpolate
export inv!
export invariant_factors
export inverse_fn
export inverse_image_fn
export inverse_mat
export invmod
export is_compatible
export is_constant
export is_degree
export is_diagonal
export is_divisible_by
export is_domain_type
export is_exact_type
export is_finiteorder
export is_gen
export is_hermitian
export is_hessenberg
export is_hnf
export is_homogeneous
export is_invertible
export is_invertible_with_inverse
export is_isomorphic
export is_lower_triangular
export is_monic
export is_monomial
export is_monomial_recursive
export is_negative
export is_popov
export is_positive
export is_reverse
export is_rimhook
export is_rref
export is_skew_symmetric
export is_snf
export is_square
export is_submodule
export is_symmetric
export is_term
export is_term_recursive
export is_trivial
export is_unicode_allowed
export is_unit
export is_univariate
export is_upper_triangular
export is_weak_popov
export is_zero_column
export is_zero_divisor
export is_zero_divisor_with_annihilator
export is_zero_entry
export is_zero_row
export kernel
export kronecker_product
export laurent_ring
export laurent_series
export laurent_series_field
export laurent_series_ring
export laurent_polynomial_ring
export lcm
export leading_coefficient
export leading_exponent_vector
export leading_exponent_word
export leading_monomial
export leading_term
export left_kernel
export leglength
export length
export lift
export lower_triangular_matrix
export lu
export lu!
export main_variable
export main_variable_extract
export main_variable_insert
export map_coefficients
export map_entries
export map_entries!
export map_from_func
export map_with_preimage_from_func
export map_with_retraction
export map_with_retraction_from_func
export map_with_section
export map_with_section_from_func
export map1
export map2
export mat
export matrix
export matrix_repr
export matrix_space
export MatrixAlgebra
export MatrixElem
export max_fields
export max_precision
export minors
export minpoly
export mod
export module_homomorphism
export module_isomorphism
export ModuleHomomorphism
export ModuleIsomorphism
export modulus
export monomial
export monomial_iszero
export monomial_set!
export monomial_to_newton!
export monomial!
export monomials
export mpoly_type
export mpoly_ring_type
export MPolyBuildCtx
export mul_classical
export mul_karatsuba
export mul_ks
export mul_red!
export mulhigh_n
export mullow
export mullow_karatsuba
export mulmod
export multiply_column
export multiply_column!
export multiply_row
export multiply_row!
export newton_to_monomial!
export norm
export normal_form
export normalise
export nullspace
export num_coeff
export O
export one
export order
export ordering
export parent_type
export parity
export Partition
export partitionseq
export perm
export Perm
export permtype
export pfaffian
export pfaffians
export pol_length
export polcoeff
export poly
export poly_ring
export PolyCoeffs
export polynomial
export polynomial_ring
export polynomial_to_power_sums
export PolynomialElem
export PolyRing
export popov
export popov_with_transform
export pow_multinomial
export power_series_ring
export power_sums_to_polynomial
export powers
export ppio
export precision
export preimage
export preimage_map
export prime
export primpart
export pseudo_inv
export pseudodivrem
export pseudorem
export puiseux_series_field
export puiseux_series_ring
export push_term!
export quo
export rand_ordering
export randmat_triu
export randmat_with_rank
export rank
export rank_profile_popov
export rational_function_field
export reduce!
export rel_series
export rel_series_type
export RelPowerSeriesRing
export rels
export remove
export renest
export renormalize!
export rescale!
export residue_field
export residue_ring
export resultant
export resultant_ducos
export resultant_euclidean
export resultant_lehmer
export resultant_subresultant
export resultant_sylvester
export resx
export retraction_map
export reverse
export reverse_cols
export reverse_cols!
export reverse_rows
export reverse_rows!
export right_kernel
export rising_factorial
export rising_factorial2
export rowlength
export rref
export rref_rational
export rref_rational!
export rref!
export scalar_matrix
export section_map
export set_attribute!
export set_coefficient!
export set_exponent_vector!
export set_exponent_word!
export set_length!
export set_limit!
export set_precision!
export set_valuation!
export setcoeff!
export setpermstyle
export shift_left
export shift_right
export similarity!
export size
export SkewDiagram
export snf
export snf_kb
export snf_kb_with_transform
export snf_kb!
export snf_with_transform
export solve
export solve_ff
export solve_left
export solve_rational
export solve_triu
export solve_with_det
export sort_terms!
export SparsePolynomialRing
export Strassen
export strictly_lower_triangular_matrix
export strictly_upper_triangular_matrix
export sub
export subst
export summands
export supermodule
export swap_cols
export swap_cols!
export swap_rows
export swap_rows!
export sylvester_matrix
export symbols
export SymmetricGroup
export tail
export term
export terms
export to_univariate
export total_degree
export total_ring_of_fractions
export tr
export trailing_coefficient
export truncate
export typed_hcat
export typed_hvcat
export unit
export UniversalPolynomialRing
export upper_triangular_matrix
export upscale
export use_karamul
export valuation
export var
export var_index
export vars
export vector_space
export VectorSpace
export with_unicode
export weak_popov
export weak_popov_with_transform
export weights
export YoungTableau
export zero
export zero_matrix

# Moved from Hecke into Misc
export divexact_low
export divhigh
export is_monic
export Loc
export localization
export LocElem
export mulhigh_n
export PolyCoeffs
export roots
export sturm_sequence

################################################################################
#
#   Parent constructors
#
################################################################################

function SkewDiagram(lambda::Generic.Partition, mu::Generic.Partition)
  Generic.SkewDiagram(lambda, mu)
end

function YoungTableau(part::Generic.Partition, fill::Vector{T}=collect(1:part.n)) where T <: Integer
   Generic.YoungTableau(part, fill)
end

@doc raw"""
    sub(m::Module{T}, subs::Vector{<:Generic.Submodule{T}}) where T <: RingElement

Return the submodule `S` of the module `m` generated by the union of the given
submodules of $m$, and a map which is the canonical injection from `S` to `m`.
"""
function sub(m::Module{T}, subs::Vector{<:Generic.Submodule{T}}) where T <: RingElement
   Generic.sub(m, subs)
end

# Handles empty vector of submodules
function sub(m::Module{T}, subs::Vector{<:Generic.Submodule{U}}) where {T <: RingElement, U <: Any}
   Generic.sub(m, subs)
end

export Generic

###############################################################################
#
#   misc
#
###############################################################################

include("misc/ProductIterator.jl")
include("misc/Evaluate.jl")
include("misc/VarNames.jl")

###############################################################################
#
#   Polynomial Ring S, x = R[:x] syntax
#
###############################################################################

getindex(R::NCRing, s::VarName) = polynomial_ring(R, s)
# `R[:x, :y]` returns `S, [x, y]` instead of `S, x, y`
getindex(R::NCRing, s::VarName, ss::VarName...) =
   polynomial_ring(R, [Symbol(x) for x in (s, ss...)])

# syntax: Rxy, y = R[:x][:y]
getindex(R::Union{Tuple{PolyRing, PolyRingElem}, Tuple{NCPolyRing, NCPolyRingElem}}, s::VarName) = polynomial_ring(R[1], s)

###############################################################################
#
#   Syntax S[i] for all parents S as a shortcut for gen(S, i)
#
################################################################################

getindex(S::Set, i::Int) = gen(S, i)

###############################################################################
#
#   Matrix M = R[...] syntax
#
################################################################################

VERSION >= v"1.7" && (Base.typed_hvncat(R::NCRing, args...) = _matrix(R, hvncat(args...)))
Base.typed_hvcat(R::NCRing, args...) = _matrix(R, hvcat(args...))
Base.typed_hcat(R::NCRing, args...) = _matrix(R, hcat(args...))
Base.typed_vcat(R::NCRing, args...) = _matrix(R, vcat(args...))
_matrix(R::NCRing, a::AbstractVector) = matrix(R, length(a), isempty(a) ? 0 : 1, a)
_matrix(R::NCRing, a::AbstractMatrix) = matrix(R, a)

###############################################################################
#
#   Load error objects
#
###############################################################################

include("error.jl")

###############################################################################
#
#   Load Groups/Rings/Fields etc.
#
###############################################################################

include("Groups.jl")
include("Rings.jl")

# Generic and specific rings and fields
include("julia/Integer.jl")
include("julia/Rational.jl")
include("julia/Float.jl")
include("julia/GF.jl")

include("Fields.jl")
include("Factor.jl")

# Generic functions to be defined after all rings
include("polysubst.jl")

include("NCRings.jl")

include("broadcasting.jl")

################################################################################
#
#   Further functionality for Julia matrices
#
################################################################################

include("julia/Matrix.jl")

################################################################################
#
#   Deprecations
#
################################################################################

include("Deprecations.jl")

################################################################################
#
#   Stuff moved from Nemo
#
################################################################################

include("NemoStuff.jl")

###############################################################################
#
#   Array creation functions
#
###############################################################################

Array(R::NCRing, r::Int...) = Array{elem_type(R)}(undef, r)

function zeros(R::NCRing, r::Int...)
   T = elem_type(R)
   A = Array{T}(undef, r)
   for i in eachindex(A)
      A[i] = R()
   end
   return A
end

###############################################################################
#
#   Set domain for ZZ, QQ
#
###############################################################################

const ZZ = JuliaZZ
const QQ = JuliaQQ

###############################################################################
#
#   Set domain for RealField
#
###############################################################################

const RealField = JuliaRealField

###############################################################################
#
#   Generic algorithms defined on generic types
#
###############################################################################

include("algorithms/MPolyEvaluate.jl")
include("algorithms/MPolyFactor.jl")
include("algorithms/MPolyNested.jl")
include("algorithms/DensePoly.jl")

###############################################################################
#
#  For backwards compatibility
#
###############################################################################

include("Aliases.jl")

###############################################################################
#
#   Test code
#
###############################################################################

function test_module(x, y)
   julia_exe = Base.julia_cmd()
   pkgdir = realpath(joinpath(dirname(@__FILE__), ".."))
   test_file = joinpath(pkgdir, "test/$x/")
   test_file = test_file * "$y-test.jl";
   test_function_name = "test_"

   x == "generic"
   if y == "RelSeries"
      test_function_name *= "gen_rel_series"
   elseif y == "AbsSeries"
      test_function_name *= "gen_abs_series"
   elseif y == "Matrix"
      test_function_name *= "gen_mat"
   elseif y == "Fraction"
      test_function_name *= "gen_frac"
   elseif y == "Residue"
      test_function_name *= "gen_res"
   else
      test_function_name *= "gen_$(lowercase(y))"
   end

   cmd = """
         using Test
         using AbstractAlgebra
         include("test/rand.jl")
         include("$test_file")
         $test_function_name()
         """
   @info("spawning ", `$julia_exe --project=$(Base.active_project()) -e $cmd`)
   run(`$julia_exe -e $cmd`)
end

end # module
