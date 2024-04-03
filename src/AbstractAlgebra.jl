@doc raw"""
 AbstractAlgebra is a pure Julia package for computational abstract algebra.

 For more information see https://github.com/Nemocas/AbstractAlgebra.jl
"""
module AbstractAlgebra

using Random: Random, AbstractRNG, GLOBAL_RNG, SamplerTrivial
using RandomExtensions: RandomExtensions, make, Make, Make2, Make3, Make4

using InteractiveUtils: InteractiveUtils, subtypes

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

import LinearAlgebra

import LinearAlgebra: det
import LinearAlgebra: dot
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
# Do not import div, divrem, exp, inv, log, sqrt, numerator and denominator
# as we have our own
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
import Base: conj!
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
import Base: delete!
import Base: empty
import Base: expm1
import Base: exponent
import Base: fill
import Base: floor
import Base: gcd
import Base: gcdx
import Base: get
import Base: getindex
import Base: getkey
import Base: hash
import Base: hcat
import Base: hex
import Base: hypot
import Base: intersect
import Base: invmod
import Base: isempty
import Base: isequal
import Base: isfinite
import Base: isless
import Base: isone
import Base: isqrt
import Base: isreal
import Base: iszero
import Base: iterate
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
import Base: pop!
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
import Base: sizehint!
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

include("exports.jl")

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
const MatrixElem{T} = Union{MatElem{T}, MatRingElem{T}}

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

using .PrettyPrinting

import .PrettyPrinting: expressify

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
include("MatRing.jl")
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
import .Generic: EuclideanRingResidueFieldElem
import .Generic: EuclideanRingResidueRing
import .Generic: EuclideanRingResidueRingElem
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
import .Generic: internal_ordering
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
import .Generic: universal_polynomial_ring
import .Generic: upscale
import .Generic: weights
import .Generic: zero

# Moved from Hecke into Misc
import .Generic: LocalizedEuclideanRing
import .Generic: localization
import .Generic: LocalizedEuclideanRingElem
import .Generic: roots
import .Generic: sturm_sequence

###############################################################################
#
#   Linear solving submodule
#
###############################################################################

include("Solve.jl")

using .Solve

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

################################################################################
#
#   generic stubs for sub, quo, ...
#
################################################################################

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

function canonical_injections(D)
  return [canonical_injection(D, i) for i=1:_number_of_direct_product_factors(D)]
end

function canonical_projections(D)
  return [canonical_projections(D, i) for i=1:_number_of_direct_product_factors(D)]
end

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
