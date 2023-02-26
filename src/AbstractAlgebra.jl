@doc Markdown.doc"""
 AbstractAlgebra is a pure Julia package for computational abstract algebra.

 For more information see https://github.com/Nemocas/AbstractAlgebra.jl
"""
module AbstractAlgebra

using Random: SamplerTrivial, GLOBAL_RNG
using RandomExtensions: RandomExtensions, make, Make, Make2, Make3, Make4

using Markdown

using InteractiveUtils

using Test # for "interface-conformance" functions

import GroupsCore
import GroupsCore: gens, ngens, order, mul!, istrivial

# A list of all symbols external packages should not import from AbstractAlgebra
import_exclude = [:import_exclude, :QQ, :ZZ,
                  :RealField, :number_field,
                  :AbstractAlgebra,
                  :inv, :log, :exp, :sqrt, :div, :divrem,
                  :numerator, :denominator,
                  :promote_rule,
                  :Set, :Module, :Ring, :Group, :Field]

# If you want to add methods to functions in LinearAlgebra they should be
# imported here and in Generic.jl, and exported below.
# They should not be imported/exported anywhere else.

import LinearAlgebra: det
import LinearAlgebra: hessenberg
import LinearAlgebra: ishermitian
import LinearAlgebra: issymmetric
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
export error_dim_negative
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
export GroupElem
export hgcd
export Ideal
export IdealSet
export IdentityMap
export is_irreducible
export is_squarefree
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
export NCRingElem
export NotImplementedError
export NotInvertibleError
export NumField
export NumFieldElem
export parent_type
export PolyRing
export PolyRingElem
export qq
export QQ
export RDF
export RealField
export RelPowerSeriesRingElem
export ResElem
export ResidueRing
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
export zero!
export zeros
export zz
export ZZ

include("Attributes.jl")
include("AliasMacro.jl")

# alternative names for some functions from Base
export is_empty, is_equal, is_finite, is_inf, is_integer, is_less, is_one, is_real, is_subset, is_valid, is_zero

@alias is_empty isempty
@alias is_equal isequal
@alias is_finite isfinite
@alias is_inf isinf
@alias is_integer isinteger
@alias is_less isless
@alias is_one isone
@alias is_real isreal
@alias is_subset issubset
@alias is_valid isvalid
@alias is_zero iszero

# alternative names for some functions from GroupsCore
export is_trivial

@alias is_trivial istrivial

# alternative names for some functions from LinearAlgebra
export is_hermitian, is_symmetric, is_upper_triangular

@alias is_hermitian ishermitian
@alias is_symmetric issymmetric
@alias is_upper_triangular istriu

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
    if get(i, :compact, false) &&
       s !== nothing
      print(i, s)
      return
    end
  end )
end

function find_name(A, M = Main)
  for a = names(M)
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
    throw(error("coercion not possible"))
  else
    return nothing
  end
end

#to allow +(a::T, b::T) where a, b have different parents, but
# a common over structure
# designed(?) to be minimally invasive in AA and Nemo, but filled with
# content in Hecke/Oscar
function force_op(op::Function, throw_error::Type{Val{T}}, a...) where {T}
  if throw_error === Val{true}
    throw(error("no common overstructure for the arguments found"))
  end
  return false
end

function force_op(op::Function, a...)
  return force_op(op, Val{true}, a...)
end

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
import .PrettyPrinting: canonicalize
import .PrettyPrinting: expr_to_latex_string
import .PrettyPrinting: expr_to_string
import .PrettyPrinting: expressify
import .PrettyPrinting: get_html_as_latex
import .PrettyPrinting: get_syntactic_sign_abs
import .PrettyPrinting: is_syntactic_one
import .PrettyPrinting: is_syntactic_zero
import .PrettyPrinting: obj_to_latex_string
import .PrettyPrinting: obj_to_string
import .PrettyPrinting: obj_to_string_wrt_times
import .PrettyPrinting: print_integer_string
import .PrettyPrinting: print_obj
import .PrettyPrinting: printer
import .PrettyPrinting: set_html_as_latex
import .PrettyPrinting: show_obj
import .PrettyPrinting: show_via_expressify

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
include("FunctionField.jl")
include("Residue.jl")
include("ResidueField.jl")
include("NumberField.jl")
include("Fraction.jl")
include("TotalFraction.jl")
include("MPoly.jl")
include("UnivPoly.jl")
include("FreeAssAlgebra.jl")
include("LaurentMPoly.jl")

###############################################################################
#
#   Generic submodule
#
###############################################################################

include("Generic.jl")

# Do not import div, divrem, exp, inv, log, sqrt, numerator and denominator
# as we have our own
import .Generic: @perm_str
import .Generic: abs_series
import .Generic: abs_series_type
import .Generic: base_field
import .Generic: basis
import .Generic: character
import .Generic: check_composable
import .Generic: collength
import .Generic: combine_like_terms!
import .Generic: cycles
import .Generic: defining_polynomial
import .Generic: degrees
import .Generic: dense_matrix_type
import .Generic: dense_poly_type
import .Generic: dim
import .Generic: disable_cache!
import .Generic: downscale
import .Generic: enable_cache!
import .Generic: exp_gcd
import .Generic: exponent
import .Generic: exponent_vector
import .Generic: exponent_word
import .Generic: finish
import .Generic: fit!
import .Generic: gcd
import .Generic: gcdx
import .Generic: has_bottom_neighbor
import .Generic: has_left_neighbor
import .Generic: hash
import .Generic: hooklength
import .Generic: identity_map
import .Generic: image_fn
import .Generic: image_map
import .Generic: intersection
import .Generic: inv!
import .Generic: inverse_fn
import .Generic: inverse_image_fn
import .Generic: inverse_mat
import .Generic: invmod
import .Generic: is_compatible
import .Generic: is_degree
import .Generic: is_divisible_by
import .Generic: is_homogeneous
import .Generic: is_isomorphic
import .Generic: is_reverse
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
import .Generic: map_from_func
import .Generic: map_with_preimage_from_func
import .Generic: map_with_retraction
import .Generic: map_with_retraction_from_func
import .Generic: map_with_section
import .Generic: map_with_section_from_func
import .Generic: map1
import .Generic: map2
import .Generic: mat
import .Generic: matrix_repr
import .Generic: max_fields
import .Generic: mod
import .Generic: monomial
import .Generic: monomial_iszero
import .Generic: monomial_set!
import .Generic: monomial!
import .Generic: monomials
import .Generic: mpoly_type
import .Generic: MPolyBuildCtx
import .Generic: mullow_karatsuba
import .Generic: ngens
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
import .Generic: preimage
import .Generic: preimage_map
import .Generic: prime
import .Generic: push_term!
import .Generic: rand_ordering
import .Generic: reduce!
import .Generic: rel_series
import .Generic: rel_series_type
import .Generic: rels
import .Generic: rescale!
import .Generic: retraction_map
import .Generic: reverse
import .Generic: reverse_rows
import .Generic: reverse_rows!
import .Generic: right_kernel
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

# Do not export inv, div, divrem, exp, log, sqrt, numerator and denominator as we define our own
export @perm_str
export abs_series
export abs_series_type
export addmul_delayed_reduction!
export addmul!
export base_field
export base_ring
export basis
export can_solve_left_reduced_triu
export canonical_unit
export change_base_ring
export change_coefficient_ring
export character
export chebyshev_t
export chebyshev_u
export check_composable
export check_parent
export collength
export combine_like_terms!
export cycles
export defining_polynomial
export dense_matrix_type
export dense_poly_type
export det
export discriminant
export elem_type
export exponent
export exponent_vector
export exponent_word
export finish
export fit!
export gcd
export gcd_with_cofactors
export gcdinv
export gcdx
export gen
export gens
export has_bottom_neighbor
export has_left_neighbor
export hash
export interpolate
export intersection
export inv!
export inverse_image_fn
export inverse_mat
export invmod
export is_compatible
export is_degree
export is_divisible_by
export is_domain_type
export is_exact_type
export is_gen
export is_homogeneous
export is_isomorphic
export is_monomial
export is_monomial_recursive
export is_negative
export is_reverse
export is_submodule
export is_symmetric
export is_term_recursive
export is_unit
export is_zero_divisor
export is_zero_divisor_with_annihilator
export laurent_series
export lcm
export length
export main_variable
export main_variable_extract
export main_variable_insert
export mat
export matrix_repr
export MatrixElem
export max_fields
export mod
export monomial
export monomial_iszero
export monomial_set!
export monomial_to_newton!
export monomial!
export monomials
export mpoly_type
export MPolyBuildCtx
export mul_ks
export mul_red!
export mullow_karatsuba
export mulmod
export newton_to_monomial!
export ngens
export normal_form
export normalise
export nullspace
export num_coeff
export one
export order
export ordering
export parent_type
export parity
export partitionseq
export perm
export permtype
export polcoeff
export poly
export poly_ring
export polynomial
export PolynomialElem
export pow_multinomial
export ppio
export precision
export preimage
export prime
export push_term!
export rand_ordering
export rank
export reduce!
export rel_series
export rel_series_type
export rels
export renormalize!
export resultant
export resultant_ducos
export resultant_euclidean
export resultant_subresultant
export resultant_sylvester
export resx
export reverse
export rowlength
export set_exponent_vector!
export set_exponent_word!
export setcoeff!
export setpermstyle
export size
export sort_terms!
export subst
export summands
export supermodule
export sylvester_matrix
export term
export terms
export total_degree
export trailing_coefficient
export truncate
export weights
export zero

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

function YoungTableau(part::Generic.Partition, tab::Matrix{Int})
   Generic.YoungTableau(part, tab)
end

function YoungTableau(part::Generic.Partition, fill::Vector{Int}=collect(1:part.n))
   Generic.YoungTableau(part, fill)
end

function number_field(a::Generic.Poly{Rational{BigInt}}, s::AbstractString, t = "\$"; cached = true)
   return Generic.number_field(a, Symbol(s), t; cached=cached)
end

function number_field(a::Generic.Poly{Rational{BigInt}}, s::Char, t = "\$"; cached = true)
   return Generic.number_field(a, Symbol(s), t; cached=cached)
end

function number_field(a::Generic.Poly{Rational{BigInt}}, s::Symbol, t = "\$"; cached = true)
   return Generic.number_field(a, s, t; cached=cached)
end

function FunctionField(p::Generic.Poly{Generic.RationalFunctionFieldElem{T, U}}, s::Symbol; cached::Bool=true) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   return Generic.FunctionField(p, s; cached=cached)
end

function FunctionField(p::Generic.Poly{Generic.RationalFunctionFieldElem{T, U}}, s::AbstractString; cached::Bool=true) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   return Generic.FunctionField(p, Symbol(s); cached=cached)
end

function FunctionField(p::Generic.Poly{Generic.RationalFunctionFieldElem{T, U}}, s::Char; cached::Bool=true) where {T <: FieldElement, U <: Union{PolyRingElem, MPolyRingElem}}
   return Generic.FunctionField(p, Symbol(s); cached=cached)
end

@doc Markdown.doc"""
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
#   Polynomial Ring S, x = R["x"] syntax
#
###############################################################################

getindex(R::NCRing, s::Union{String, Char, Symbol}) = polynomial_ring(R, s)
getindex(R::NCRing, s::Union{String, Char, Symbol}, ss::Union{String, Char}...) =
   polynomial_ring(R, [s, ss...])

# syntax x = R["x"]["y"]
getindex(R::Tuple{Union{Ring, NCRing}, Union{PolyRingElem, NCPolyRingElem}}, s::Union{String, Char, Symbol}) = polynomial_ring(R[1], s)

###############################################################################
#
#   Matrix M = R[...] syntax
#
################################################################################

typed_hvcat(R::Ring, dims::Dims, d...) = matrix(R, length(dims), dims[1], hvcat(dims, d...))
typed_hcat(R::Ring, d...) = matrix(R, 1, length(d), hcat(d...))
typed_vcat(R::Ring, d...) = matrix(R, length(d), 1, vcat(d...))

###############################################################################
#
#   misc
#
###############################################################################

include("misc/ProductIterator.jl")

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

################################################################################
#
#   Deprecations
#
################################################################################

include("Deprecations.jl")

###############################################################################
#
#   Array creation functions
#
###############################################################################

Array(R::Ring, r::Int...) = Array{elem_type(R)}(undef, r)

function zeros(R::Ring, r::Int...)
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
include("algorithms/DensePoly.jl")

###############################################################################
#
#  For backwards compatibility
#
###############################################################################

function is_negative end

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
