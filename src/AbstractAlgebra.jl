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
                  :RealField, :NumberField,
                  :AbstractAlgebra,
                  :inv, :log, :exp, :sqrt, :div, :divrem,
                  :numerator, :denominator,
                  :promote_rule,
                  :Set, :Module, :Ring, :Group, :Field]

# If you want to add methods to functions in LinearAlgebra they should be
# imported here and in Generic.jl, and exported below.
# They should not be imported/exported anywhere else.

import LinearAlgebra: det, ishermitian, issymmetric, istriu, norm, nullspace, rank,
                      hessenberg

import LinearAlgebra: lu, lu!, tr

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

import Base: Array, abs, acos, acosh, asin, asinh, atan, atanh, axes,
             bin, ceil, checkbounds, conj, convert, cmp, cos, cosh, cospi, cot,
             coth, dec, deepcopy, deepcopy_internal, expm1, exponent, fill,
             floor, gcd, gcdx, getindex, hash, hcat, hex, hypot, intersect,
             invmod, isequal, isfinite, isless, isone, isqrt, isreal,
             iszero, lcm, ldexp, length, log1p, mod, ndigits, oct, one,
             parent, parse, powermod,
             precision, rand, Rational, rem, reverse, setindex!,
             show, sincos, similar, sign, sin, sinh, sinpi, size, string, tan,
             tanh, trailing_zeros, transpose, truncate, typed_hvcat,
             typed_hcat, typed_vcat, vcat, xor, zero, zeros,
             +, -, *, ==, ^, &, |, <<, >>, ~, <=, >=, <, >, //, /, !=

using Random: Random, AbstractRNG, SamplerTrivial
using RandomExtensions: RandomExtensions, make, Make2

export elem_type, parent_type

export Field

export SetElem, GroupElem, AdditiveGroupElem, NCRingElem, RingElem, ModuleElem, FieldElem, RingElement,
       FieldElement, Map

export SetMap, FunctionalMap, IdentityMap

export NCPolyElem, PolyElem, SeriesElem, AbsSeriesElem, RelSeriesElem, ResElem, FracElem,
       MatElem, MatAlgElem, FinFieldElem, MPolyElem, UnivPolyElem, NumFieldElem, Ideal,
       SimpleNumFieldElem, FreeAssAlgElem

export PolyRing, SeriesRing, ResRing, FracField, MatSpace, MatAlgebra,
       FinField, MPolyRing, UnivPolyRing, NumField, SimpleNumField, IdealSet,
       FreeAssAlgebra

export ZZ, QQ, zz, qq, RealField, RDF

export zeros

export NotImplementedError, NotInvertibleError, ErrorConstrDimMismatch,
       error_dim_negative

export add!, addeq!,
       crt,
       factor, factor_squarefree,
       hgcd,
       is_irreducible, is_squarefree,
       mul!,
       sub!,
       zero!

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

const PolynomialElem{T} = Union{PolyElem{T}, NCPolyElem{T}}
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

import .PrettyPrinting: get_html_as_latex, set_html_as_latex, expressify,
                        show_via_expressify, @enable_all_show_via_expressify,
                        expr_to_string, expr_to_latex_string, canonicalize,
                        printer, print_integer_string, get_syntactic_sign_abs,
                        is_syntactic_one, is_syntactic_zero,
                        obj_to_string, obj_to_string_wrt_times,
                        obj_to_latex_string, show_obj, print_obj

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
import .Generic: abs_series, abs_series_type,
                 base_field, basis,
                 character,
                 check_composable, collength,
                 combine_like_terms!, cycles,
                 defining_polynomial, degrees,
                 dense_matrix_type, dense_poly_type,
                 dim, disable_cache!,
                 downscale,
                 enable_cache!, exp_gcd,
                 exponent, exponent_vector, exponent_word,
                 finish, fit!, gcd, gcdx,
                 has_left_neighbor, has_bottom_neighbor, hash,
                 hooklength, identity_map,
                 image_map, image_fn, intersection,
                 inverse_fn, inverse_image_fn,
                 inverse_mat, reverse_rows, reverse_rows!,
                 inv!, invmod,
                 is_compatible, is_degree, is_divisible_by,
                 is_homogeneous, is_isomorphic,
                 isone, is_reverse, is_rimhook, is_submodule,
                 is_unit,
                 laurent_ring, laurent_series, lcm,
                 leading_coefficient, leading_monomial,
                 leading_exponent_vector,
                 leading_term, leading_exponent_word, length,
                 leglength, main_variable,
                 main_variable_extract, main_variable_insert,
                 map1, map2, map_from_func,
                 map_with_preimage_from_func, map_with_retraction,
                 map_with_retraction_from_func,
                 map_with_section, map_with_section_from_func, mat,
                 matrix_repr, max_fields, mod,
                 monomial, monomial!, monomials,
                 monomial_iszero, monomial_set!,
                 mpoly_type, MPolyBuildCtx, mullow_karatsuba,
                 ngens, norm, normal_form, normalise,
                 num_coeff, one,
                 order, ordering, parity, partitionseq, Perm, perm,
                 permtype, @perm_str, polcoeff, poly, poly_ring,
                 precision, preimage, preimage_map,
                 prime, push_term!,
                 rand_ordering, reduce!,
                 rels, rel_series, rel_series_type,
                 rescale!, retraction_map, reverse,
                 right_kernel, rowlength, section_map, setcoeff!,
                 set_exponent_vector!, set_exponent_word!, set_limit!,
                 setpermstyle, size,
                 sort_terms!, summands,
                 supermodule, term, terms, total_degree,
                 to_univariate, trailing_coefficient, truncate,
                 unit, upscale, weights,
                 to_univariate, trailing_coefficient, truncate,
                 unit, upscale,
                 zero,
       # Moved from Hecke into Misc
                 Loc, Localization, LocElem,
                 roots, sturm_sequence

# Do not export inv, div, divrem, exp, log, sqrt, numerator and denominator as we define our own
export abs_series, abs_series_type,
                 addmul_delayed_reduction!, addmul!,
                 base_field, base_ring, basis,
                 canonical_unit, can_solve_left_reduced_triu,
                 change_base_ring, change_coefficient_ring, character,
                 chebyshev_t,
                 chebyshev_u, check_composable, check_parent,
                 collength, combine_like_terms!, cycles,
                 defining_polynomial,
                 dense_matrix_type, dense_poly_type, det,
                 discriminant,
                 elem_type,
                 exponent, exponent_vector, exponent_word,
                 finish, fit!, gcd, gen,
                 gens, gcdinv, gcdx,
                 has_left_neighbor, has_bottom_neighbor, hash,
                 interpolate, intersection,
                 inv!, inverse_image_fn,
                 inverse_mat, invmod,
                 is_compatible, is_degree, is_divisible_by,
                 is_domain_type, is_exact_type, is_gen,
                 is_homogeneous,
                 is_isomorphic, is_monomial, is_monomial_recursive,
                 is_negative, is_reverse,
                 is_submodule, is_symmetric,
                 is_term_recursive, is_unit,
                 lcm,
                 laurent_series, length,
                 main_variable, main_variable_extract, main_variable_insert,
                 mat, matrix_repr, max_fields, mod,
                 monomial, monomial!, monomials,
                 monomial_iszero, monomial_set!, monomial_to_newton!,
                 mpoly_type, MPolyBuildCtx,
                 mul_ks, mul_red!, mullow_karatsuba, mulmod,
                 newton_to_monomial!, ngens,
                 normal_form, normalise, nullspace, num_coeff,
                 one, order, ordering,
                 parent_type, parity, partitionseq, Perm, perm, permtype,
                 @perm_str, polcoeff, polynomial, poly,
                 poly_ring, pow_multinomial,
                 ppio, precision, preimage, prime,
                 push_term!, rank,
                 rand_ordering, reduce!,
                 renormalize!, rel_series, rel_series_type, rels,
                 resultant, resultant_ducos,
                 resultant_euclidean, resultant_subresultant,
                 resultant_sylvester, resx, reverse, rowlength,
                 setcoeff!, set_exponent_vector!, set_exponent_word!,
                 setpermstyle,
                 size, sort_terms!, subst, summands, supermodule,
                 sylvester_matrix, term, terms,
                 total_degree, trailing_coefficient, truncate, weights,
                 zero,
                 MatrixElem, PolynomialElem,
       # Moved from Hecke into Misc
                 divexact_low, divhigh,
                 is_monic, Loc, Localization, LocElem, mulhigh_n,
                 PolyCoeffs, roots, sturm_sequence

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

function NumberField(a::Generic.Poly{Rational{BigInt}}, s::AbstractString, t = "\$"; cached = true)
   return Generic.NumberField(a, Symbol(s), t; cached=cached)
end

function NumberField(a::Generic.Poly{Rational{BigInt}}, s::Char, t = "\$"; cached = true)
   return Generic.NumberField(a, Symbol(s), t; cached=cached)
end

function NumberField(a::Generic.Poly{Rational{BigInt}}, s::Symbol, t = "\$"; cached = true)
   return Generic.NumberField(a, s, t; cached=cached)
end

function FunctionField(p::Generic.Poly{Generic.Rat{T, U}}, s::Symbol; cached::Bool=true) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   return Generic.FunctionField(p, s; cached=cached)
end

function FunctionField(p::Generic.Poly{Generic.Rat{T, U}}, s::AbstractString; cached::Bool=true) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
   return Generic.FunctionField(p, Symbol(s); cached=cached)
end

function FunctionField(p::Generic.Poly{Generic.Rat{T, U}}, s::Char; cached::Bool=true) where {T <: FieldElement, U <: Union{PolyElem, MPolyElem}}
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

getindex(R::NCRing, s::Union{String, Char, Symbol}) = PolynomialRing(R, s)
getindex(R::NCRing, s::Union{String, Char, Symbol}, ss::Union{String, Char}...) =
   PolynomialRing(R, [s, ss...])

# syntax x = R["x"]["y"]
getindex(R::Tuple{Union{Ring, NCRing}, Union{PolyElem, NCPolyElem}}, s::Union{String, Char, Symbol}) = PolynomialRing(R[1], s)

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

include("algorithms/MPolyFactor.jl")
include("algorithms/DensePoly.jl")

###############################################################################
#
#  For backwards compability
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
