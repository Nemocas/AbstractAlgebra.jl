@doc Markdown.doc"""
 AbstractAlgebra is a pure Julia package for computational abstract algebra.

 For more information see https://github.com/Nemocas/AbstractAlgebra.jl
"""
module AbstractAlgebra

using RandomExtensions: RandomExtensions, make, Make, Make2, Make3, Make4

using Markdown

using InteractiveUtils

using Test # for "interface-conformance" functions

# A list of all symbols external packages should not import from AbstractAlgebra
import_exclude = [:import_exclude, :QQ, :ZZ,
                  :RealField, :NumberField,
                  :AbstractAlgebra,
                  :exp, :sqrt, :div, :divrem, :numerator, :denominator,
                  :promote_rule,
                  :Set, :Module, :Ring, :Group, :Field]

# If you want to add methods to functions in LinearAlgebra they should be
# imported here and in Generic.jl, and exported below.
# They should not be imported/exported anywhere else.

import LinearAlgebra: det, issymmetric, istriu, norm, nullspace, rank,
                      transpose!, hessenberg

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
const Base_import_exclude = [:exp, :log, :sqrt, :div, :divrem, :numerator,
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

function sqrt(a::T) where T
  return Base.sqrt(a)
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

import Base: Array, abs, acos, acosh, adjoint, asin, asinh, atan, atanh, axes,
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

export SetElem, GroupElem, AdditiveGroupElem, NCRingElem, RingElem, ModuleElem, FieldElem, RingElement,
       FieldElement, Map, AccessorNotSetError

export SetMap, FunctionalMap, IdentityMap

export NCPolyElem, PolyElem, SeriesElem, AbsSeriesElem, RelSeriesElem, ResElem, FracElem,
       MatElem, MatAlgElem, FinFieldElem, MPolyElem, NumFieldElem, SimpleNumFieldElem

export PolyRing, SeriesRing, ResRing, FracField, MatSpace, MatAlgebra,
       FinField, MPolyRing, NumField, SimpleNumField

export ZZ, QQ, zz, qq, RealField, RDF

export SymmetricGroup

export create_accessors, get_handle, package_handle, zeros,
       Array, sig_exists

export error_dim_negative, ErrorConstrDimMismatch

export crt, factor, factor_squarefree

function expressify
end

function show_via_expressify
end

###############################################################################
# Macros for fancy printing and extending objects when desired
# fancy printing (and extending)
# change
#
# struct bla..
# ..
# end
#
# to
#
# struct bla ..
# @declare_other
# ...
# end
#
# Then, in the show function, start with
# @show_name(io, obj)
# If the user assigned a name to the object (in the REPL mainly) by doing
# A = bla...
# then, in the compact printing only the name "A" is printed
# also adding
# @show_special(io, obj)
# allows, if present to call a different printing function for this instance
# See FreeModule for an example
#
###############################################################################

macro declare_other()
   esc(quote other::Dict{Symbol, Any} end )
end

function set_name!(G, name::String)
   set_special(G, :name => name)
end

function hasspecial(G)
   if !isdefined(G, :other)
      return false, nothing
   else
     return true, G.other
   end
end

function get_special(G, s::Symbol)
   fl, D = hasspecial(G)
   fl && return get(D, s, nothing)
   nothing
end

function set_name!(G)
   s = get_special(G, :name)
   s === nothing || return
   sy = find_name(G)
   sy === nothing && return
   set_name!(G, string(sy))
end

function set_special(G, data::Pair{Symbol, <:Any}...)
  if !isdefined(G, :other)
    D = G.other = Dict{Symbol, Any}()
  else
    D = G.other
  end

  for d in data
    push!(D, d)
  end
end

extra_name(G) = nothing

macro show_name(io, O)
  return :( begin
    local i = $(esc(io))
    local o = $(esc(O))
    s = get_special(o, :name)
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
  for a = names(Main)
    a === :ans && continue
    d = Meta.parse("$M.$a")
    try
      z = eval(d);
      if z === A
        return a
      end
    catch e
    end
  end
end

macro show_special(io, O)
  return :( begin
    local i = $(esc(io))
    local o = $(esc(O))
    s = get_special(o, :show)
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
    s = get_special(o, :show_elem)
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

@static if false #VERSION >= v"1.6"
  const CacheDictType = WeakValueDict
else
  const CacheDictType = Dict
end

function get_cached!(default::Base.Callable, dict::AbstractDict,
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
#  Version information
#
################################################################################

function versioninfo()
  print("AbstractAlgebra version 0.3.0\n")
  abstractalgebrarepo = dirname(dirname(@__FILE__))

  print("AbstractAlgebra: ")
  prepo = Base.LibGit2.GitRepo(abstractalgebrarepo)
  Base.LibGit2.with(LibGit2.head(prepo)) do phead
    print("commit: ")
    print(string(LibGit2.Oid(phead))[1:8])
    print(" date: ")
    commit = Base.LibGit2.get(Base.LibGit2.GitCommit, prepo, LibGit2.Oid(phead))
    print(Base.Dates.unix2datetime(Base.LibGit2.author(commit).time))
    print(")\n")
  end

  finalize(prepo)

  return nothing
end

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

###############################################################################
#
#   Generic submodule
#
###############################################################################

include("Generic.jl")

# Do not import div, divrem, exp, inv, log, sqrt, numerator and denominator
# as we have our own
import .Generic: abs_series, abs_series_type, add!, addeq!, 
                 base_field, basis, cached,
                 character,
                 check_composable,
                 codomain, coeff, combine_like_terms!, cycles,
                 data, defining_polynomial, degrees,
                 dense_matrix_type, dense_poly_type,
                 dim, disable_cache!,
                 domain, downscale,
                 enable_cache!, exp_gcd,
                 exponent, exponent_vector, exponent_vectors,
                 finish, fit!, gcd,
                 get_field, gcdx,
                 has_left_neighbor, has_bottom_neighbor, hash,
                 hooklength, identity_map, image,
                 image_map, image_fn,
                 invariant_factors,
                 inverse_fn, inverse_image_fn,
                 inverse_mat, reverse_rows, reverse_rows!,
                 inv!, invmod,
                 iscompatible, isconstant, isdegree,
                 ishomogeneous, isisomorphic,
                 isone, isreverse, isrimhook, issubmodule,
                 isunit, isunivariate,
                 laurent_ring, lcm, leading_coefficient, leading_monomial,
		 leading_term, length,
                 leglength, main_variable,
                 main_variable_extract, main_variable_insert,
                 map1, map2, map_from_func,
                 map_with_preimage_from_func, map_with_retraction,
                 map_with_retraction_from_func,
                 map_with_section, map_with_section_from_func, mat,
                 matrix_repr, max_fields, max_precision, mod,
                 monomial, monomial!, monomials,
                 monomial_iszero, monomial_set!,
                 MPolyBuildCtx, mul!, mullow_karatsuba,
                 ngens, norm, normalise,
		 num_coeff, O, one,
                 order, ordering, parity, partitionseq, Perm, perm,
                 permtype, @perm_str, polcoeff, pol_length,
                 precision, preimage, preimage_map,
		 prime, push_term!,
                 rand_ordering, reduce!,
                 renormalize!, rels, rel_series, rel_series_type,
		 rescale!, retraction_map, reverse,
                 right_kernel, section_map, setcoeff!,
                 set_exponent_vector!, set_field!, set_length!, set_limit!,
                 setpermstyle, set_precision!, set_valuation!, size,
                 sort_terms!, sub, summands,
                 supermodule, term, terms, total_degree,
                 to_univariate, trailing_coefficient,
		 truncate, upscale, var_index, vars,
                 zero,
                 @PolynomialRing, MatrixElem,
       # Moved from Hecke into Misc
		 Loc, Localization, LocElem,
		 polynomial_to_power_sums,
		 power_sums_to_polynomial, roots, sturm_sequence

# Do not export inv, div, divrem, exp, log, sqrt, numerator and denominator as we define our own
export abs_series, abs_series_type, add!, addeq!,
                 addmul_delayed_reduction!, addmul!,
		 base_field, base_ring, cached,
                 canonical_unit, can_solve_left_reduced_triu,
                 change_base_ring, character,
                 chebyshev_t,
                 chebyshev_u, check_composable, check_parent,
                 codomain, coeff, combine_like_terms!, cycles,
                 data, defining_polynomial, degrees,
                 dense_matrix_type, dense_poly_type, det,
                 dim, disable_cache!,
                 discriminant, domain, downscale,
                 elem_type, enable_cache!, exp_gcd,
                 exponent, exponent_vector, exponent_vectors,
                 finish, fit!, gcd, gen,
                 gens, get_field, gcdinv, gcdx,
                 has_left_neighbor, has_bottom_neighbor, hash,
                 hooklength, identity_map, image,
                 image_map, image_fn, interpolate,
                 inv!, invariant_factors,
                 inverse_fn, inverse_image_fn,
                 inverse_mat, invmod,
                 iscompatible, isconstant, isdegree,
                 isdomain_type, isexact_type, isgen,
                 ishomogeneous,
                 isisomorphic, ismonomial, ismonomial_recursive,
                 isnegative, isone, isreverse,
                 isrimhook, issubmodule, issymmetric,
                 isterm_recursive, isunit, isunivariate, iszero, laurent_ring,
                 lcm, leading_coefficient, leading_monomial, leading_term,
		 leglength, length,
                 main_variable, main_variable_extract, main_variable_insert,
                 map1, map2, map_from_func,
		 map_with_preimage_from_func,
                 map_with_retraction, map_with_retraction_from_func,
                 map_with_section, map_with_section_from_func,
                 mat, matrix_repr, max_fields,
                 max_precision, mod,
                 monomial, monomial!, monomials,
                 monomial_iszero, monomial_set!, monomial_to_newton!,
                 MPolyBuildCtx, mul!,
                 mul_ks, mul_red!, mullow_karatsuba, mulmod,
                 needs_parentheses, newton_to_monomial!, ngens, norm,
                 normalise, nullspace, num_coeff,
		 O, one, order, ordering,
                 parent_type, parity, partitionseq, Perm, perm, permtype,
                 @perm_str, polcoeff, pol_length, polynomial, pow_multinomial,
                 ppio, precision, preimage,
                 preimage_map, prime,
		 push_term!, rank,
                 rand_ordering, reduce!,
                 renormalize!, rel_series, rel_series_type, rels,
		 resultant, resultant_ducos, rescale!,
                 resultant_euclidean, resultant_subresultant,
                 resultant_sylvester, resx, retraction_map, reverse,
                 section_map, setcoeff!,
                 set_exponent_vector!, set_field!, set_length!, set_limit!,
                 setpermstyle, set_precision!, set_valuation!,
                 size, sort_terms!, sub, subst, summands, supermodule,
                 sylvester_matrix, term, terms, to_univariate,
                 total_degree, trailing_coefficient, truncate,
                 upscale, var_index, vars, zero,
                 @PolynomialRing, MatrixElem,
       # Moved from Hecke into Misc
                 divexact_low, divhigh,
		 ismonic, Loc, Localization, LocElem, mulhigh_n,
		 polynomial_to_power_sums, PolyCoeffs,
		 power_sums_to_polynomial, roots, sturm_sequence

# TODO remove these two once removed from dependent packages (Hecke)
export displayed_with_minus_in_front, show_minus_one

################################################################################
#
#   Parent constructors
#
################################################################################

function SymmetricGroup(n::T) where T
  Generic.SymmetricGroup(n)
end

function AllPerms(n::T) where T
  Generic.AllPerms(n)
end

function Partition(part::AbstractVector{T}, check::Bool=true) where T
  Generic.Partition(part, check)
end

function AllParts(n::T) where T
  Generic.AllParts(n)
end

function SkewDiagram(lambda::Generic.Partition, mu::Generic.Partition)
  Generic.SkewDiagram(lambda, mu)
end

function SkewDiagram(lambda::Vector{T}, mu::Vector{T}) where T
  Generic.SkewDiagram(lambda, mu)
end

function YoungTableau(part::Generic.Partition, tab::Array{Int, 2})
   Generic.YoungTableau(part, tab)
end

function YoungTableau(part::Generic.Partition, fill::Vector{Int}=collect(1:part.n))
   Generic.YoungTableau(part, fill)
end

function YoungTableau(p::Vector{Int})
   Generic.YoungTableau(p)
end

@doc (@doc Generic.PowerSeriesRing)
PowerSeriesRing(R::Ring, prec::Int, s::Union{Char, AbstractString}; cached=true, model=:capped_relative)

function PowerSeriesRing(R::Ring, prec::Int, s::AbstractString; cached=true, model=:capped_relative)
   Generic.PowerSeriesRing(R, prec, s; cached=cached, model=model)
end

function PowerSeriesRing(R::Ring, prec::Int, s::Char; cached=true, model=:capped_relative)
   PowerSeriesRing(R, prec, string(s); cached=cached, model=model)
end

function PowerSeriesRing(R::AbstractAlgebra.Ring, prec::Vector{Int}, s::Vector{T}; cached=true, model=:capped_absolute) where T <: AbstractString
   Generic.PowerSeriesRing(R, prec, s; cached=cached, model=model)
end

function PowerSeriesRing(R::AbstractAlgebra.Ring, prec::Int, s::Vector{T}; cached=true, model=:capped_absolute) where T <: AbstractString
   Generic.PowerSeriesRing(R, prec, s; cached=cached, model=model)
end

@doc (@doc Generic.LaurentSeriesRing)
LaurentSeriesRing(R::Ring, prec::Int, s::Union{Char, AbstractString}; cached=true)

function LaurentSeriesRing(R::Ring, prec::Int, s::AbstractString; cached=true)
   Generic.LaurentSeriesRing(R, prec, s; cached=cached)
end

function LaurentSeriesRing(R::Ring, prec::Int, s::Char; cached=true)
   LaurentSeriesRing(R, prec, string(s); cached=cached)
end

function LaurentSeriesRing(R::Field, prec::Int, s::AbstractString; cached=true)
   Generic.LaurentSeriesField(R, prec, s; cached=cached)
end

function LaurentSeriesRing(R::Field, prec::Int, s::Char; cached=true)
   LaurentSeriesField(R, prec, string(s); cached=cached)
end

function LaurentSeriesField(R::Field, prec::Int, s::AbstractString; cached=true)
   Generic.LaurentSeriesField(R, prec, s; cached=cached)
end

function LaurentSeriesField(R::Field, prec::Int, s::Char; cached=true)
   LaurentSeriesField(R, prec, string(s); cached=cached)
end

function PuiseuxSeriesRing(R::Ring, prec::Int, s::AbstractString; cached=true)
   Generic.PuiseuxSeriesRing(R, prec, s; cached=cached)
end

function PuiseuxSeriesRing(R::Ring, prec::Int, s::Char; cached=true)
   PuiseuxSeriesRing(R, prec, string(s); cached=cached)
end

function PuiseuxSeriesRing(R::Field, prec::Int, s::AbstractString; cached=true)
   Generic.PuiseuxSeriesField(R, prec, s; cached=cached)
end

function PuiseuxSeriesRing(R::Field, prec::Int, s::Char; cached=true)
   PuiseuxSeriesField(R, prec, string(s); cached=cached)
end

function PuiseuxSeriesField(R::Field, prec::Int, s::AbstractString; cached=true)
   Generic.PuiseuxSeriesField(R, prec, s; cached=cached)
end

function PuiseuxSeriesField(R::Field, prec::Int, s::Char; cached=true)
   PuiseuxSeriesField(R, prec, string(s); cached=cached)
end

@doc Markdown.doc"""
    PolynomialRing(R::AbstractAlgebra.Ring, s::Vector{T}; cached::Bool = true, ordering::Symbol = :lex) where T <: Union{String, Char}

Given a base ring `R` and an array of strings `s` specifying how the
generators (variables) should be printed, return a tuple `T, (x1, x2, ...)`
representing the new polynomial ring $T = R[x1, x2, ...]$ and the generators
$x1, x2, ...$ of the polynomial ring. By default the parent object `T` will
depend only on `R` and `x1, x2, ...` and will be cached. Setting the optional
argument `cached` to `false` will prevent the parent object `T` from being
cached. `S` is a symbol corresponding to the ordering of the polynomial and
can be one of `:lex`, `:deglex` or `:degrevlex`.
"""
PolynomialRing(R::Ring, s::Union{Vector{String}, Vector{Char}}; cached::Bool = true, ordering::Symbol = :lex)

function PolynomialRing(R::Ring, s::Array{String, 1}; cached::Bool = true, ordering::Symbol = :lex)
   Generic.PolynomialRing(R, s; cached=cached, ordering=ordering)
end

function PolynomialRing(R::Ring, s::Array{Char, 1}; cached::Bool = true, ordering::Symbol = :lex)
   PolynomialRing(R, string.(s); cached=cached, ordering=ordering)
end

function PolynomialRing(R::AbstractAlgebra.Ring, n::Int, s::String="x";
                                 cached::Bool = false, ordering::Symbol = :lex)
   PolynomialRing(R, ["$s$i" for i=1:n]; cached = cached, ordering = ordering)
end

function PolynomialRing(R::AbstractAlgebra.Ring, n::Int, s::Char;
                                 cached::Bool = false, ordering::Symbol = :lex)
   PolynomialRing(R, n, string(s); cached = cached, ordering = ordering)
end

@doc (@doc Generic.SparsePolynomialRing)
SparsePolynomialRing(R::Ring, s::Union{Char, String}; cached::Bool = true)

function SparsePolynomialRing(R::Ring, s::String; cached::Bool = true)
   Generic.SparsePolynomialRing(R, s; cached=cached)
end

function SparsePolynomialRing(R::Ring, s::Char; cached::Bool = true)
   SparsePolynomialRing(R, string(s); cached=cached)
end

@doc (@doc Generic.LaurentPolynomialRing)
LaurentPolynomialRing(R::Ring, s::AbstractString) = Generic.LaurentPolynomialRing(R, s)

@doc (@doc Generic.FractionField)
FractionField(R::Ring; cached=true)

function FractionField(R::Ring; cached=true)
   Generic.FractionField(R; cached=cached)
end

@doc (@doc Generic.ResidueRing)
ResidueRing(R::Ring, a::Union{RingElement, Integer}; cached::Bool = true)

function ResidueRing(R::Ring, a::Union{RingElement, Integer}; cached::Bool = true)
   Generic.ResidueRing(R, a; cached=cached)
end

@doc (@doc Generic.ResidueField)
ResidueField(R::Ring, a::Union{RingElement, Integer}; cached::Bool = true)

function ResidueField(R::Ring, a::Union{RingElement, Integer}; cached::Bool = true)
   Generic.ResidueField(R, a; cached=cached)
end

function NumberField(a::AbstractAlgebra.Generic.Poly{Rational{BigInt}}, s::AbstractString, t = "\$"; cached = true)
   Generic.NumberField(a, s, t; cached=cached)
end

function NumberField(a::AbstractAlgebra.Generic.Poly{Rational{BigInt}}, s::Char, t = "\$"; cached = true)
   NumberField(a, string(s), t; cached=cached)
end

function RationalFunctionField(k::Field, s::AbstractString; cached = true)
   Generic.RationalFunctionField(k, s; cached=cached)
end

function FunctionField(p::Generic.Poly{Generic.Rat{T}}, s::AbstractString; cached::Bool=true) where T <: FieldElement
   Generic.FunctionField(p, s; cached=cached)
end

@doc Markdown.doc"""
    FreeModule(R::NCRing, rank::Int; cached::Bool = true)

Return the free module over the ring $R$ with the given rank.
"""
function FreeModule(R::NCRing, rank::Int; cached::Bool = true)
   Generic.FreeModule(R, rank; cached=cached)
end

function free_module(R::NCRing, rank::Int; cached::Bool = true)
   Generic.FreeModule(R, rank; cached=cached)
end

@doc Markdown.doc"""
    VectorSpace(R::Field, dim::Int; cached::Bool = true)

Return the vector space over the field $R$ with the given dimension.
"""
function VectorSpace(R::Field, dim::Int; cached::Bool = true)
   Generic.FreeModule(R, dim; cached=cached)
end

function vector_space(R::Field, dim::Int; cached::Bool = true)
   Generic.FreeModule(R, dim; cached=cached)
end

@doc Markdown.doc"""
    sub(m::Module{T}, gens::Vector{<:ModuleElem{T}}) where T <: RingElement

Return the submodule `S` of the module `m` generated by the given generators,
given as elements of `m`, and a map which is the canonical injection from `S`
to `m`.
"""
function sub(m::Module{T}, gens::Vector{<:ModuleElem{T}}) where T <: RingElement
   Generic.sub(m, gens)
end

# Handles empty vector of elements
function sub(m::Module{T}, gens::Vector{Any}) where T <: RingElement
   Generic.sub(m, gens)
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

@doc Markdown.doc"""
    quo(m::Module{T}, N::Module{T}) where T <: RingElement

Return the quotient `Q` of the module `m` by the submodule `N` of `m`, and a
map which is a lift of elements of `Q` to `m`.
"""
function quo(m::Module{T}, N::Module{T}) where T <: RingElement
   Generic.quo(m, N)
end

@doc Markdown.doc"""
    DirectSum(m::Vector{<:Module{T}}) where T <: RingElement

Return a tuple $M, f, g$ consisting of $M$ the direct sum of the modules `m`
(supplied as a vector of modules), a vector $f$ of the canonical injections
of the $m[i]$ into $M$ and a vector $g$ of the canonical projections from
$M$ onto the $m[i]$.
"""
function DirectSum(m::Vector{<:Module{T}}) where T <: RingElement
   Generic.DirectSum(m)
end

function direct_sum(m::Vector{<:Module{T}}) where T <: RingElement
   Generic.DirectSum(m)
end

@doc Markdown.doc"""
    DirectSum(m::Module{T}...) where T <: RingElement

Return a tuple $M, f, g$ consisting of $M$ the direct sum of the given
modules, a vector $f$ of the canonical injections of the $m[i]$ into $M$
and a vector $g$ of the canonical projections from $M$ onto the $m[i]$.
"""
function DirectSum(m::Module{T}...) where T <: RingElement
   Generic.DirectSum(m...)
end

function direct_sum(m::Module{T}...) where T <: RingElement
   Generic.DirectSum(m...)
end

function ModuleHomomorphism(M1::AbstractAlgebra.Module, M2::AbstractAlgebra.Module, A...)
   Generic.ModuleHomomorphism(M1, M2, A...)
end

function module_homomorphism(M1::AbstractAlgebra.Module, M2::AbstractAlgebra.Module, m::MatElem)
   Generic.ModuleHomomorphism(M1, M2, m)
end

function ModuleIsomorphism(M1::AbstractAlgebra.Module, M2::AbstractAlgebra.Module, m::MatElem)
   Generic.ModuleIsomorphism(M1, M2, m)
end

function module_isomorphism(M1::AbstractAlgebra.Module, M2::AbstractAlgebra.Module, m::MatElem)
   Generic.ModuleIsomorphism(M1, M2, m)
end

export PowerSeriesRing, PolynomialRing, SparsePolynomialRing, LaurentPolynomialRing,
       MatrixSpace, MatrixAlgebra, FractionField, ResidueRing, Partition, SymmetricGroup,
       YoungTableau, AllParts, SkewDiagram, AllPerms, Perm, LaurentSeriesRing,
       LaurentSeriesField, ResidueField, NumberField, RationalFunctionField,
       FunctionField, PuiseuxSeriesRing,
       PuiseuxSeriesField, FreeModule, VectorSpace, ModuleHomomorphism, sub,
       quo, DirectSum, ModuleIsomorphism, free_module, vector_space,
       module_homomorphism, direct_sum, module_isomorphism, basis

export Generic

###############################################################################
#
#   Polynomial Ring S, x = R["x"] syntax
#
###############################################################################

getindex(R::NCRing, s::Union{String, Char}) = PolynomialRing(R, s)
getindex(R::NCRing, s::Union{String, Char}, ss::Union{String, Char}...) =
   PolynomialRing(R, [s, ss...])

# syntax x = R["x"]["y"]
getindex(R::Tuple{Union{Ring, NCRing}, Union{PolyElem, NCPolyElem}}, s::Union{String, Char}) = PolynomialRing(R[1], s)

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
#   Printing
#
################################################################################

include("PrettyPrinting.jl")

################################################################################
#
#   Deprecations
#
################################################################################

include("Deprecations.jl")

###############################################################################
#
#   Package handle creation
#
###############################################################################

const package_handle = [1]

function get_handle()
   package_handle[1] += 1
   return package_handle[1] - 1
end

###############################################################################
#
#   Auxilliary data accessors
#
###############################################################################

mutable struct AccessorNotSetError <: Exception
end

function create_accessors(T, S, handle)
   get = function(a, error::Bool = true)
      if handle > length(a.auxilliary_data) ||
         !isassigned(a.auxilliary_data, handle)
        if error
          throw(AccessorNotSetError())
        else
          return nothing
        end
      end
      return a.auxilliary_data[handle]
   end
   set = function(a, b)
      if handle > length(a.auxilliary_data)
         resize!(a.auxilliary_data, handle)
      end
      a.auxilliary_data[handle] = b
   end
   return get, set
end

###############################################################################
#
#   Promote rule helpers
#
###############################################################################

function sig_exists(T::Type{Tuple{U, V, W}}, sig_table::Array{X, 1}) where {U, V, W, X}
   for s in sig_table
      if s === T
         return true
      end
   end
   return false
end

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

needs_parentheses(x) = false

function isnegative end

# TODO remove these two once removed from dependent packages (Hecke)
@noinline function displayed_with_minus_in_front(x)
  Base.depwarn("This function is deprecated",
               :displayed_with_minus_in_front)
  return string(x)[1] == "-"
end

@noinline function show_minus_one(x)
  Base.depwarn("This function is deprecated",
               :show_minus_one)
  return true
end

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

   cmd = "using Test; using AbstractAlgebra; include(\"$test_file\"); $test_function_name();"
   @info("spawning ", `$julia_exe -e \"$cmd\"`)
   run(`$julia_exe -e $cmd`)
end

end # module
