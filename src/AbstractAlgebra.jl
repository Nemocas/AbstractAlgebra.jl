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
include("MPoly.jl")

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
                 exponent, exponent_vector,
                 finish, fit!, gcd, gcdx,
                 has_left_neighbor, has_bottom_neighbor, hash,
                 hooklength, identity_map,
                 image_map, image_fn,
                 inverse_fn, inverse_image_fn,
                 inverse_mat, reverse_rows, reverse_rows!,
                 inv!, invmod,
                 iscompatible, isdegree,
                 ishomogeneous, isisomorphic,
                 isone, isreverse, isrimhook, issubmodule,
                 isunit,
                 laurent_ring, lcm, leading_coefficient, leading_monomial,
		 leading_term, length,
                 leglength, main_variable,
                 main_variable_extract, main_variable_insert,
                 map1, map2, map_from_func,
                 map_with_preimage_from_func, map_with_retraction,
                 map_with_retraction_from_func,
                 map_with_section, map_with_section_from_func, mat,
                 matrix_repr, max_fields, mod,
                 monomial, monomial!, monomials,
                 monomial_iszero, monomial_set!,
                 MPolyBuildCtx, mullow_karatsuba,
                 ngens, norm, normalise,
		 num_coeff, one,
                 order, ordering, parity, partitionseq, Perm, perm,
                 permtype, @perm_str, polcoeff, poly, poly_ring,
                 precision, preimage, preimage_map,
		 prime, push_term!,
                 rand_ordering, reduce!,
                 rels, rel_series, rel_series_type,
		 rescale!, retraction_map, reverse,
                 right_kernel, rowlength, section_map, setcoeff!,
                 set_exponent_vector!, set_limit!,
                 setpermstyle, size,
                 sort_terms!, summands,
                 supermodule, term, terms, total_degree,
                 to_univariate, trailing_coefficient,
		 truncate, upscale, var_index,
                 zero,
       # Moved from Hecke into Misc
		 Loc, Localization, LocElem,
		 polynomial_to_power_sums,
		 power_sums_to_polynomial, roots, sturm_sequence

# Do not export inv, div, divrem, exp, log, sqrt, numerator and denominator as we define our own
export abs_series, abs_series_type,
                 addmul_delayed_reduction!, addmul!,
                 base_field, base_ring, basis,
                 canonical_unit, can_solve_left_reduced_triu,
                 change_base_ring, change_coefficient_ring, character,
                 chebyshev_t,
                 chebyshev_u, check_composable, check_parent,
                 collength, coefficient_ring, combine_like_terms!, cycles,
                 defining_polynomial, degrees,
                 dense_matrix_type, dense_poly_type, det,
                 discriminant,
                 elem_type,
                 exponent, exponent_vector,
                 finish, fit!, gcd, gen,
                 gens, gcdinv, gcdx,
                 has_left_neighbor, has_bottom_neighbor, hash,
                 interpolate,
                 inv!, inverse_image_fn,
                 inverse_mat, invmod,
                 iscompatible, isdegree,
                 isdomain_type, isexact_type, isgen,
                 ishomogeneous,
                 isisomorphic, ismonomial, ismonomial_recursive,
                 isnegative, isone, isreverse,
                 issubmodule, issymmetric,
                 isterm_recursive, isunit, iszero,
                 lcm, leading_coefficient, leading_monomial, leading_term,
		 length,
                 main_variable, main_variable_extract, main_variable_insert,
                 mat, matrix_repr, max_fields, mod,
                 monomial, monomial!, monomials,
                 monomial_iszero, monomial_set!, monomial_to_newton!,
                 MPolyBuildCtx,
                 mul_ks, mul_red!, mullow_karatsuba, mulmod,
                 needs_parentheses, newton_to_monomial!, ngens,
                 normalise, nullspace, num_coeff,
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
		 setcoeff!, set_exponent_vector!,
                 setpermstyle,
                 size, sort_terms!, subst, summands, supermodule,
                 sylvester_matrix, term, terms, to_univariate,
                 total_degree, trailing_coefficient, truncate,
                 var_index, zero,
                 MatrixElem, PolynomialElem,
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

function SkewDiagram(lambda::Generic.Partition, mu::Generic.Partition)
  Generic.SkewDiagram(lambda, mu)
end

function YoungTableau(part::Generic.Partition, tab::Array{Int, 2})
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

function FunctionField(p::Generic.Poly{Generic.Rat{T}}, s::Symbol; cached::Bool=true) where T <: FieldElement
   return Generic.FunctionField(p, s; cached=cached)
end

function FunctionField(p::Generic.Poly{Generic.Rat{T}}, s::AbstractString; cached::Bool=true) where T <: FieldElement
   return Generic.FunctionField(p, Symbol(s); cached=cached)
end

function FunctionField(p::Generic.Poly{Generic.Rat{T}}, s::Char; cached::Bool=true) where T <: FieldElement
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
