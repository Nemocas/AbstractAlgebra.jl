module AbstractAlgebra

using Markdown

using InteractiveUtils

# A list of all symbols external packages should not import from AbstractAlgebra
import_exclude = [:import_exclude, :QQ, :ZZ,
                  :RealField, :NumberField,
                  :AbstractAlgebra,
                  :exp, :sqrt, :div, :divrem, :numerator, :denominator,
                  :promote_rule,
                  :Set, :Module, :Ring, :Group, :Field]

import LinearAlgebra: det, norm, nullspace, rank, transpose!, hessenberg

import LinearAlgebra: lu, lu!, tr

export nullspace

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
const Base_import_exclude = [:exp, :sqrt, :div, :divrem, :numerator, :denominator]

################################################################################
#
#  Functions that we do not import from Base
#
################################################################################

function exp(a::T) where T
   return Base.exp(a)
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

function numerator(a::T, canonicalise::Bool=true) where T
  return Base.numerator(a, canonicalise)
end

function denominator(a::T, canonicalise::Bool=true) where T
  return Base.denominator(a, canonicalise)
end

import Base: Array, abs, acos, acosh, adjoint, asin, asinh, atan, atanh, bin,
             ceil, checkbounds, conj, convert, cmp, cos, cosh, cospi, cot,
             coth, dec, deepcopy, deepcopy_internal, expm1, exponent, fill,
             floor, gcd, gcdx, getindex, hash, hcat, hex, hypot, intersect,
             inv, invmod, isequal, isfinite, isless, isone, isqrt, isreal,
             iszero, lcm, ldexp, length, log, log1p, mod, ndigits, oct, one,
             parent, parse, precision, rand, Rational, rem, reverse, setindex!,
             show, sincos, similar, sign, sin, sinh, sinpi, size, string, tan,
             tanh, trailing_zeros, transpose, truncate, typed_hvcat,
             typed_hcat, typed_vcat, vcat, xor, zero, zeros,
             +, -, *, ==, ^, &, |, <<, >>, ~, <=, >=, <, >, //, /, !=

using Random: Random, AbstractRNG

export elem_type, parent_type

export SetElem, GroupElem, NCRingElem, RingElem, ModuleElem, FieldElem, RingElement,
       FieldElement, Map, AccessorNotSetError

export SetMap, FunctionalMap, IdentityMap

export NCPolyElem, PolyElem, SeriesElem, AbsSeriesElem, RelSeriesElem, ResElem, FracElem,
       MatElem, MatAlgElem, FinFieldElem, MPolyElem, NumFieldElem, SimpleNumFieldElem

export PolyRing, SeriesRing, ResRing, FracField, MatSpace, MatAlgebra,
       FinField, MPolyRing, NumField, SimpleNumField

export ZZ, QQ, zz, qq, RealField, RDF

export PermutationGroup

export create_accessors, get_handle, package_handle, zeros,
       Array, sig_exists

export error_dim_negative, ErrorConstrDimMismatch

export crt

include("AbstractTypes.jl")

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
#   Generic submodule
#
###############################################################################

include("Generic.jl")

# Do not import numerator and denominator as we have our own
import .Generic: add!, addeq!, addmul!, add_column, add_column!, add_row, add_row!, base_ring, cached,
                 canonical_unit, can_solve_left_reduced_triu, change_base_ring,
                 character, characteristic, charpoly, charpoly_danilevsky!,
                 charpoly_danilevsky_ff!, charpoly_hessenberg!, chebyshev_t,
                 chebyshev_u, _check_dim, check_composable, check_parent,
                 codomain, coeff, coeffs, ncols,
                 combine_like_terms!, compose, content, cycles,
                 data, deflate, deflation, degree, degrees,
                 dense_matrix_type, derivative, det, det_clow,
                 det_df, det_fflu, det_popov, dim, disable_cache!,
                 discriminant, displayed_with_minus_in_front,
                 divexact, divexact_left, divexact_right, divides,
                 domain, downscale,
                 elem_type, enable_cache!, evaluate, exp_gcd,
                 exponent, exponent_vector, exponent_vectors,
                 extended_weak_popov, extended_weak_popov_with_transform,
                 finish, fflu!,
                 fflu, find_pivot_popov, fit!, gcd, gen,
                 gens, get_field, gcdinv, gcdx,
                 gram, has_left_neighbor, has_bottom_neighbor, hash,
                 hessenberg!, hnf, hnf_cohen, hnf_cohen_with_transform,
                 hnf_kb, hnf_kb_with_transform,
                 hnf_minors, hnf_minors_with_transform,
                 hnf_with_transform, hnf_via_popov,
                 hnf_via_popov_with_transform,
                 hooklength, identity_map, identity_matrix, image,
                 image_map, image_fn, inflate, integral, interpolate, inv,
                 inv!, invariant_factors,
                 inverse_fn, inverse_image_fn,
                 inverse_mat, reverse_rows, reverse_rows!,
                 reverse_cols, reverse_cols!,
                 invmod, involves_at_most_one_variable,
                 iscompatible, isconstant, isdegree,
                 isdomain_type, isexact_type, isgen, ishessenberg,
                 ishnf, ishomogeneous, isisomorphic, ismonomial,
                 isone, isreverse, isrimhook,
                 isrref, issquare, issubmodule, isterm,
                 isunit, iszero_row, iszero_column,
                 kernel, kronecker_product,
                 laurent_ring, lc, lcm, lead, left_kernel, length,
                 leglength, lm, lt, main_variable,
                 main_variable_extract, main_variable_insert,
                 map1, map2, map_from_func, map_coeffs, map_entries, map_entries!,
                 map_with_preimage_from_func, map_with_retraction,
                 map_with_retraction_from_func,
                 map_with_section, map_with_section_from_func, mat, matrix,
                 matrix_repr, max_fields, max_precision, minors, minpoly, mod,
                 modulus, monomial, monomial!, monomials,
                 monomial_iszero, monomial_set!, monomial_to_newton!,
                 MPolyBuildCtx, mul!, mul_classical, mul_karatsuba, mul_ks,
                 mullow, mulmod, multiply_column, multiply_column!,
                 multiply_row, multiply_row!, needs_parentheses,
                 newton_to_monomial!, ngens, normalise, nrows, nvars, O, one,
                 order, ordering, parent_type, parity, partitionseq, Perm, perm,
                 permtype, @perm_str, polcoeff, pol_length, powmod,
                 pow_multinomial, popov, popov_with_transform, 
                 precision, preimage, preimage_map, primpart, pseudodivrem,
                 pseudo_inv, pseudorem, push_term!, rank, randmat_triu,
                 randmat_with_rank, rand_ordering, rank_profile_popov, remove,
                 renormalize!, rels, rescale!, resultant, resultant_ducos,
                 resultant_euclidean, resultant_subresultant,
                 resultant_sylvester, resx, retraction_map, reverse,
                 right_kernel, rref, rref!, section_map, setcoeff!,
                 set_exponent_vector!, set_field!, set_length!, set_limit!,
                 setpermstyle, set_prec!, set_val!, size, shift_left,
                 shift_right, show_minus_one, similarity!, snf, snf_kb,
                 snf_kb_with_transform, snf_with_transform, solve, solve_left,
                 solve_rational, solve_triu, sort_terms!, sub, subst, summands,
                 supermodule, swap_cols, swap_cols!, swap_rows, swap_rows!,
                 sylvester_matrix, symbols, term, terms, total_degree,
                 to_univariate, trail, truncate, typed_hcat, typed_hvcat,
                 upscale, valuation, var, var_index, vars, weak_popov,
                 weak_popov_with_transform, zero, zero!, zero_matrix,
                 @PolynomialRing, MatrixElem

# Do not export divrem, exp, sqrt, numerator and denominator as we define our own
export add!, addeq!, addmul!, addmul_delayed_reduction!, addmul!, add_column, add_column!, add_row, add_row!, base_ring, cached,
                 canonical_unit, can_solve_left_reduced_triu,
                 change_base_ring, character,
                 characteristic, charpoly, charpoly_danilevsky!,
                 charpoly_danilevsky_ff!, charpoly_hessenberg!, chebyshev_t,
                 chebyshev_u, _check_dim, check_composable, check_parent,
                 codomain, coeff, coeffs, ncols,
                 combine_like_terms!, compose, content, cycles,
                 data, deflate, deflation, degree, degrees,
                 dense_matrix_type, derivative, det, det_clow,
                 det_df, det_fflu, det_popov, dim, disable_cache!,
                 discriminant, displayed_with_minus_in_front,
                 divexact, divexact_left, divexact_right, divides,
                 domain, downscale,
                 elem_type, enable_cache!, evaluate, exp_gcd,
                 exponent, exponent_vector, exponent_vectors,
                 extended_weak_popov,
                 extended_weak_popov_with_transform, fflu!,
                 fflu, find_pivot_popov, finish, fit!, gcd, gen,
                 gens, get_field, gcdinv, gcdx,
                 gram, has_left_neighbor, has_bottom_neighbor, hash,
                 hessenberg!, hessenberg, hnf,
                 hnf_cohen, hnf_cohen_with_transform,
                 hnf_kb, hnf_kb_with_transform,
                 hnf_minors, hnf_minors_with_transform,
                 hnf_with_transform, hnf_via_popov,
                 hnf_via_popov_with_transform,
                 hooklength, identity_map, identity_matrix, image,
                 image_map, image_fn, inflate, integral, interpolate, inv,
                 inv!, invariant_factors,
                 inverse_fn, inverse_image_fn,
                 inverse_mat, invmod, reverse_rows,
                 reverse_rows!, reverse_cols, reverse_cols!,
                 iscompatible, isconstant, isdegree,
                 isdomain_type, isexact_type, isgen, ishessenberg,
                 ishnf, ishomogeneous, isisomorphic, ismonomial,
                 isnegative, isone, isreverse,
                 isrimhook, isrref, issquare, issubmodule,
                 isterm, isunit, iszero,
                 iszero_row, iszero_column, kernel,
                 kronecker_product, laurent_ring,
                 lc, lcm, lead, left_kernel, leglength,
                 length, lm, lu, lu!, lt,
                 main_variable, main_variable_extract, main_variable_insert,
                 map1, map2, map_from_func, map_coeffs, map_entries, map_entries!, map_with_preimage_from_func,
                 map_with_retraction, map_with_retraction_from_func,
                 map_with_section, map_with_section_from_func,
                 mat, matrix, matrix_repr, max_fields,
                 max_precision, minors, minpoly, mod,
                 modulus, monomial, monomial!, monomials,
                 monomial_iszero, monomial_set!, monomial_to_newton!,
                 MPolyBuildCtx, mul!, mul_classical,
                 mul_karatsuba, mul_ks, mul_red!, mullow, mulmod,
                 multiply_column, multiply_column!, multiply_row,
                 multiply_row!, needs_parentheses, newton_to_monomial!, ngens,
                 normalise, nrows, nvars, O, one, order, ordering, parent_type,
                 parity, partitionseq, Perm, perm, permtype, @perm_str, polcoeff,
                 pol_length, powmod, pow_multinomial, popov,
                 popov_with_transform, powers, ppio, precision, preimage,
                 preimage_map, primpart, pseudo_inv, pseudodivrem, pseudorem,
                 push_term!, rank, randmat_triu, randmat_with_rank,
                 rand_ordering, rank_profile_popov, reduce!, remove,
                 renormalize!, rels, resultant, resultant_ducos, rescale!,
                 resultant_euclidean, resultant_subresultant,
                 resultant_sylvester, resx, retraction_map, reverse,
                 right_kernel, rref, rref!, section_map, setcoeff!,
                 set_exponent_vector!, set_field!, set_length!, set_limit!,
                 setpermstyle, set_prec!, set_val!, shift_left, shift_right,
                 show_minus_one, similarity!, size, snf, snf_kb,
                 snf_kb_with_transform, snf_with_transform, solve, solve_left,
                 solve_rational, solve_triu, sort_terms!, sub, subst, summands,
                 supermodule, swap_rows, swap_rows!, swap_cols, swap_cols!,
                 sylvester_matrix, symbols, term, terms, to_univariate,
                 total_degree, tr, trail, truncate, typed_hcat, typed_hvcat,
                 upscale, valuation, var, var_index, vars, weak_popov,
                 weak_popov_with_transform, zero, zero!, zero_matrix,
                 @PolynomialRing, MatrixElem

################################################################################
#
#   Parent constructors
#
################################################################################

function PermGroup(n::T) where T
  Generic.PermGroup(n)
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

function PowerSeriesRing(R::Ring, prec::Int, s::AbstractString; cached=true, model=:capped_relative)
   Generic.PowerSeriesRing(R, prec, s; cached=cached, model=model)
end

function PowerSeriesRing(R::Ring, prec::Int, s::Char; cached=true, model=:capped_relative)
   PowerSeriesRing(R, prec, string(s); cached=cached, model=model)
end

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

function PolynomialRing(R::Ring, s::AbstractString; cached::Bool = true)
   Generic.PolynomialRing(R, s; cached=cached)
end

function PolynomialRing(R::Ring, s::Char; cached::Bool = true)
   PolynomialRing(R, string(s); cached=cached)
end

function PolynomialRing(R::NCRing, s::AbstractString; cached::Bool = true)
   Generic.PolynomialRing(R, s; cached=cached)
end

function PolynomialRing(R::NCRing, s::Char; cached::Bool = true)
   PolynomialRing(R, string(s); cached=cached)
end

function PolynomialRing(R::Ring, s::Array{String, 1}; cached::Bool = true, ordering::Symbol = :lex)
   Generic.PolynomialRing(R, s; cached=cached, ordering=ordering)
end

function PolynomialRing(R::Ring, s::Array{Char, 1}; cached::Bool = true, ordering::Symbol = :lex)
   PolynomialRing(R, string.(s); cached=cached, ordering=ordering)
end

function SparsePolynomialRing(R::Ring, s::String; cached::Bool = true)
   Generic.SparsePolynomialRing(R, s; cached=cached)
end

function SparsePolynomialRing(R::Ring, s::Char; cached::Bool = true)
   SparsePolynomialRing(R, string(s); cached=cached)
end

function MatrixSpace(R::Ring, r::Int, c::Int, cached::Bool = true)
   Generic.MatrixSpace(R, r, c, cached)
end

function MatrixAlgebra(R::Ring, n::Int, cached::Bool = true)
   Generic.MatrixAlgebra(R, n, cached)
end

function FractionField(R::Ring; cached=true)
   Generic.FractionField(R; cached=cached)
end

function ResidueRing(R::Ring, a::Union{RingElement, Integer}; cached::Bool = true)
   Generic.ResidueRing(R, a; cached=cached)
end

function ResidueField(R::Ring, a::Union{RingElement, Integer}; cached::Bool = true)
   Generic.ResidueField(R, a; cached=cached)
end

function NumberField(a::AbstractAlgebra.Generic.Poly{Rational{BigInt}}, s::AbstractString, t = "\$"; cached = true)
   Generic.NumberField(a, s, t; cached=cached)
end

function NumberField(a::AbstractAlgebra.Generic.Poly{Rational{BigInt}}, s::Char, t = "\$"; cached = true)
   NumberField(a, string(s), t; cached=cached)
end

@doc Markdown.doc"""
    FreeModule(R::NCRing, rank::Int; cached::Bool = true)
> Return the free module over the ring $R$ with the given rank.
"""
function FreeModule(R::NCRing, rank::Int; cached::Bool = true)
   Generic.FreeModule(R, rank; cached=cached)
end

function free_module(R::NCRing, rank::Int; cached::Bool = true)
   Generic.FreeModule(R, rank; cached=cached)
end

@doc Markdown.doc"""
    VectorSpace(R::Field, dim::Int; cached::Bool = true)
> Return the vector space over the field $R$ with the given dimension.
"""
function VectorSpace(R::Field, dim::Int; cached::Bool = true)
   Generic.FreeModule(R, dim)
end

function vector_space(R::Field, dim::Int; cached::Bool = true)
   Generic.FreeModule(R, dim)
end

@doc Markdown.doc"""
    sub(m::Module{T}, gens::Vector{<:ModuleElem{T}}) where T <: RingElement
> Return the submodule `S` of the module `m` generated by the given generators,
> given as elements of `m`, and a map which is the canonical injection from `S`
> to `m`.
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
> Return the submodule `S` of the module `m` generated by the union of the given
> submodules of $m$, and a map which is the canonical injection from `S` to `m`.
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
> Return the quotient `Q` of the module `m` by the submodule `N` of `m`, and a
> map which is a lift of elements of `Q` to `m`.
"""
function quo(m::Module{T}, N::Module{T}) where T <: RingElement
   Generic.quo(m, N)
end

@doc Markdown.doc"""
    DirectSum(m::Vector{<:Module{T}}) where T <: RingElement
> Return a tuple $M, f, g$ consisting of $M$ the direct sum of the modules `m`
> (supplied as a vector of modules), a vector $f$ of the canonical injections
> of the $m[i]$ into $M$ and a vector $g$ of the canonical projections from
> $M$ onto the $m[i]$.
"""
function DirectSum(m::Vector{<:Module{T}}) where T <: RingElement
   Generic.DirectSum(m)
end

function direct_sum(m::Vector{<:Module{T}}) where T <: RingElement
   Generic.DirectSum(m)
end

@doc Markdown.doc"""
    DirectSum(m::Module{T}...) where T <: RingElement
> Return a tuple $M, f, g$ consisting of $M$ the direct sum of the given
> modules, a vector $f$ of the canonical injections of the $m[i]$ into $M$
> and a vector $g$ of the canonical projections from $M$ onto the $m[i]$.
"""
function DirectSum(m::Module{T}...) where T <: RingElement
   Generic.DirectSum(m...)
end

function direct_sum(m::Module{T}...) where T <: RingElement
   Generic.DirectSum(m...)
end

function ModuleHomomorphism(M1::AbstractAlgebra.Module, M2::AbstractAlgebra.Module, m::MatElem)
   Generic.ModuleHomomorphism(M1, M2, m)
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

# add empty functions so that Singular, Nemo and Hecke can import and extend.
function crt end

export PowerSeriesRing, PolynomialRing, SparsePolynomialRing, MatrixSpace,
       MatrixAlgebra, FractionField, ResidueRing, Partition, PermGroup,
       YoungTableau, AllParts, SkewDiagram, AllPerms, Perm, LaurentSeriesRing,
       LaurentSeriesField, ResidueField, NumberField, PuiseuxSeriesRing,
       PuiseuxSeriesField, FreeModule, VectorSpace, ModuleHomomorphism, sub,
       quo, DirectSum, ModuleIsomorphism, free_module, vector_space,
       module_homomorphism, direct_sum, module_isomorphism

export Generic

###############################################################################
#
#   Polynomial Ring S, x = R["x"] syntax
#
###############################################################################

getindex(R::NCRing, s::Union{String, Char}) = PolynomialRing(R, s)

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
#   Set domain for PermutationGroup
#
###############################################################################

const PermutationGroup = PermGroup

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
#  For backwards compability
#
###############################################################################

displayed_with_minus_in_front(x) = isnegative(x)

isnegative(x) = displayed_with_minus_in_front(x)

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
