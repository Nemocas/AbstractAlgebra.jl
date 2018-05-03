VERSION >= v"0.4.0-dev+6521" && __precompile__()

module AbstractAlgebra

import Base: Array, abs, acos, acosh, asin, asinh, atan, atan2, atanh, base,
             bin, ceil, checkbounds, conj, convert, cmp, contains, cos, cosh,
             cospi, cot, coth, dec, deepcopy, deepcopy_internal, deserialize,
             det, div, divrem, expm1, eye, floor, gamma, gcd, gcdx, getindex,
             hash, hcat, hex, hypot, intersect, inv, invmod, isequal,
             isfinite, isless, isqrt, isreal, iszero, lcm, ldexp, length,
             lgamma, log, log1p, lufact, lufact!, mod, ndigits, nextpow2, norm,
             nullspace, numerator, oct, one, parent, parse, precision,
             prevpow2, rand, rank, Rational, rem, reverse, serialize,
             setindex!, show, similar, sign, sin, sinh, sinpi, size, string,
             tan, tanh, trace, trailing_zeros, transpose, transpose!, truncate,
             typed_hvcat, typed_hcat, var, vcat, xor, zero, zeros, +, -, *, ==, ^,
             &, |, <<, >>, ~, <=, >=, <, >, //, /, !=

export elem_type, parent_type

export SetElem, GroupElem, RingElem, ModuleElem, FieldElem, RingElement,
       FieldElement, Map, AccessorNotSetError

export SetMap, FunctionalMap, IdentityMap

export PolyElem, SeriesElem, AbsSeriesElem, RelSeriesElem, ResElem, FracElem,
       MatElem, FinFieldElem, MPolyElem

export PolyRing, SeriesRing, ResRing, FracField, MatSpace,
       FinField, MPolyRing

export ZZ, QQ, zz, qq, RealField, RDF

export PermutationGroup

export create_accessors, get_handle, package_handle, zeros,
       Array, sig_exists

export error_dim_negative, ErrorConstrDimMismatch

export FiniteField, crt

if VERSION >= v"0.6.0-dev.2024" # julia started exporting iszero (again?)
   import Base: iszero
end

if VERSION < v"0.6-"
   import Base: isprime, factor, parity, sub, call
end

if VERSION >= v"0.7.0-DEV.264" # julia started exporting sincos
   import Base: sincos
end

if VERSION >= v"0.7.0-DEV.1144"
    import Base: isone
end

include("AbstractTypes.jl")

###############################################################################
#
#   Set up environment / load libraries
#
###############################################################################

function __init__()
   println("")
   println("Welcome to AbstractAlgebra version 0.0.6")
   println("")
   println("AbstractAlgebra comes with absolutely no warranty whatsoever")
   println("")
end

###############################################################################
#
#  Version information
#
################################################################################

function versioninfo()
  print("AbstractAlgebra version 0.0.4\n")
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

import .Generic: add!, addeq!, addmul!, base_ring, cached, canonical_unit,
                 character, characteristic, charpoly, charpoly_danilevsky!,
                 charpoly_danilevsky_ff!, charpoly_hessenberg!, chebyshev_t,
                 chebyshev_u, _check_dim, check_parent, codomain, coeff, cols,
                 compose, content, cycles, data, deflate, degree,
                 denominator, derivative, det, det_clow,
                 det_df, det_fflu, det_popov, dim, disable_cache!, discriminant,
                 divexact, divides, divrem, domain, downscale,
                 elem_type, elements, enable_cache!, evaluate, exp_gcd,
                 extended_weak_popov, extended_weak_popov_with_trafo, fflu!,
                 fflu, find_pivot_popov, fit!, gcd, gen,
                 gens, get_field, gcdinv, gcdx,
                 gram, has_left_neighbor, has_bottom_neighbor, hash,
                 hessenberg!, hessenberg, hnf, hnf_cohen, hnf_cohen_with_trafo,
                 hnf_kb, hnf_kb_with_trafo, hnf_minors, hnf_minors_with_trafo,
                 hnf_with_trafo, hnf_via_popov, hnf_via_popov_with_trafo,
                 hooklength, identity_map, identity_matrix,
                 image_map, image_fn, inflate, inskewdiag, integral,
                 interpolate, inv, inv!, inverse_fn, invmod, isconstant, isdegree,
                 isdomain_type, isexact_type, isgen, ishessenberg, ismonomial,
                 isnegative, isone, isreverse, isrimhook,
                 isrref, issquare, isterm,
                 isunit, laurent_ring, lcm, lead, length,
                 leglength, main_variable,
                 main_variable_extract, main_variable_insert,
                 map1, map2, map_from_func,
                 map_with_preimage_from_func, map_with_retraction,
                 map_with_retraction_from_func,
                 map_with_section, map_with_section_from_func, matrix,
                 matrix_repr, max_degrees, max_precision, minpoly, mod,
                 modulus, monomial_iszero, monomial_set!, monomial_to_newton!,
                 mul!, mul_classical, mul_karatsuba, mul_ks, mullow, mulmod,
                 needs_parentheses, newton_to_monomial!,
                 normalise, nvars, numerator,
                 O, one, order, ordering, parent_type, parity, partitionseq,
                 perm, polcoeff, pol_length, powmod, pow_multinomial, popov,
                 popov_with_trafo, powers, precision,
                 preimage_map, primpart, pseudodivrem, pseudorem, randmat_triu,
                 randmat_with_rank, rand_ordering, rank_profile_popov, remove,
                 renormalize!, rescale!, resultant,
                 resultant_ducos, resultant_euclidean,
                 resultant_subresultant, resultant_sylvester, resx, 
                 retraction_map, reverse,
                 rows, rref, rref!, section_map, setcoeff!,
                 set_field!, set_length!, set_limit!, setpermstyle,
                 set_prec!, set_val!, size, shift_left,
                 shift_right, show_minus_one,
                 similarity!, snf, snf_kb, snf_kb_with_trafo, snf_with_trafo,
                 solve, solve_rational, solve_triu, sub, subst, swap_rows,
                 swap_rows!, trail, truncate, typed_hcat, typed_hvcat, upscale,
                 valuation, var, vars, weak_popov, weak_popov_with_trafo, zero,
                 zero!, zero_matrix, kronecker_product

export add!, addeq!, addmul!, base_ring, cached, canonical_unit, character,
                 characteristic, charpoly, charpoly_danilevsky!,
                 charpoly_danilevsky_ff!, charpoly_hessenberg!, chebyshev_t,
                 chebyshev_u, _check_dim, check_parent, codomain, coeff, cols,
                 compose, content, cycles, data, deflate, degree,
                 denominator, derivative, det, det_clow,
                 det_df, det_fflu, det_popov, dim, disable_cache!, discriminant,
                 divexact, divides, divrem, domain, downscale,
                 elem_type, elements, enable_cache!, evaluate, exp_gcd,
                 extended_weak_popov, extended_weak_popov_with_trafo, fflu!,
                 fflu, find_pivot_popov, fit!, gcd, gen, 
                 gens, get_field, gcdinv, gcdx,
                 gram, has_left_neighbor, has_bottom_neighbor, hash,
                 hessenberg!, hessenberg, hnf, hnf_cohen, hnf_cohen_with_trafo,
                 hnf_kb, hnf_kb_with_trafo, hnf_minors, hnf_minors_with_trafo,
                 hnf_with_trafo, hnf_via_popov, hnf_via_popov_with_trafo,
                 hooklength, identity_map, identity_matrix,
                 image_map, image_fn, inflate, inskewdiag, integral,
                 interpolate, inv, inv!, inverse_fn, invmod, isconstant, isdegree,
                 isdomain_type, isexact_type, isgen, ishessenberg,
                 ismonomial, isnegative, isone, isreverse,
                 isrimhook, isrref, issquare,
                 isterm, isunit, iszero, laurent_ring,
                 lcm, lead, leglength, length,
                 main_variable, main_variable_extract, main_variable_insert,
                 map1, map2, map_from_func, map_with_preimage_from_func,
                 map_with_retraction, map_with_retraction_from_func,
                 map_with_section, map_with_section_from_func,
                 matrix, matrix_repr, max_degrees, max_precision, minpoly, mod,
                 modulus, monomial_iszero, monomial_set!, monomial_to_newton!,
                 mul!, mul_classical, mul_karatsuba, mul_ks, mullow, mulmod,
                 needs_parentheses, newton_to_monomial!, normalise, nvars,
                 numerator, O, one, order, ordering, parent_type, parity,
                 partitionseq, perm, polcoeff, pol_length,
                 powmod, pow_multinomial,
                 popov, popov_with_trafo, powers, ppio, 
                 precision, preimage_map, primpart,
                 pseudodivrem, pseudorem, randmat_triu, randmat_with_rank,
                 rand_ordering, rank_profile_popov, remove,
                 renormalize!, resultant,
                 resultant_ducos, rescale!, resultant_euclidean,
                 resultant_subresultant,
                 resultant_sylvester, resx, retraction_map,
                 reverse, rows, rref, rref!,
                 section_map, setcoeff!, set_field!, set_length!, 
                 set_limit!, setpermstyle, set_prec!, set_val!,
                 shift_left, shift_right, show_minus_one, similarity!, size,
                 snf, snf_kb, snf_kb_with_trafo, snf_with_trafo, solve,
                 solve_rational, solve_triu, sub, subst, swap_rows, swap_rows!,
                 trail, truncate, typed_hcat, typed_hvcat,
                 upscale, valuation, var,
                 vars, weak_popov, weak_popov_with_trafo, zero, zero!,
                 zero_matrix, kronecker_product

function exp(a::T) where T
   return Base.exp(a)
end

function sqrt(a::T) where T
  return Base.sqrt(a)
end

function PermGroup(n::T) where T
  Generic.PermGroup(n)
end

function AllPerms(n::T) where T
  Generic.AllPerms(n)
end

function perm(n::T) where T
  Generic.perm(n)
end

function perm(a::Array{T, 1}) where T
  Generic.perm(a)
end

function Partition(part::Vector{T}, check::Bool=true) where T
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

function LaurentSeriesRing(R::Ring, prec::Int, s::AbstractString; cached=true)
   Generic.LaurentSeriesRing(R, prec, s; cached=cached)
end

function LaurentSeriesRing(R::Field, prec::Int, s::AbstractString; cached=true)
   Generic.LaurentSeriesField(R, prec, s; cached=cached)
end

function LaurentSeriesField(R::Field, prec::Int, s::AbstractString; cached=true)
   Generic.LaurentSeriesField(R, prec, s; cached=cached)
end

function PuiseuxSeriesRing(R::Ring, prec::Int, s::AbstractString; cached=true)
   Generic.PuiseuxSeriesRing(R, prec, s; cached=cached)
end

function PuiseuxSeriesRing(R::Field, prec::Int, s::AbstractString; cached=true)
   Generic.PuiseuxSeriesField(R, prec, s; cached=cached)
end

function PuiseuxSeriesField(R::Field, prec::Int, s::AbstractString; cached=true)
   Generic.PuiseuxSeriesField(R, prec, s; cached=cached)
end

function PolynomialRing(R::Ring, s::AbstractString; cached::Bool = true)
   Generic.PolynomialRing(R, s; cached=cached)
end

function PolynomialRing(R::Ring, s::Array{String, 1}; cached::Bool = true, ordering::Symbol = :lex)
   Generic.PolynomialRing(R, s; cached=cached, ordering=ordering)
end

function SparsePolynomialRing(R::Ring, s::String; cached::Bool = true)
   Generic.SparsePolynomialRing(R, s; cached=cached)
end

function MatrixSpace(R::Ring, r::Int, c::Int, cached::Bool = true)
   Generic.MatrixSpace(R, r, c, cached)
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

#add empty functions so that Singular, Nemo and Hecke can import and extend.
function crt(A...)
  return AbstractAlgebra.crt(A...)
end

export PowerSeriesRing, PolynomialRing, SparsePolynomialRing, MatrixSpace,
       FractionField, ResidueRing, Partition, PermGroup, YoungTableau,
       AllParts, SkewDiagram, AllPerms, perm, LaurentSeriesRing,
       LaurentSeriesField, ResidueField, NumberField, PuiseuxSeriesRing,
       PuiseuxSeriesField

export Generic

###############################################################################
#
#   Polynomial Ring S, x = R["x"] syntax
#
###############################################################################

getindex(R::Ring, s::String) = PolynomialRing(R, s)

getindex(R::Tuple{Ring, T}, s::String) where {T} = PolynomialRing(R[1], s)

###############################################################################
#
#   Matrix M = R[...] syntax
#
################################################################################

function typed_hvcat(R::Ring, dims, d...)
   T = elem_type(R)
   r = length(dims)
   c = dims[1]
   A = Array{T}(r, c)
   for i = 1:r
      dims[i] != c && throw(ArgumentError("row $i has mismatched number of columns (expected $c, got $(dims[i]))"))
      for j = 1:c
         A[i, j] = R(d[(i - 1)*c + j])
      end
   end
   S = matrix(R, A)
   return S
end

function typed_hcat(R::Ring, d...)
   T = elem_type(R)
   r = length(d)
   A = Array{T}(1, r)
   for i = 1:r
      A[1, i] = R(d[i])
   end
   S = matrix(R, A)
   return S
end

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

###########################################################
#
#   Package handle creation
#
###########################################################

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
   get = function(a)
      if handle > length(a.auxilliary_data) ||
         !isassigned(a.auxilliary_data, handle)
        throw(AccessorNotSetError())
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

if VERSION >= v"0.5.0-dev+3171"

function sig_exists(T::Type{Tuple{U, V, W}}, sig_table::Array{X, 1}) where {U, V, W, X}
   for s in sig_table
      if s === T
         return true
      end
   end
   return false
end

else

function sig_exists(T::Type{Tuple{U, V, W}}, sig_table::Array{X, 1}) where {U, V, W, X}
   return false
end

end # if VERSION

###############################################################################
#
#   Array creation functions
#
###############################################################################

Array(R::Ring, r::Int...) = Array{elem_type(R)}(r)

function zeros(R::Ring, r::Int...)
   T = elem_type(R)
   A = Array{T}(r)
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

PermutationGroup = PermGroup

###############################################################################
#
#   Set domain for ZZ, QQ
#
###############################################################################

ZZ = JuliaZZ
QQ = JuliaQQ

###############################################################################
#
#   Set domain for RealField
#
###############################################################################

RealField = JuliaRealField

###############################################################################
#
#   Test code
#
###############################################################################

function test_module(x, y)
   julia_exe = Base.julia_cmd()
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

   cmd = "using Base.Test; using AbstractAlgebra; include(\"$test_file\"); $test_function_name();"
   info("spawning ", `$julia_exe -e \"$cmd\"`)
   run(`$julia_exe -e $cmd`)
end

end # module
