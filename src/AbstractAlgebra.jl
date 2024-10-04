@doc raw"""
 AbstractAlgebra is a pure Julia package for computational abstract algebra.

 For more information see https://github.com/Nemocas/AbstractAlgebra.jl
"""
module AbstractAlgebra

include("imports.jl")

# A list of all symbols external packages should not import from AbstractAlgebra
const import_exclude = [:import_exclude, :QQ, :ZZ,
                  :RealField, :GF,
                  :AbstractAlgebra,
                  :inv, :log, :exp, :sqrt, :div, :divrem,
                  :numerator, :denominator,
                  :promote_rule,
                  :Set, :Module, :Group,
                 ]


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


include("exports.jl")
include("AliasMacro.jl")
include("Aliases.jl") # needs to be included after AliasMacro.jl
include("Assertions.jl")

include("Attributes.jl")
include("PrintHelper.jl")

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
function force_coerce(a, b, ::Val{throw_error} = Val(true)) where {throw_error}
  if throw_error
    error("coercion not possible")
  end
  return nothing
end

#to allow +(a::T, b::T) where a, b have different parents, but
# a common over structure
# designed(?) to be minimally invasive in AA and Nemo, but filled with
# content in Hecke/Oscar
function force_op(op::Function, ::Val{throw_error}, a...) where {throw_error}
  if throw_error
    error("no common overstructure for the arguments found")
  end
  return false
end

function force_op(op::Function, a...)
  return force_op(op, Val(true), a...)
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

include("julia/JuliaTypes.jl")

# Unions of AbstactAlgebra abstract types and Julia types
const RingElement   = Union{RingElem,   Integer, Rational, AbstractFloat}
const NCRingElement = Union{NCRingElem, Integer, Rational, AbstractFloat}

const FieldElement = Union{FieldElem, Rational, AbstractFloat}

###############################################################################
#
#   Fundamental interface for AbstractAlgebra
#
###############################################################################

include("fundamental_interface.jl")
include("misc/VarNames.jl")

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

function check_parent(a, b, throw::Bool = true)
   flag = parent(a) === parent(b)
   flag || !throw || error("parents do not match")
   return flag
end

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
include("FreeAssociativeAlgebra.jl")
include("LaurentMPoly.jl")
include("MatrixNormalForms.jl")


include("Groups.jl")
include("Rings.jl")
include("NCRings.jl")
include("Fields.jl")
include("Factor.jl")

# More functionality for Julia types
include("julia/Integer.jl")
include("julia/Rational.jl")
include("julia/Float.jl")
include("julia/GF.jl")
include("julia/Matrix.jl")

###############################################################################
#
#   Generic submodule
#
###############################################################################

include("Generic.jl")

using .Generic

###############################################################################
#
#   Linear solving submodule
#
###############################################################################

include("Solve.jl")

using .Solve

################################################################################
#
#   Number fields (some stuff moved from Nemo)
#
################################################################################

include("NumFields.jl")

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
  return [canonical_projection(D, i) for i=1:_number_of_direct_product_factors(D)]
end

###############################################################################
#
#   misc
#
###############################################################################

include("misc/ProductIterator.jl")
include("misc/Evaluate.jl")

###############################################################################
#
#   Syntax S[i] for all parents S as a shortcut for gen(S, i)
#
################################################################################

getindex(S::Set, i::Int) = gen(S, i)

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


# Generic functions to be defined after all rings
include("polysubst.jl")

include("broadcasting.jl")

################################################################################
#
#   Deprecations
#
################################################################################

include("Deprecations.jl")


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
include("algorithms/coprime_base.jl")

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

###############################################################################
#
#   Utils
#
###############################################################################
include("utils.jl")

end # module
