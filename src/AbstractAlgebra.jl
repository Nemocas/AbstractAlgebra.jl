@doc raw"""
AbstractAlgebra is a pure Julia package for computational abstract algebra.

Repository: <https://github.com/Nemocas/AbstractAlgebra.jl>

Documentation: <https://nemocas.github.io/AbstractAlgebra.jl>
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
                  :InfiniteDimensionError, # remove in next breaking release, see #2135
                  :identity_map, # see #5188 in Oscar.jl
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
include("PrettyOrdering.jl")

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

pretty_lt(x::Number, y::Number) = isless(x, y)
pretty_eq(x::Number, y::Number) = (x == y)

include("julia/JuliaTypes.jl")

# Unions of AbstractAlgebra abstract types and Julia types
const JuliaRingElement = Union{Integer, Rational, AbstractFloat}
const JuliaFieldElement = Union{Rational, AbstractFloat}
const JuliaExactRingElement = Union{Integer, Rational}

const RingElement = Union{RingElem, JuliaRingElement}
const NCRingElement = Union{NCRingElem, JuliaRingElement}

const FieldElement = Union{FieldElem, JuliaFieldElement}

include("ConcreteTypes.jl")

# Type aliases
const UniversalPolyRing{T} = UniversalRing{<:MPolyRingElem, T}
const UniversalPolyRingElem{T} = UniversalRingElem{<:MPolyRingElem, T}
const UnivPoly{T} = UniversalPolyRingElem{T}

###############################################################################
#
#   Fundamental interface for AbstractAlgebra
#
###############################################################################

# prevent the default Julia methods for iszero and isone from being invoked on
# one of our parent objects: the default implementations are
#   iszero(x) = x == zero(x)    and    isone(x) = x == one(x)
# which for a parent object will always return false, which usually is not
# what the user hoped for. So let's rather throw an error and leave it to
# subtypes to implement better methods (e.g. `iszero(::Ring)` or potentially
# isone(::Group) come to mind).
iszero(S::Set) = throw(NotImplementedError(:iszero, S))
isone(S::Set) = throw(NotImplementedError(:isone, S))

include("fundamental_interface.jl")
include("misc/VarNames.jl")
include("Infinity.jl")

################################################################################
#
#   Printing
#
################################################################################

include("PrettyPrinting.jl")

using .PrettyPrinting

import .PrettyPrinting: expressify


################################################################################
#
#   Conformance tests (function stubs for TestExt)
#
################################################################################

include("ConformanceTests.jl")

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

function check_base_ring(a, b, throw::Bool = true)
   flag = base_ring(a) === base_ring(b)
   flag || !throw || error("base rings do not match")
   return flag
end

include("algorithms/LaurentPoly.jl")
include("algorithms/FinField.jl")
include("algorithms/GenericFunctions.jl")

include("Groups.jl")
include("Rings.jl")
include("NCRings.jl")
include("Fields.jl")
include("Factor.jl")

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
include("UniversalLaurentPoly.jl")
include("UniversalRing.jl")
include("UnivPoly.jl")
include("FreeAssociativeAlgebra.jl")
include("LaurentMPoly.jl")
include("MatrixNormalForms.jl")

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

@doc raw"""
    is_pairwise(f, a::AbstractVector)

Return whether `f(a[i], a[j])` is true for all `i != j`. It is assumed that `f`
is symmetric.

# Examples

```jldoctest
julia> is_pairwise(!is_associated, [2, -3, 5])
true

julia> is_pairwise(!is_associated, [2, -3, -2])
false
```
"""
function is_pairwise(fun, arr::AbstractVector)
  for i in 1:length(arr)
    for j in 1:(i - 1)
      if !fun(arr[i], arr[j])
        return false
      end
    end
  end
  return true
end


###############################################################################
#
#   methods to inquire about known properties of objects
#
###############################################################################

include("KnownProperties.jl")

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


# Generic functions to be defined after all rings
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
   rand_file = joinpath(pkgdir, "test/rand.jl")

   cmd = """
         using Test
         using AbstractAlgebra
         include("$rand_file")
         include("$test_file")
         """
   @info("spawning ", `$julia_exe --project=$(Base.active_project()) -e $cmd`)
   run(`$julia_exe -e $cmd`)
end

function doctestsetup()
  return :(using AbstractAlgebra; AbstractAlgebra.set_current_module(@__MODULE__))
end

###############################################################################
#
#   Utils
#
###############################################################################
include("utils.jl")

using PrecompileTools: @setup_workload, @compile_workload    # this is a small dependency

@setup_workload begin
  # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
  # precompile file and potentially make loading faster.
  @compile_workload begin
    # Julia.Integers.conformance_tests
    ConformanceTests.exercise_Ring_interface_recursive(ZZ)
    ConformanceTests.exercise_EuclideanRing_interface(ZZ)

    # Generic.Poly.conformance
    R, x = polynomial_ring(ZZ, "x")
    ConformanceTests.exercise_Poly_interface(R)
    R, x = polynomial_ring(QQ, "x")
    ConformanceTests.exercise_Poly_interface(R)
    R, x = polynomial_ring(GF(5), "x")
    ConformanceTests.exercise_Poly_interface(R)

    # EuclideanRingResidueRingElem.conformance_tests
    R, = residue_ring(ZZ, 1)
    ConformanceTests.exercise_Ring_interface(R)   # is_gen fails on polys

    R, = residue_ring(ZZ, -4)
    ConformanceTests.exercise_Ring_interface_recursive(R)

    R, = residue_ring(ZZ, 16453889)
    ConformanceTests.exercise_Ring_interface_recursive(R)

    S, x = polynomial_ring(R, "x")
    R, = residue_ring(S, x^3 + 3x + 1)
    ConformanceTests.exercise_Ring_interface_recursive(R)

    S, x = polynomial_ring(QQ, "x")
    R, = residue_ring(S, x^2 + 1)
    ConformanceTests.exercise_Ring_interface_recursive(R)
  end
end


end # module
