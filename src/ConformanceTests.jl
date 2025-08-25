module ConformanceTests

using ..AbstractAlgebra

# helper

function equality(a, b)
  if is_exact_type(typeof(a)) && is_exact_type(typeof(b))
     return a == b
  else
     return isapprox(a, b)
  end
end

function equality(a::Tuple, b::Tuple)
  return all(equality.(a, b))
end


function equality_up_to_units(a, b)
  iszero(a) && return iszero(b)
  iszero(b) && return iszero(a)
  return divides(a, b)[1] && divides(b, a)[1]
end

const default_adhoc_partner_rings = [
    AbstractAlgebra.Integers{BigInt}(),
    AbstractAlgebra.Integers{Int}(),
    AbstractAlgebra.Integers{UInt}(),
    AbstractAlgebra.Integers{UInt8}(),
  ]

adhoc_partner_rings(R::NCRing) = default_adhoc_partner_rings

"""
    generate_element(::AbstractAlgebra.NCRing)

Generates some random, representative element of the given set for use
in the conformance tests.

It is supposed to be implemented in the src file of the respective type.
"""
function generate_element end


###############################################################################
#
#  `implements` trait
#
###############################################################################

# Calling `_implements(T, f)` checks whether a "sensible" method for the unary
# function `f` is implemented for inputs of type `T`. The argument order is
# meant to be similar to e.g. `isa`, and thus indicates `T implements f`.
#
# For example, `_implements(MyRingElem, is_unit)` should return true if
# invoking `is_unit` on elements of type `MyRingElem` is supported.
#
# The generic fallback uses `hasmethod`. However, this may return `true` in
# cases where it shouldn't, as we often provide generic methods for that rely
# on other methods being implemented -- either for the same type, or for types
# derived from it. For example the `is_nilpotent(::PolyElem{T})` method needs
# `is_nilpotent(::T)` in order to work.
#
# To reflect this, additional `_implements` methods need to be provided.
# We currently do this for at least the following functions:
# - factor
# - is_irreducible
# - is_nilpotent
# - is_squarefree
# - is_unit
# - is_zero_divisor
#
_implements(::Type{T}, f::Any) where {T} = hasmethod(f, Tuple{T})

# Alternatively, the first argument can be a concrete object. By default we
# then redispatch to the type based version. But one may also choose to
# implement custom methods for this: certain operations will only work for
# *some* instances. E.g. for `Z/nZ` it may happen that for `n` a prime we can
# perform a certain operation, but not if `n` is composite.
#
# In that case the recommendation is that `_implements` invoked on the type
# returns `false`, but invoked on a concrete instance of a type, it may use
# specifics of the instance to also return `true` if appropriate.
function _implements(x::T, f::Any) where {T}
  @assert !(x isa Type) # paranoia
  return _implements(T, f)
end

# helper for `_implements` which checks if `f` has a method explicitly for
# a concrete type `T` (i.e. not a generic method that can be specialized to `T`
# but really one that is implement for `T` and `T` only).
function _implements_directly(::Type{T}, f::Any) where {T}
  isconcretetype(T) || return false  # TODO: drop this?
  meth = methods(f, Tuple{T})
  # TODO: deal with type parameters: if `T` is `FreeAssociativeAlgebraElem{ZZRingElem}`
  # and `f` has a method for `FreeAssociativeAlgebraElem` then we should still consider
  # this a match.
  return any(m -> m.sig == Tuple{typeof(f), T}, meth)
end


###############################################################################
#
# The following function stubs' actual implementations are in the folder `ext/TestExt/`.
#
###############################################################################

function test_iterate end

# Groups-conformance-tests.jl
function test_Group_interface end
function test_GroupElem_interface end

# Mutating-ops.jl
function test_mutating_op_like_zero end
function test_mutating_op_like_neg end
function test_mutating_op_like_add end
function test_mutating_op_like_addmul end

# Rings-conformance-tests.jl
function test_NCRing_interface end
function test_Ring_interface end
function test_Field_interface end
function test_EuclideanRing_interface end
function test_Poly_interface end
function test_MPoly_interface end
function test_MatSpace_interface end
function test_MatAlgebra_interface end
function test_Ring_interface_recursive end
function test_Field_interface_recursive end


end # module
