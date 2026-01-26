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
function test_MatRing_interface end
function test_Ring_interface_recursive end
function test_Field_interface_recursive end

# deprecated during 0.48.*
@deprecate test_MatAlgebra_interface test_MatRing_interface

end # module
