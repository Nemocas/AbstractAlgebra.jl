# This file is just for backwards compatibility and should be removed in the future.

# Assumptions (as that was already required with the previous setup):
# - `AbstractAlgebra` is already loaded
# - `Test` is available in the environment

# load TestExt
using Test

import .ConformanceTests: test_elem as test_elem
import .ConformanceTests: adhoc_partner_rings as adhoc_partner_rings
import .ConformanceTests: equality as equality

test_mutating_op_like_zero = ConformanceTests.test_mutating_op_like_zero
test_mutating_op_like_neg = ConformanceTests.test_mutating_op_like_neg
test_mutating_op_like_add = ConformanceTests.test_mutating_op_like_add
test_mutating_op_like_addmul = ConformanceTests.test_mutating_op_like_addmul

test_NCRing_interface = ConformanceTests.test_NCRing_interface
test_Ring_interface = ConformanceTests.test_Ring_interface
test_Field_interface = ConformanceTests.test_Field_interface
test_EuclideanRing_interface = ConformanceTests.test_EuclideanRing_interface
test_Poly_interface = ConformanceTests.test_Poly_interface
test_MPoly_interface = ConformanceTests.test_MPoly_interface
test_MatSpace_interface = ConformanceTests.test_MatSpace_interface
test_MatAlgebra_interface = ConformanceTests.test_MatAlgebra_interface
test_Ring_interface_recursive = ConformanceTests.test_Ring_interface_recursive
test_Field_interface_recursive = ConformanceTests.test_Field_interface_recursive
