# This file is just for backwards compatibility and should be removed in the future.

# Assumptions (as that was already required with the previous setup):
# - `AbstractAlgebra` is already loaded
# - `Test` is available in the environment

# load TestExt
using Test

import .ConformanceTests: generate_element as test_elem
import .ConformanceTests: adhoc_partner_rings
import .ConformanceTests: equality

import .ConformanceTests: test_mutating_op_like_zero
import .ConformanceTests: test_mutating_op_like_neg
import .ConformanceTests: test_mutating_op_like_add
import .ConformanceTests: test_mutating_op_like_addmul

import .ConformanceTests: test_NCRing_interface
import .ConformanceTests: test_Ring_interface
import .ConformanceTests: test_Field_interface
import .ConformanceTests: test_EuclideanRing_interface
import .ConformanceTests: test_Poly_interface
import .ConformanceTests: test_MPoly_interface
import .ConformanceTests: test_MatSpace_interface
import .ConformanceTests: test_MatAlgebra_interface
import .ConformanceTests: test_Ring_interface_recursive
import .ConformanceTests: test_Field_interface_recursive
