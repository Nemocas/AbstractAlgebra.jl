# This file is just for backwards compatibility and should be removed in the future.

# Assumptions (as that was already required with the previous setup):
# - `AbstractAlgebra` is already loaded
# - `Test` is available in the environment

# load TestExt
using Test

test_Group_interface = ConformanceTests.test_Group_interface
test_GroupElem_interface = ConformanceTests.test_GroupElem_interface
