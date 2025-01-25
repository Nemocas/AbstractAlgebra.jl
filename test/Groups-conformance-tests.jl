# This file is just for backwards compatibility and should be removed in the future.

# Assumptions (as that was already required with the previous setup):
# - `AbstractAlgebra` is already loaded
# - `Test` is available in the environment

# load TestExt
using Test

import .ConformanceTests: test_Group_interface
import .ConformanceTests: test_GroupElem_interface
