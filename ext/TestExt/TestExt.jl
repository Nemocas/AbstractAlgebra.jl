module TestExt

using AbstractAlgebra
using Test

using .ConformanceTests:
  equality,
  equality_up_to_units,
  adhoc_partner_rings,
  generate_element

import .ConformanceTests: test_mutating_op_like_zero
import .ConformanceTests: test_mutating_op_like_neg
import .ConformanceTests: test_mutating_op_like_add
import .ConformanceTests: test_mutating_op_like_addmul

import .ConformanceTests: test_Group_interface
import .ConformanceTests: test_GroupElem_interface

import .ConformanceTests: test_NCRing_interface
import .ConformanceTests: test_Ring_interface
import .ConformanceTests: test_Field_interface
import .ConformanceTests: test_EuclideanRing_interface
import .ConformanceTests: test_Poly_interface
import .ConformanceTests: test_MPoly_interface
import .ConformanceTests: test_MatSpace_interface
import .ConformanceTests: test_MatRing_interface
import .ConformanceTests: test_Ring_interface_recursive
import .ConformanceTests: test_Field_interface_recursive

include("Groups-conformance-tests.jl")
include("Mutating-ops.jl")
include("Rings-conformance-tests.jl")


function ConformanceTests.test_iterate(F::FinField)
  elts = collect(Iterators.take(F, 20))
  @test elts isa Vector{elem_type(F)}
  @test allunique(elts)
  @test length(elts) == min(order(F), 20)
  if order(F) < 100
    elts = collect(F)
    @test elts isa Vector{elem_type(F)}
    @test allunique(elts)
    @test length(elts) == order(F)
  end
  return elts
end

end # module
