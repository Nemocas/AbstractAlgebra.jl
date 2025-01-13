module TestExt

using AbstractAlgebra
isdefined(Base, :get_extension) ? (using Test) : (using ..Test)

using .ConformanceTests:
  equality,
  equality_up_to_units,
  adhoc_partner_rings,
  test_elem

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
