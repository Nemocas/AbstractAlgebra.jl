module TestExt

using AbstractAlgebra
isdefined(Base, :get_extension) ? (using Test) : (using ..Test)


include("Groups-conformance-tests.jl")


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
