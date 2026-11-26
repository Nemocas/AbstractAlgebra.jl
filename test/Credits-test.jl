module CreditsTest
  # Import under a different name: the macros must not need `AbstractAlgebra` by name.
  import AbstractAlgebra as AA

  const SOFTA = AA.Credit("softA", "@misc{softA, title = {Software A}}")
  const PAPERB = AA.Credit("paperB", "@article{paperB, title = {Algorithm B}}")
  a(x) = (AA.@credit SOFTA; 2x)
  b(x) = (AA.@credit PAPERB; a(x) + 1)
end

# Test that @collect_credits records the hits of @credit, and nothing otherwise.
@testset "Credits" begin
  using .CreditsTest: SOFTA, PAPERB, a, b
  COLLECTION = AbstractAlgebra.Credits.COLLECTION

  @test b(1) == 3
  @test COLLECTION[] === nothing

  cred = @collect_credits (b(1), a(1))
  @test keys(cred) == Set([SOFTA, PAPERB])
  @test only(values(cred[SOFTA])) == 2
  @test only(values(cred[PAPERB])) == 1
  @test occursin("CreditsTest at ", only(keys(cred[SOFTA])))
  @test occursin("Credits-test.jl:", only(keys(cred[SOFTA])))
  @test AbstractAlgebra.Credits.bibtex(cred) ==
        "@article{paperB, title = {Algorithm B}}\n\n@misc{softA, title = {Software A}}"

  @test isempty(@collect_credits 1 + 1)

  cred = @collect_credits r = b(2)
  @test r == 5
  @test keys(cred) == Set([SOFTA, PAPERB])

  function g()
    cred = @collect_credits r = a(3)
    return r, length(cred)
  end
  @test g() == (6, 1)

  local inner
  outer = @collect_credits begin
    a(1)
    inner = @collect_credits b(1)
  end
  @test keys(inner) == Set([SOFTA, PAPERB])
  @test only(values(outer[SOFTA])) == 2

  cred = @collect_credits fetch.([Threads.@spawn a(i) for i in 1:1000])
  @test only(values(cred[SOFTA])) == 1000

  local tasks
  cred = @collect_credits tasks = [Threads.@spawn (sleep(0.1); a(1)) for _ in 1:4]
  foreach(wait, tasks)
  @test isempty(cred)

  @test_throws ErrorException @collect_credits (a(1); error("boom"))
  @test COLLECTION[] === nothing
end

# Test hits inside `with_logger`. JULIA_1_10: on Julia 1.10, the ScopedValues.jl
# backport keeps the scope in the current logger, and `with_logger` hides it.
@testset "Credits and loggers (Julia $VERSION)" begin
  using .CreditsTest: SOFTA, a
  using Base.CoreLogging: current_logger, with_logger, NullLogger
  outside = current_logger()
  cred = @collect_credits begin
    inside = current_logger()
    with_logger(() -> a(1), NullLogger())
    @test_logs (:info, "hi") (@info "hi"; a(1))   # @test_logs uses with_logger
  end
  if VERSION < v"1.11"   # JULIA_1_10: remove this branch, keep the else branch
    @test inside !== outside
    @test isempty(cred)
    ScopedValues = parentmodule(AbstractAlgebra.Credits.ScopedValue)
    cred = @collect_credits ScopedValues.with_logger(() -> a(1), NullLogger())
    @test only(values(cred[SOFTA])) == 1
  else
    @test inside === outside
    @test only(values(cred[SOFTA])) == 2
  end
end
