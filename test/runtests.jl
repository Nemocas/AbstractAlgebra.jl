using ParallelTestRunner
using AbstractAlgebra

# artificially low cutoffs for testing purposes
function AbstractAlgebra.hgcd_prefers_basecase(a::T, b::T) where T <: PolyRingElem
   return degree(b) < 5
end

function AbstractAlgebra.mat22_mul_prefers_classical(
   a11::T, a12::T, a21::T, a22::T,
   b11::T, b12::T, b21::T, b22::T
) where T <: PolyRingElem

   return degree(a11) + degree(a22) < 4 || degree(b11) + degree(b22) < 4
end

# disable until we encounter GC problems again
if false # VERSION >= v"1.8.0"
  GC.enable_logging(true)

  # print gc settings
  jlmax = @ccall jl_gc_get_max_memory()::UInt64
  totalmem = @ccall uv_get_total_memory()::UInt64
  constrmem = @ccall uv_get_constrained_memory()::UInt64
  println("Memory:   max: ", Base.format_bytes(jlmax))
  println("        total: ", Base.format_bytes(totalmem))
  println("       constr: ", Base.format_bytes(constrmem))

#= FIXME/TODO: in the future we may wish to experiment with limiting the GC heap here
  if VERSION >= v"1.10.0-"
    # adjust heap size hint
    memenv = parse(Int, get(ENV, "OSCARCI_MAX_MEM_GB", "5")) * 2^30
    println("Setting heap size hint to ", Base.format_bytes(memenv))
    @ccall jl_gc_set_max_memory(memenv::UInt64)::Cvoid
  end
=#
end

# Start with autodiscovered tests
testsuite = find_tests(pwd())

# Parse arguments
args = parse_args(ARGS)

if filter_tests!(testsuite, args)
  # Remove Banners test packages
  delete!(testsuite, "utils/Banners/ModA/src/ModA")
  delete!(testsuite, "utils/Banners/ModA/src/ModA")
  delete!(testsuite, "utils/Banners/ModA/src/ModA")
end

runtests(AbstractAlgebra, args; testsuite)
