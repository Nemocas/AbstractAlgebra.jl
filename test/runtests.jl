using AbstractAlgebra

using Test


# disable until we encounter GC problems again
if false VERSION >= v"1.8.0"
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

include("Aqua.jl")
include("rand.jl")
include("AbstractAlgebra-test.jl")
