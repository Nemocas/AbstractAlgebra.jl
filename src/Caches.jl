################################################################################
#
#  Caches
#
################################################################################

# A registered cache is found through a hidden constant of this type next to
# it: constants survive precompilation, entries pushed into a registry owned by
# AbstractAlgebra while another package precompiles do not.
struct RegisteredCache
  cache
  mod::Core.Module
  name::Symbol
  source::LineNumberNode
end

@doc raw"""
    AbstractAlgebra.@register_cache const name = cache

Define the constant `name` and register `cache`, which must support `empty!`
and `length`, to be emptied by [`AbstractAlgebra.clear_cache`](@ref). Use it
at the top level of any module:

    AbstractAlgebra.@register_cache const _my_cache = Dict{Int, BigInt}()
"""
macro register_cache(ex)
  if !(Meta.isexpr(ex, :const, 1) && Meta.isexpr(ex.args[1], :(=), 2) && ex.args[1].args[1] isa Symbol)
    throw(ArgumentError("expected `@register_cache const name = cache`"))
  end
  name = ex.args[1].args[1]
  entry = Symbol("#", name, "#cache")
  return quote
    $(esc(ex))
    const $(esc(entry)) = $RegisteredCache($(esc(name)), $__module__, $(QuoteNode(name)), $(QuoteNode(__source__)))
    nothing
  end
end

function _registered_caches()
  caches = RegisteredCache[]
  seen = Base.Set{Core.Module}()
  for m in Base.loaded_modules_array()
    _registered_caches!(caches, seen, m)
  end
  return caches
end

function _registered_caches!(caches::Vector{RegisteredCache}, seen::Base.Set{Core.Module}, m::Core.Module)
  m in seen && return
  push!(seen, m)
  for s in names(m; all = true)
    isdefined(m, s) || continue
    Base.isdeprecated(m, s) && continue
    x = getglobal(m, s)
    if x isa RegisteredCache
      push!(caches, x)
    elseif x isa Core.Module && parentmodule(x) === m
      _registered_caches!(caches, seen, x)
    end
  end
end

@doc raw"""
    AbstractAlgebra.clear_cache()

Empty all caches registered with [`AbstractAlgebra.@register_cache`](@ref) in
loaded modules.
"""
function clear_cache()
  for c in _registered_caches()
    empty!(c.cache)
  end
  return nothing
end

@doc raw"""
    AbstractAlgebra.print_caches(io::IO = stdout)

Print name, length and definition site of each cache registered with
[`AbstractAlgebra.@register_cache`](@ref) in loaded modules.
"""
function print_caches(io::IO = stdout)
  for c in _registered_caches()
    println(io, c.mod, ".", c.name, ": length ", length(c.cache), ", ", c.source.file, ":", c.source.line)
  end
end
