# Credits: record on demand which software and papers a computation uses.
module Credits

# JULIA_1_10: Julia 1.10 has no Base.ScopedValues. The package ScopedValues.jl
# is a backport for Julia 1.10; on Julia ≥ 1.11 it only re-exports
# Base.ScopedValues. To drop Julia 1.10:
#   1. Remove all code marked JULIA_1_10 (`grep -rn JULIA_1_10`).
#   2. Here, use `using Base.ScopedValues: ScopedValue, @with`.
#   3. Remove ScopedValues from Project.toml.
using ScopedValues: ScopedValue, @with

export Credit, @credit, @collect_credits

"""
    Credit(name, bibtex)

Software or a paper to cite. Define each credit once, as a `const`.
"""
struct Credit
  name::String
  bibtex::String
end

"""
Result of `@collect_credits`: credit => (location of `@credit` => number of hits).
"""
const CreditCollection = Dict{Credit, Dict{String, Int}}

# The collection of the innermost running `@collect_credits`, or `nothing`.
# A task inherits this value from the task that starts it.
const COLLECTION = ScopedValue{Union{Nothing, CreditCollection}}(nothing)

# Guards all access to collections: tasks on several threads can record into
# the same collection.
const LOCK = ReentrantLock()

function record_credit(c::Credit, site::String, n::Int = 1)
  col = COLLECTION[]
  col === nothing && return
  @lock LOCK begin
    sites = get!(Dict{String, Int}, col, c)
    sites[site] = get(sites, site, 0) + n
  end
  return
end

"""
    @credit C

Record a hit of `C::Credit` if a `@collect_credits` runs. Otherwise do nothing.
"""
macro credit(c)
  site = "$(__module__) at $(__source__.file):$(__source__.line)"
  return :(record_credit($(esc(c)), $site))
end

"""
    cred = @collect_credits expr

Evaluate `expr` and return the credits it hits. Assignments in `expr` stay
visible: `@collect_credits r = f()` sets `r`. An enclosing `@collect_credits`
also receives the credits.
"""
macro collect_credits(ex)
  # `@with`, not `with(() -> ex, ...)`: no closure, so assignments stay visible
  return quote
    col = CreditCollection()
    @with COLLECTION => col $(esc(ex))
    finish(col)
  end
end

# Pass the credits in `col` on to the enclosing collection, if any, and return
# a copy: tasks started in `expr` can still run and record into `col`.
function finish(col::CreditCollection)
  @lock LOCK begin
    for (c, sites) in col, (site, n) in sites
      record_credit(c, site, n)
    end
    CreditCollection(c => copy(sites) for (c, sites) in col)
  end
end

"""
    bibtex(cred::CreditCollection)

The BibTeX entries of all credits in `cred`.
"""
bibtex(cred::CreditCollection) = join(sort!([c.bibtex for c in keys(cred)]), "\n\n")

end # module Credits
