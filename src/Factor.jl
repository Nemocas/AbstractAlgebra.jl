#################################################################################
#
#   Factor.jl : Factorization
#
#################################################################################

export Fac, factors, unit

################################################################################
#
#   Type
#
################################################################################

type Fac{T <: RingElem}
   unit::T
   fac::Dict{T, Int}
end

doc"""
    unit(a::Fac{T}) -> T

> Returns the unit of the factorization.
"""
unit(a::Fac) = a.unit

#primes(a::Fac) = collect(keys(a.fac))

################################################################################
#
#   Syntax sugar
#
################################################################################

doc"""
    in(a::T, b::Fac{T})

> Test whether `a` is a factor of `b`.
"""
function Base.in{T <: RingElem}(a::T, b::Fac{T})
   a in keys(b.fac)
end

doc"""
    getindex(a::Fac{T}, b::T) -> Int

> If `b` is a factor of `a`, the corresponding exponent is returned. Otherwise
> an error is thrown.
"""
function getindex{T}(a::Fac{T}, b::T)
  if haskey(a.fac, b)
    return a.fac[b]
  else
    error("$b is not a factor of $a")
  end
end

################################################################################
#
#   String I/O
#
################################################################################

function Base.show(io::IO, a::Fac)
  print(io, "$(a.unit)")
  for (p, e) in a.fac
    bracket = needs_parentheses(p)
    print(io, " * ")
    bracket ? print(io, "($p)") : print(io, "$p")
    a[p] == 1 ? nothing : print(io, "^$e")
  end
end

################################################################################
#
#   Make Factor objects iterable
#
################################################################################

Base.start(a::Fac) = Base.start((a.fac))

Base.done(a::Fac, state) = Base.done((a.fac), state)

Base.next(a::Fac, state) = Base.next((a.fac), state)

doc"""
    length(a::Fac) -> Int

> Returns the number of factors of `a`, not including the unit.
"""
Base.length(a::Fac) = Base.length(a.fac)

################################################################################
#
#   Special functions for Fac{fmpz}
#
################################################################################

Base.in(b::Integer, a::Fac{fmpz}) = Base.in(fmpz(b), a)

Base.getindex(a::Fac{fmpz}, b::Integer) = getindex(a, fmpz(b))

