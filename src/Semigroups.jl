###############################################################################
#
#   Semigroup and monoid interface
#
###############################################################################

Base.parent(x::SemigroupElem) = throw(NotImplementedError(:parent, x))

Base.:*(x::T, y::T) where {T<:SemigroupElem} = throw(NotImplementedError(:*, x, y))

Base.one(M::Monoid) = throw(NotImplementedError(:one, M))

Base.one(x::MonoidElem) = one(parent(x))

Base.isone(x::MonoidElem) = x == one(x)

function Base.:^(x::SemigroupElem, n::Integer)
  n > 0 || throw(DomainError(n, "semigroup powers require a positive exponent"))
  n == 1 && return deepcopy(x)
  return internal_power(x, n)
end

function Base.:^(x::MonoidElem, n::Integer)
  n >= 0 || throw(DomainError(n, "monoid powers require a nonnegative exponent"))
  n == 0 && return one(parent(x))
  n == 1 && return deepcopy(x)
  return internal_power(x, n)
end

# Respect the power interface even if a representation happens to supply an
# ambient identity or a generalized inverse. Neither establishes a monoid or
# group power operation for the specified parent.
Base.literal_pow(::typeof(^), x::SemigroupElem, ::Val{n}) where {n} = x^n

Base.broadcastable(x::SemigroupElem) = Ref(x)
