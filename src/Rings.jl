###############################################################################
#
#   Rings.jl : Generic commutative rings
#
###############################################################################

function isequal(a::RingElem, b::RingElem)
   return parent(a) == parent(b) && a == b
end

# Implement `isapprox` for ring elements via equality by default. On the one
# hand, we need isapprox methods to be able to conformance test series rings.
# On the other hand this is essentially the only sensible thing to do in
# positive characteristic so we might as well do it in a generic method.
function Base.isapprox(x::NCRingElem, y::NCRingElem;
                       atol::Real=0, rtol::Real=0,
                       nans::Bool=false, norm::Function=abs)
  if is_exact_type(typeof(x)) && is_exact_type(typeof(y))
    @req is_zero(atol) "non-zero atol not supported"
    @req is_zero(rtol) "non-zero rtol not supported"
    return x == y
  end
  throw(NotImplementedError(:isapprox, x, y))
end

"""
    divexact(x, y; check::Bool=true)

Return an exact quotient of `x` by `y`, i.e. an element
`z` such that `x == yz`; when `x` and `y` do not belong to the same ring,
they are first coerced into a common ring.
By default if no exact division is possible, an exception is raised. If
`check=false` this check may be omitted for performance reasons and the
behaviour of the function undefined if the division is not exact.
"""
function divexact end

function divexact(x::RingElem, y::RingElem; check::Bool=true)
  xx, yy = promote(x, y)
  # - if divexact is not implemented, we need to break the recursion
  #   we assume that promotion returns identical operands if no proper
  #   promotion can be performed
  # - add type checks to make it constant-fold away in most cases
  if typeof(xx) === typeof(x) && typeof(y) === typeof(yy) && (xx, yy) === (x, y)
    throw(NotImplementedError(:divexact, x, y))
  end
  return divexact(xx, yy; check=check)
end

divexact(x::RingElem, y::RingElement; check::Bool=true) = divexact(x, parent(x)(y); check=check)

divexact(x::RingElement, y::RingElem; check::Bool=true) = divexact(parent(y)(x), y; check=check)

divexact_left(x::T, y::T; check::Bool=true) where T <: RingElement = divexact(x, y; check=check)

divexact_right(x::T, y::T; check::Bool=true) where T <: RingElement = divexact(x, y; check=check)

Base.:/(x::ModuleElem, y::RingElement) = divexact(x, y; check=true)
Base.:/(x::RingElem, y::RingElem) = divexact(x, y; check=true)
Base.:/(x::RingElem, y::JuliaRingElement) = divexact(x, y; check=true)
Base.:/(x::JuliaRingElement, y::RingElem) = divexact(x, y; check=true)

Base.inv(x::RingElem) = divexact(one(parent(x)), x)

@doc raw"""
    is_divisible_by(x::T, y::T) where T <: RingElem

Check if `x` is divisible by `y`, i.e. if $x = zy$ for some $z$.
"""
function is_divisible_by(x::T, y::T) where T <: RingElem
   if iszero(y)
      return iszero(x)
   end
   return divides(x, y)[1]
end

@doc raw"""
    is_associated(x::T, y::T) where T <: RingElem

Check if `x` and `y` are associated, i.e. if `x` is a unit times `y`.
"""
function is_associated(x::T, y::T) where T <: RingElem
   return is_divisible_by(x, y) && is_divisible_by(y, x)
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(x::PolyRingElem{T}, y::Integer) where T <: RingElem
   return evaluate(x, base_ring(x)(y))
end

function evaluate(x::MPolyRingElem{T}, y::Integer) where T <: RingElem
   return evaluate(x, base_ring(x)(y))
end

###############################################################################
#
#   Type traits
#
###############################################################################

# Type can only represent elements of an exact ring
# true unless explicitly specified
#
# implementors should only implement this trait for RingElem subtypes, but for
# convenience we support calling this also on Ring subtypes as well as Ring
# and RingElem instances
is_exact_type(R::Type{T}) where T <: RingElem = true

is_exact_type(x) = is_exact_type(typeof(x))
is_exact_type(x::Type{<:Ring}) = is_exact_type(elem_type(x))
is_exact_type(T::DataType) = throw(MethodError(is_exact_type, (T,)))

# Type can only represent elements of domains, i.e. without zero divisors
# false unless explicitly specified
#
# implementors should only implement this trait for RingElem subtypes, but for
# convenience we support calling this also on Ring subtypes as well as Ring
# and RingElem instances
is_domain_type(R::Type{T}) where T <: NCRingElem = false

is_domain_type(x) = is_domain_type(typeof(x))
is_domain_type(x::Type{<:NCRing}) = is_domain_type(elem_type(x))
is_domain_type(T::DataType) = throw(MethodError(is_domain_type, (T,)))

###############################################################################
#
#   Exponential function for generic rings
#
###############################################################################

function Base.exp(a::RingElem)
   a != 0 && error("Exponential of nonzero element")
   return one(parent(a))
end

###############################################################################
#
#   Coprime bases
#
###############################################################################

# Bernstein, "Factoring into coprimes in essentially linear time"
# ppio(a,b) = (c,n) where v_p(c) = v_p(a) if v_p(b) != 0, 0 otherwise
# c*n = a or c = gcd(a, b^infty), n = div(a, c).
# This is used in various Euclidean domains for Chinese remaindering.

function ppio(a::E, b::E) where E <: RingElem
   if iszero(a)
     return one(parent(a)), a
   end
   c = gcd(a, b)
   n = div(a, c)
   g = gcd(c, n)
   while !isone(g)
      c *= g
      n = div(n, g)
      g = gcd(c, n)
   end
   return c, n
end

################################################################################
#
#   Squares
#
################################################################################

function Base.sqrt(a::FieldElem; check::Bool=true)
  R = parent(a)
  R, t = polynomial_ring(R, :t; cached = false)
  f = factor(t^2 - a)
  for (p, e) in f
    if !check || degree(p) == 1
      return -divexact(coeff(p, 0), coeff(p, 1); check=check)
    end
  end
  error("Element $a does not have a square root")
end

# assumes the existence of is_square and sqrt for input
function is_square_with_sqrt(a::RingElem)
  if is_square(a)
     return true, sqrt(a)
  else
     return false, parent(a)()
  end
end

###############################################################################
#
#   Ring properties
#
###############################################################################

@doc raw"""
    is_trivial(R::NCRing)

Test whether the ring $R$ is trivial. A ring is trivial if it consists
of a single element, or equivalently if its characteristic is 1. Such
rings are also called zero rings.
"""
is_trivial(R::NCRing) = !is_domain_type(elem_type(R)) && iszero(one(R))
is_known(::typeof(is_trivial), R::NCRing) = is_domain_type(elem_type(R))

@doc raw"""
    is_zero(R::NCRing)

Test whether the ring $R$ is trivial. However, the recommended method
for doing this is [`is_trivial`](@ref).

```jldoctest
julia> R, (x,) = polynomial_ring(QQ, [:x]);

julia> S=quo(R,1)[1]
Residue ring of R modulo 1

julia> is_trivial(S)
true

julia> is_zero(S)
true
```
"""
is_zero(R::NCRing) = is_trivial(R)

@doc raw"""
    is_perfect(F::Field)

Test whether the field $F$ is perfect, that is, whether the characteristic is zero or
else whether every element of $F$ admits a $p$-th root, where $p > 0$ is the characteristic of $F$.
"""
is_perfect(F::Field) = characteristic(F) == 0 || throw(NotImplementedError(:is_perfect, F))
is_known(::typeof(is_perfect), F::Field) = is_known(characteristic, F) && characteristic(F) == 0

is_perfect(F::FinField) = true
is_known(::typeof(is_perfect), F::FinField) = true

is_finite(F::FinField) = true
is_known(::typeof(is_finite), F::FinField) = true

function is_finite(R::NCRing)
  c = characteristic(R)
  c == 0 && return false
  c == 1 && return true
  throw(NotImplementedError(:is_finite, R))
end
is_known(::typeof(is_finite), R::NCRing) = is_known(characteristic, R) && characteristic(R) <= 1

@doc raw"""
    is_noetherian(R::Ring)

Check if the ring $R$ is Noetherian.

# Examples
```jldoctest
julia> R, x = polynomial_ring(ZZ, [:x]);

julia> is_noetherian(R)
true
```
"""
function is_noetherian(R::Ring)
  throw(NotImplementedError(:is_noetherian, R))
end

is_noetherian(::Field) = true
is_noetherian(::Integers) = true

is_noetherian(R::Union{PolyRing, MPolyRing, LaurentPolyRing, LaurentMPolyRing}) = is_noetherian(coefficient_ring(R))
is_noetherian(R::Union{MSeriesRing, SeriesRing}) = is_noetherian(base_ring(R))
is_noetherian(R::ResidueRing) = is_noetherian(base_ring(R)) || throw(NotImplementedError(:is_noetherian, R))

@doc raw"""
    krull_dim(R::Ring)

Return the Krull dimension of the ring $R$. The method is currently only supported for certain commutative Noetherian Rings.

# Examples
```jldoctest
julia> R, x = polynomial_ring(ZZ, [:x]);

julia> krull_dim(R)
2
```
"""
function krull_dim(R::Ring)
  @req is_noetherian(R) "Krull dimension is only supported for Noetherian rings."
  throw(NotImplementedError(:krull_dim, R))
end

krull_dim(::Field) = 0
krull_dim(::Integers) = 1

krull_dim(R::Union{MPolyRing, PolyRing, LaurentMPolyRing, LaurentPolyRing}) = krull_dim(coefficient_ring(R)) + nvars(R)
krull_dim(R::Union{SeriesRing, MSeriesRing}) = krull_dim(base_ring(R)) + nvars(R)


########################################################################
# locality checks
########################################################################
@doc raw"""
    is_local(R::Ring)

Return whether a given ring `R` is a local ring.
"""
function is_local(R::Ring)
  throw(NotImplementedError(:is_local, R))
end

# In general we can not assume it to be known whether a given ring is local
is_known(::typeof(is_local), R::Ring) = false

is_local(::Field) = true

is_local(R::MPolyRing{<:FieldElem}) = is_zero(ngens(R))
is_known(::typeof(is_local), R::MPolyRing{<:FieldElem}) = true

