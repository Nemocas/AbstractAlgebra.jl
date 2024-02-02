###############################################################################
#
#   Rings.jl : Generic commutative rings
#
###############################################################################

function isequal(a::RingElem, b::RingElem)
   return parent(a) == parent(b) && a == b
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

divexact(x::RingElem, y::RingElem; check::Bool=true) = divexact(promote(x, y)...; check=check)

divexact(x::RingElem, y::RingElement; check::Bool=true) = divexact(x, parent(x)(y); check=check)

divexact(x::RingElement, y::RingElem; check::Bool=true) = divexact(parent(y)(x), y; check=check)

divexact_left(x::T, y::T; check::Bool=true) where T <: RingElement = divexact(x, y; check=check)

divexact_right(x::T, y::T; check::Bool=true) where T <: RingElement = divexact(x, y; check=check)

Base.:/(x::ModuleElem, y::RingElement) = divexact(x, y; check=true)
Base.:/(x::RingElem, y::RingElem) = divexact(x, y; check=true)
Base.:/(x::RingElem, y::Union{Integer, Rational, AbstractFloat}) = divexact(x, y; check=true)
Base.:/(x::Union{Integer, Rational, AbstractFloat}, y::RingElem) = divexact(x, y; check=true)

Base.inv(x::RingElem) = divexact(one(parent(x)), x)

function is_divisible_by(x::T, y::T) where T <: RingElem
   if iszero(y)
      return iszero(x)
   end
   return divides(x, y)[1]
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
is_exact_type(R::Type{T}) where T <: RingElem = true

# Type can only represent elements of domains
# false unless explicitly specified
is_domain_type(R::Type{T}) where T <: RingElem = false

###############################################################################
#
#   Exponential function for generic rings
#
###############################################################################

function Base.exp(a::RingElem)
   a != 0 && error("Exponential of nonzero element")
   return one(parent(a))
end

################################################################################
#
#   Transpose for ring elements
#
################################################################################

transpose(x::T) where {T <: RingElem} = deepcopy(x)

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

@doc raw"""
    sqrt(a::FieldElem)

Return the square root of the element `a`. By default the function will
throw an exception if the input is not square. If `check=false` this test is
omitted.
"""
function Base.sqrt(a::FieldElem; check::Bool=true)
  R = parent(a)
  R, t = polynomial_ring(R, "t", cached = false)
  f = factor(t^2 - a)
  for (p, e) in f
    if !check || degree(p) == 1
      return -divexact(coeff(p, 0), coeff(p, 1); check=check)
    end
  end
  error("Element $a does not have a square root")
end

# assumes the existence of sqrt without check argument for input
function Base.sqrt(a::RingElem; check::Bool=true)
  s = sqrt(a)
  if check
    s != a^2 && error("Element $a does not have a square root")
  end
  return s
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

is_perfect(F::Field) = characteristic(F) == 0 || F isa FinField ||
                                                 throw(NotImplementedError(:is_perfect, F))

is_finite(F::FinField) = true

is_finite(F::Field) = characteristic(F) != 0 && throw(NotImplementedError(:is_finite, F))

characteristic(F::NumField) = 0
