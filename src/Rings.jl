###############################################################################
#
#   Rings.jl : Generic commutative rings
#
###############################################################################

function isequal(a::RingElem, b::RingElem)
   return parent(a) == parent(b) && a == b
end

################################################################################
#
#   Promotion system
#
# The promote_rule functions are not extending Base.promote_rule. The
# AbstractAlgebra promotion system is orthogonal to the built-in julia promotion
# system. The julia system assumes that whenever you have a method signature of
# the form Base.promote_rule(::Type{T}, ::Type{S}) = R, then there is also a
# corresponding Base.convert(::Type{R}, ::T) and similar for S. Since we
# cannot use the julia convert system (we need an instance of the type and not
# the type), we cannot use the julia promotion system.
#
# The AbstractAlgebra promotion system is used to define catch all functions for
# arithmetic between arbitrary ring elements.
#
################################################################################

promote_rule(::Type{T}, ::Type{T}) where T <: RingElement = T

function promote_rule_sym(::Type{T}, ::Type{S}) where {T, S}
   U = promote_rule(T, S)
   if U !== Union{}
      return U
   else
      UU = promote_rule(S, T)
      return UU
   end
end

@inline function try_promote(x::S, y::T) where {S <: RingElem, T <: RingElem}
   U = promote_rule_sym(S, T)
   if S === U
      return true, x, parent(x)(y)
   elseif T === U
      return true, parent(y)(x), y
   else
      return false, x, y
   end
end

function Base.promote(x::S, y::T) where {S <: RingElem, T <: RingElem}
  fl, u, v = try_promote(x, y)
  if fl
    return u, v
  else
    error("Cannot promote to common type")
  end
end

###############################################################################
#
#   Generic catchall functions
#
###############################################################################

+(x::RingElem, y::RingElem) = +(promote(x, y)...)

+(x::RingElem, y::RingElement) = x + parent(x)(y)

+(x::RingElement, y::RingElem) = parent(y)(x) + y

-(x::RingElem, y::RingElem) = -(promote(x, y)...)

-(x::RingElem, y::RingElement) = x - parent(x)(y)

-(x::RingElement, y::RingElem) = parent(y)(x) - y

*(x::RingElem, y::RingElem) = *(promote(x, y)...)

*(x::RingElem, y::RingElement) = x*parent(x)(y)

*(x::RingElement, y::RingElem) = parent(y)(x)*y

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

Base.inv(x::RingElem) = divexact(one(parent(x)), x)

function isdivisible_by(x::T, y::T) where T <: RingElem
   if iszero(y)
      return iszero(x)
   end
   r = rem(x, y)
   return iszero(r)
end

function ==(x::RingElem, y::RingElem)
  fl, u, v = try_promote(x, y)
  if fl
    return u == v
  else
    return false
  end
end

==(x::RingElem, y::RingElement) = x == parent(x)(y)

==(x::RingElement, y::RingElem) = parent(y)(x) == y

function addmul!(z::T, x::T, y::T, c::T) where {T <: RingElem}
   c = mul!(c, x, y)
   z = addeq!(z, c)
   return z
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(x::AbstractAlgebra.PolyElem{T}, y::Integer) where T <: RingElem
   return evaluate(x, base_ring(x)(y))
end

function evaluate(x::AbstractAlgebra.MPolyElem{T}, y::Integer) where T <: RingElem
   return evaluate(x, base_ring(x)(y))
end

###############################################################################
#
#   Delayed reduction
#
###############################################################################

Base.broadcastable(m::RingElem) = Ref(m)

###############################################################################
#
#   Delayed reduction
#
###############################################################################

# Fall back to ordinary multiplication
function mul_red!(a::T, b::T, c::T, flag::Bool) where T <: RingElement
   return mul!(a, b, c)
end

# Define addmul_delayed_reduction! for all ring elem types
function addmul_delayed_reduction!(a::T, b::T, c::T, d::T) where T <: RingElement
   d = mul_red!(d, b, c, false)
   return addeq!(a, d)
end

# Fall back to nop
function reduce!(a::RingElement)
   return a
end

###############################################################################
#
#   Type traits
#
###############################################################################

# Type can only represent elements of an exact ring
# true unless explicitly specified
isexact_type(R::Type{T}) where T <: RingElem = true

# Type can only represent elements of domains
# false unless explicitly specified
isdomain_type(R::Type{T}) where T <: RingElem = false

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

adjoint(x::T) where {T <: MatElem} = transpose(x)

adjoint(x::T) where {T <: RingElem} = deepcopy(x)

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

@doc Markdown.doc"""
    sqrt(a::FieldElem)

Return the square root of the element `a`. By default the function will
throw an exception if the input is not square. If `check=false` this test is
omitted.
"""
function Base.sqrt(a::FieldElem; check::Bool=true)
  R = parent(a)
  R, t = PolynomialRing(R, "t", cached = false)
  f = factor(t^2 - a)
  for (p, e) in f
    if !check || degree(p) == 1
      return -divexact(coeff(p, 0), coeff(p, 1); check=check)
    end
  end
  throw(error("Element $a does not have a square root"))
end

# assumes the existence of sqrt without check argument for input
function Base.sqrt(a::RingElem; check::Bool=true)
  s = sqrt(a)
  if check
    s != a^2 && error("Element $a does not have a square root")
  end
  return s
end  

# assumes the existence of issquare and sqrt for input  
function issquare_with_sqrt(a::RingElem)
  if issquare(a)
     return true, sqrt(a)
  else
     return false, parent(a)()
  end
end

###############################################################################
#
#   Generic and specific rings and fields
#
###############################################################################

include("julia/Integer.jl")

include("julia/Rational.jl")

include("julia/Float.jl")

include("Fields.jl")

include("Factor.jl")

###############################################################################
#
#   Generic functions to be defined after all rings
#
###############################################################################

include("polysubst.jl")

