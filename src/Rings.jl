###############################################################################
#
#   Rings.jl : Generic commutative rings
#
###############################################################################

function isequal(a::RingElem, b::RingElem)
   return parent(a) == parent(b) && a == b
end

###############################################################################
#
#  Element/Parent types for instances
#
###############################################################################

@doc Markdown.doc"""
    elem_type(parent)
    elem_type(parent_type)
Given a parent object (or its type), return the type of its elements.
"""
elem_type(x)  = elem_type(typeof(x))
elem_type(T::DataType) = throw(MethodError(elem_type, (T,)))

@doc Markdown.doc"""
    parent_type(element)
    parent_type(element_type)
Given an element (or its type), return the type of its parent object.
"""
parent_type(x) = parent_type(typeof(x))
parent_type(T::DataType) = throw(MethodError(parent_type, (T,)))

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

###############################################################################
#
#   Generic catchall functions
#
###############################################################################

function +(x::S, y::T) where {S <: RingElem, T <: RingElem}
   if S == promote_rule(S, T)
      +(x, parent(x)(y))
   else
      +(parent(y)(x), y)
   end
end

+(x::RingElem, y::RingElement) = x + parent(x)(y)

+(x::RingElement, y::RingElem) = parent(y)(x) + y

function -(x::S, y::T) where {S <: RingElem, T <: RingElem}
   if S == promote_rule(S, T)
      -(x, parent(x)(y))
   else
      -(parent(y)(x), y)
   end
end

-(x::RingElem, y::RingElement) = x - parent(x)(y)

-(x::RingElement, y::RingElem) = parent(y)(x) - y

function *(x::S, y::T) where {S <: RingElem, T <: RingElem}
   if S == promote_rule(S, T)
      *(x, parent(x)(y))
   else
      *(parent(y)(x), y)
   end
end

*(x::RingElem, y::RingElement) = x*parent(x)(y)

*(x::RingElement, y::RingElem) = parent(y)(x)*y


"""
    divexact(x, y)

Return an exact quotient of `x` by `y`, i.e. an element
`z` such that `x == yz`; when `x` and `y` do not belong to the same ring,
they are first coerced into a common ring.
If no exact division is possible, an exception is raised.
"""
function divexact end

function divexact(x::S, y::T) where {S <: RingElem, T <: RingElem}
   if S == promote_rule(S, T)
      divexact(x, parent(x)(y))
   else
      divexact(parent(y)(x), y)
   end
end

divexact(x::RingElem, y::RingElement) = divexact(x, parent(x)(y))

divexact(x::RingElement, y::RingElem) = divexact(parent(y)(x), y)

divexact_left(x::T, y::T) where T <: RingElement = divexact(x, y)

divexact_right(x::T, y::T) where T <: RingElement = divexact(x, y)

Base.inv(x::RingElem) = divexact(one(parent(x)), x)

function divides(x::T, y::T) where {T <: RingElem}
   if iszero(y)
      return iszero(x), y
   end
   q, r = divrem(x, y)
   return iszero(r), q
end

function ==(x::S, y::T) where {S <: RingElem, T <: RingElem}
   R = promote_rule(S, T)
   if S == R
      ==(x, parent(x)(y))
   elseif T == R
      ==(parent(y)(x), y)
   else
      false
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

Return the square root of the element `a`.
"""
function Base.sqrt(a::FieldElem)
  R = parent(a)
  R, t = PolynomialRing(R, "t", cached = false)
  f = factor(t^2 - a)
  for (p, e) in f
    if degree(p) == 1
      return -divexact(coeff(p, 0), coeff(p, 1))
    end
  end
  throw(error("Element $a does not have a square root"))
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

