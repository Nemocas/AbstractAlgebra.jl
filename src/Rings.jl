###############################################################################
#
#   Rings.jl : Generic commutative rings
#
###############################################################################

elem_type(::T) where {T <: Ring} = elem_type(T)

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

function divexact(x::S, y::T) where {S <: RingElem, T <: RingElem}
   if S == promote_rule(S, T)
      divexact(x, parent(x)(y))
   else
      divexact(parent(y)(x), y)
   end
end

divexact(x::RingElem, y::RingElement) = divexact(x, parent(x)(y))

divexact(x::RingElement, y::RingElem) = divexact(parent(y)(x), y)

function divides(x::T, y::T) where {T <: RingElem}
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

Base.literal_pow(::typeof(^), x::T, ::Val{p}) where {p, T <: RingElem} = x^p

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
#   Baby-steps giant-steps powering
#
###############################################################################

function powers(a::T, d::Int) where {T <: RingElement}
   d <= 0 && throw(DomainError())
   S = parent(a)
   A = Array{T}(undef, d + 1)
   A[1] = one(S)
   if d > 0
      c = a
      A[2] = a
      for i = 2:d
         c *= a
         A[i + 1] = c
      end
   end
   return A
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

function exp(a::T) where {T <: RingElem}
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
#   One and zero
#
###############################################################################

one(x::T) where {T <: RingElem} = one(parent(x))

zero(x::T) where {T <: RingElem} = zero(parent(x))

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
#  Change base ring
#
################################################################################

@doc Markdown.doc"""
    change_base_ring(::PolyElem{T}, g::Any) where T <: RingElement
> Return the polynomial obtained by applying `g` to the coefficients. The new
> base ring is defined by the image of `0`. 
"""
function change_base_ring(p::PolyElem{T}, g) where T <: RingElement
   z = zero(base_ring(parent(p)))
   new_base_ring = parent(g(z))
   new_var_name = string(var(parent(p)))
   P, _ = PolynomialRing(new_base_ring, new_var_name)
   return change_base_ring(p, g, P)
end

@doc Markdown.doc"""
    change_base_ring(::PolyElem{T}, g::Any, Rx::PolyRing) where T <: RingElement
> Return the polynomial obtained by applying `g` to the coefficients. The result
> will have parent `Rx`.
"""
function change_base_ring(p::PolyElem{T}, g, Rx::PolyRing) where T <: RingElement
   new_coefficients = [g(coeff(p, i)) for i in 0:degree(p)]
   return Rx(new_coefficients)
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
