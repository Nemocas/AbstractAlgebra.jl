###############################################################################
#
#   NCRings.jl : Generic not necessarily commutative rings
#
###############################################################################

elem_type(::T) where {T <: NCRing} = elem_type(T)

include("Rings.jl")

###############################################################################
#
#   Promotion system
#
###############################################################################

promote_rule(::Type{T}, ::Type{T}) where T <: NCRingElem = T

###############################################################################
#
#   Generic catchall functions
#
###############################################################################

function +(x::S, y::T) where {S <: NCRingElem, T <: NCRingElem}
   if S == promote_rule(S, T)
      +(x, parent(x)(y))
   else
      +(parent(y)(x), y)
   end
end

+(x::NCRingElem, y::RingElement) = x + parent(x)(y)

+(x::RingElement, y::NCRingElem) = parent(y)(x) + y

function -(x::S, y::T) where {S <: NCRingElem, T <: NCRingElem}
   if S == promote_rule(S, T)
      -(x, parent(x)(y))
   else
      -(parent(y)(x), y)
   end
end

-(x::NCRingElem, y::RingElement) = x - parent(x)(y)

-(x::RingElement, y::NCRingElem) = parent(y)(x) - y

function *(x::S, y::T) where {S <: NCRingElem, T <: NCRingElem}
   if S == promote_rule(S, T)
      *(x, parent(x)(y))
   else
      *(parent(y)(x), y)
   end
end

*(x::NCRingElem, y::RingElement) = x*parent(x)(y)

*(x::RingElement, y::NCRingElem) = parent(y)(x)*y

function divexact_left(x::S, y::T) where {S <: NCRingElem, T <: NCRingElem}
   if S == promote_rule(S, T)
      divexact_left(x, parent(x)(y))
   else
      divexact_left(parent(y)(x), y)
   end
end

function divexact_right(x::S, y::T) where {S <: NCRingElem, T <: NCRingElem}
   if S == promote_rule(S, T)
      divexact_right(x, parent(x)(y))
   else
      divexact_right(parent(y)(x), y)
   end
end

function divexact_left(x::NCRingElem, y::Union{Integer, Rational, AbstractFloat})
   return divexact_left(x, parent(x)(y))
end

function divexact_right(x::NCRingElem, y::Union{Integer, Rational, AbstractFloat})
   return divexact_right(x, parent(x)(y))
end

function ==(x::S, y::T) where {S <: NCRingElem, T <: NCRingElem}
   if S == promote_rule(S, T)
      ==(x, parent(x)(y))
   else
      ==(parent(y)(x), y)
   end
end

==(x::NCRingElem, y::RingElement) = x == parent(x)(y)

==(x::RingElement, y::NCRingElem) = parent(y)(x) == y

###############################################################################
#
#   Baby-steps giant-steps powering
#
###############################################################################

function powers(a::T, d::Int) where {T <: NCRingElem}
   d <= 0 && throw(DomainError())
   S = parent(a)
   A = Array{T}(undef, d + 1)
   A[1] = one(S)
   if d > 1
      c = a
      A[2] = a
      for i = 2:d
         c *= a
         A[i + 1] = c
      end
   end
   return A
end

