###############################################################################
#
#   GenericFunctions.jl : Functions for Generic types
#
###############################################################################

function internal_power(a, n)
   @assert n > 1
   while iseven(n)
      a = a*a
      n >>= 1
   end
   z = a
   while !iszero(n >>= 1)
      a = a*a
      if isodd(n)
         z = z*a
      end
   end
   return z
end

function ^(a::T, n::Integer) where T <: RingElem
   if n > 1
      return internal_power(a, n)
   elseif n == 1
      return deepcopy(a)
   elseif n == 0
      return one(parent(a))
   elseif n == -1
      return inv(a)
   else
      return internal_power(inv(a), -widen(n))
   end
end

###############################################################################
#
# Euclidean interface has 13 functions, 12 of which follow from divrem
#
###############################################################################

function mod(a::T, b::T) where T <: RingElem
   return divrem(a, b)[2]
end

function Base.div(a::T, b::T) where T <: RingElem
   return divrem(a, b)[1]
end

function mulmod(a::T, b::T, m::T) where T <: RingElement
   return mod(a*b, m)
end

function internal_powermod(a, n, m)
   @assert n > 1
   while iseven(n)
      a = mulmod(a, a, m)
      n >>= 1
   end
   z = a
   while !iszero(n >>= 1)
      a = mulmod(a, a, m)
      if isodd(n)
         z = mulmod(z, a, m)
      end
   end
   return z
end

function powermod(a::T, n::Integer, m::T) where T <: RingElem
   parent(a) == parent(m) || error("Incompatible parents")
   if n > 1
      return internal_powermod(a, n, m)
   elseif n == 1
      return mod(a, m)
   elseif n == 0
      return mod(one(parent(a)), m)
   elseif n == -1
      return invmod(a, m)
   else
      return internal_powermod(invmod(a, m), -widen(n), m)
   end
end

function invmod(a::T, m::T) where T <: RingElem
   g, s = gcdinv(a, m)
   isone(g) || throw(NotInvertibleError(a, m))
   return mod(s, m)  # gcdinv has no canonicity requirement on s
end

function divides(a::T, b::T) where T <: RingElem
   parent(a) == parent(b) || error("Incompatible parents")
   if iszero(b)
      return iszero(a), b
   end
   q, r = divrem(a, b)
   return iszero(r), q
end

function remove(a::T, b::T) where T <: Union{RingElem, Number}
   parent(a) == parent(b) || error("Incompatible parents")
   if (iszero(b) || isunit(b))
      throw(ArgumentError("Second argument must be a non-zero non-unit"))
   end
   if iszero(a)
      return (0, zero(parent(a))) # questionable case, consistent with fmpz
   end
   v = 0
   while begin; (ok, q) = divides(a, b); ok; end
      a = q
      v += 1
   end
   return v, a
end

function valuation(a::T, b::T) where T <: Union{RingElem, Number}
   return remove(a, b)[1]
end

function gcd(a::T, b::T) where T <: RingElem
   parent(a) == parent(b) || error("Incompatible parents")
   while !iszero(b)
      (a, b) = (b, mod(a, b))
   end
   return iszero(a) ? a : divexact(a, canonical_unit(a))
end

function lcm(a::T, b::T) where T <: RingElem
   g = gcd(a, b)
   iszero(g) && return g
   return a*divexact(b, g)
end

function gcdx(a::T, b::T) where T <: RingElem
   parent(a) == parent(b) || error("Incompatible parents")
   R = parent(a)
   if iszero(a)
      if iszero(b)
         return zero(R), zero(R), zero(R)
      else
         t = canonical_unit(b)
         return divexact(b, t), zero(R), inv(t)
      end
   elseif iszero(b)
      t = canonical_unit(a)
      return divexact(a, t), inv(t), zero(R)
   end
   m11, m12 = one(R), zero(R)
   m21, m22 = zero(R), one(R)
   while !iszero(b)
      (q, b), a = divrem(a, b), b
      m11, m12 = m12, m11 - q*m12
      m21, m22 = m22, m21 - q*m22
   end
   t = canonical_unit(a)
   return divexact(a, t), divexact(m11, t), divexact(m21, t)
end

function gcdinv(a::T, b::T) where T <: RingElem
   g, s, t = gcdx(a, b)
   return (g, s)
end


# TODO: Move from CRT from Hecke to AbstractAlgebra?
# Currently no implementation, only example on how the arbitrary inputs `crt`
# should look like.
# @doc Markdown.doc"""
#     crt(r::AbstractVector{T}, m::AbstractVector{T}) where T
#     crt(r::T, m::T...) where T

# Return $x$ in the Euclidean domain $T$ such that $x \equiv r_i \mod m_i$
# for all $i$.
# """
function crt end

@doc Markdown.doc"""
    factor(a::T)

Return a factorization of the element $a$ as a `Fac{T}`.
"""
function factor end

@doc Markdown.doc"""
    factor_squarefree(a::T)

Return a squarefree factorization of the element $a$ as a `Fac{T}`.
"""
function factor_squarefree end

@doc Markdown.doc"""
    isirreducible(a)

Return `true` if $a$ is irreducible, else return `false`.
"""
function isirreducible end

@doc Markdown.doc"""
    issquarefree(a)

Return `true` if $a$ is squarefree, else return `false`.
"""
function issquarefree end

