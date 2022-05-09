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

@doc Markdown.doc"""
    divrem(f::T, g::T) where T <: RingElem

Return a pair `q, r` consisting of the Euclidean quotient and remainder of $f$
by $g$. A `DivideError` should be thrown if $g$ is zero.
"""
function divrem end

@doc Markdown.doc"""
    mod(f::T, g::T) where T <: RingElem

Return the Euclidean remainder of $f$ by $g$. A `DivideError` should be thrown
if $g$ is zero.

!!! note
    For best compatibility with the internal assumptions made by AbstractAlgebra,
    the Euclidean remainder function should provide unique representatives for
    the residue classes; the `mod` function should satisfy

    1. `mod(a_1, b) = mod(a_2, b)` if and only if $b$ divides $a_1 - a_2$, and
    2. `mod(0, b) = 0`.
"""
function mod(a::T, b::T) where T <: RingElem
   return divrem(a, b)[2]
end

@doc Markdown.doc"""
    div(f::T, g::T) where T <: RingElem

Return the Euclidean quotient of $f$ by $g$. A `DivideError` should be thrown
if $g$ is zero.
"""
function Base.div(a::T, b::T) where T <: RingElem
   return divrem(a, b)[1]
end

@doc Markdown.doc"""
    mulmod(f::T, g::T, m::T) where T <: RingElem

Return `mod(f*g, m)` but possibly computed more efficiently.
"""
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

@doc Markdown.doc"""
    powermod(f::T, e::Int, m::T) where T <: RingElem

Return `mod(f^e, m)` but possibly computed more efficiently.
"""
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

@doc Markdown.doc"""
    invmod(f::T, m::T) where T <: RingElem

Return an inverse of $f$ modulo $m$, meaning that `isone(mod(invmod(f,m)*f,m))`
returns `true`.

If such an inverse doesn't exist, a `NotInvertibleError` should be thrown.
"""
function invmod(a::T, m::T) where T <: RingElem
   g, s = gcdinv(a, m)
   isone(g) || throw(NotInvertibleError(a, m))
   return mod(s, m)  # gcdinv has no canonicity requirement on s
end

@doc Markdown.doc"""
    divides(f::T, g::T) where T <: RingElem

Return a pair, `flag, q`, where `flag` is set to `true` if $g$ divides $f$, in which
case `q` is set to the quotient, or `flag` is set to `false` and `q`
is set to `zero(f)`.
"""
function divides(a::T, b::T) where T <: RingElem
   parent(a) == parent(b) || error("Incompatible parents")
   if iszero(b)
      return iszero(a), b
   end
   q, r = divrem(a, b)
   return iszero(r), q
end

@doc Markdown.doc"""
    remove(f::T, p::T) where T <: RingElem

Return a pair `v, q` where $p^v$ is the highest power of $p$ dividing $f$ and $q$ is
the cofactor after $f$ is divided by this power.

See also [`valuation`](@ref), which only returns the valuation.
"""
function remove(a::T, b::T) where T <: Union{RingElem, Number}
   parent(a) == parent(b) || error("Incompatible parents")
   if (iszero(b) || is_unit(b))
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

@doc Markdown.doc"""
    valuation(f::T, p::T) where T <: RingElem

Return `v` where $p^v$ is the highest power of $p$ dividing $f$.

See also [`remove`](@ref).
"""
function valuation(a::T, b::T) where T <: Union{RingElem, Number}
   return remove(a, b)[1]
end

@doc Markdown.doc"""
    gcd(f::T, g::T) where T <: RingElem

Return a greatest common divisor of $f$ and $g$, i.e., an element $d$
which is a common divisor of $f$ and $g$, and with the property that
any other common divisor of $f$ and $g$ divides $d$.

!!! note
    For best compatibility with the internal assumptions made by
    AbstractAlgebra, the return is expected to be unit-normalized in such a
    way that if the return is a unit, that unit should be one.
"""
function gcd(a::T, b::T) where T <: RingElem
   parent(a) == parent(b) || error("Incompatible parents")
   while !iszero(b)
      (a, b) = (b, mod(a, b))
   end
   return iszero(a) ? a : divexact(a, canonical_unit(a))
end

@doc Markdown.doc"""
    lcm(f::T, g::T) where T <: RingElem

Return a least common multiple of $f$ and $g$, i.e., an element $d$
which is a common multiple of $f$ and $g$, and with the property that
any other common multiple of $f$ and $g$ is a multiple of $d$.
"""
function lcm(a::T, b::T) where T <: RingElem
   g = gcd(a, b)
   iszero(g) && return g
   return a*divexact(b, g)
end

@doc Markdown.doc"""
    gcdx(f::T, g::T) where T <: RingElem

Return a triple `d, s, t` such that $d = gcd(f, g)$ and $d = sf + tg$, with $s$
loosely reduced modulo $g/d$ and $t$ loosely reduced modulo $f/d$.
"""
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

@doc Markdown.doc"""
    gcdinv(f::T, g::T) where T <: RingElem

Return a tuple `d, s` such that $d = gcd(f, g)$ and $s = (f/d)^{-1} \pmod{g/d}$. Note
that $d = 1$ iff $f$ is invertible modulo $g$, in which case $s = f^{-1} \pmod{g}$.
"""
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
    is_irreducible(a)

Return `true` if $a$ is irreducible, else return `false`.
"""
function is_irreducible end

@doc Markdown.doc"""
    is_squarefree(a)

Return `true` if $a$ is squarefree, else return `false`.
"""
function is_squarefree end

