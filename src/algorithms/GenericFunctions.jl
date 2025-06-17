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
   elseif isone(n)
      return deepcopy(a)
   elseif iszero(n)
      return one(parent(a))
   elseif n == -1
      return inv(a)
   else
      return internal_power(inv(a), -widen(n))
   end
end

###############################################################################
#
# The whole Euclidean interface can be derived from divrem if it is available.
#
###############################################################################

@doc raw"""
    divrem(f::T, g::T) where T <: RingElem

Return a pair `q, r` consisting of the Euclidean quotient and remainder of $f$
by $g$. A `DivideError` should be thrown if $g$ is zero.
"""
function Base.divrem(a::T, b::T) where T <: RingElem
  throw(NotImplementedError(:divrem, a, b))
end

@doc raw"""
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

@doc raw"""
    div(f::T, g::T) where T <: RingElem

Return the Euclidean quotient of $f$ by $g$. A `DivideError` should be thrown
if $g$ is zero.
"""
function Base.div(a::T, b::T) where T <: RingElem
   return divrem(a, b)[1]
end

@doc raw"""
    mulmod(f::T, g::T, m::T) where T <: RingElem

Return `mod(f*g, m)` but possibly computed more efficiently.
"""
function mulmod(a::T, b::T, m::T) where T <: RingElement
   return mod(a*b, m)
end

function mulmod(a::T, b::T, m::T) where T <: Integer
   return mod(widen(a)*b, m) % T
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

@doc raw"""
    powermod(f::T, e::Int, m::T) where T <: RingElem

Return `mod(f^e, m)` but possibly computed more efficiently.
"""
function powermod(a::T, n::Integer, m::T) where T <: RingElem
   check_parent(a, m)
   if n > 1
      return internal_powermod(a, n, m)
   elseif isone(n)
      return mod(a, m)
   elseif iszero(n)
      return mod(one(parent(a)), m)
   elseif n == -1
      return invmod(a, m)
   else
      return internal_powermod(invmod(a, m), -widen(n), m)
   end
end

@doc raw"""
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

@doc raw"""
    divides(f::T, g::T) where T <: RingElem

Return a pair, `flag, q`, where `flag` is set to `true` if $g$ divides $f$, in which
case `q` is set to the quotient, or `flag` is set to `false` and `q`
is undefined.
"""
function divides(a::T, b::T) where T <: RingElem
   check_parent(a, b)
   if iszero(b)
      return iszero(a), b
   end
   q, r = divrem(a, b)
   return iszero(r), q
end

@doc raw"""
    remove(f::T, p::T) where T <: RingElem

Return a pair `v, q` where $p^v$ is the highest power of $p$ dividing $f$ and $q$ is
the cofactor after $f$ is divided by this power.

See also [`valuation`](@ref), which only returns the valuation.
"""
function remove(a::T, b::T) where T <: Union{RingElem, Number}
   check_parent(a, b)
   if (iszero(b) || is_unit(b))
      throw(ArgumentError("Second argument must be a non-zero non-unit"))
   end
   if iszero(a)
      return (0, zero(parent(a))) # questionable case, consistent with ZZRingElem
   end
   v = 0
   while begin; (ok, q) = divides(a, b); ok; end
      a = q
      v += 1
   end
   return v, a
end

@doc raw"""
    valuation(f::T, p::T) where T <: RingElem

Return `v` where $p^v$ is the highest power of $p$ dividing $f$.

See also [`remove`](@ref).
"""
function valuation(a::T, b::T) where T <: Union{RingElem, Number}
   return remove(a, b)[1]
end

@doc raw"""
    gcd(a::T, b::T) where T <: RingElem

Return a greatest common divisor of $a$ and $b$, i.e., an element $g$
which is a common divisor of $a$ and $b$, and with the property that
any other common divisor of $a$ and $b$ divides $g$.

!!! note
    For best compatibility with the internal assumptions made by
    AbstractAlgebra, the return is expected to be unit-normalized in such a
    way that if the return is a unit, that unit should be one.
"""
function gcd(a::T, b::T) where T <: RingElem
   check_parent(a, b)
   while !iszero(b)
      (a, b) = (b, mod(a, b))
   end
   return iszero(a) ? a : divexact(a, canonical_unit(a))
end

@doc raw"""
    gcd(fs::AbstractArray{<:T}) where T <: RingElem

Return a greatest common divisor of the elements in `fs`.
Requires that `fs` is not empty.
"""
function gcd(fs::AbstractArray{<:T}) where T <: RingElem
   length(fs) > 0 || error("Empty collection")
   return reduce(gcd, fs)
end

@doc raw"""
    gcd(f::T, g::T, hs::T...) where T <: RingElem

Return a greatest common divisor of $f$, $g$ and the elements in `hs`.
"""
function gcd(f::T, g::T, hs::T...) where T <: RingElem
   return gcd(f, gcd(g, hs...))
end

@doc raw"""
    gcd_with_cofactors(a::T, b::T) where T <: RingElem

Return a tuple `(g, abar, bbar)` consisting of `g = gcd(a, b)` and cofactors
`abar` and `bbar` with `a = g*abar` and `b = g*bbar`.
"""
function gcd_with_cofactors(a::T, b::T) where T <: RingElement
   g = gcd(a, b)
   if iszero(g) || isone(g)
      return (g, a, b)
   else
      return (g, divexact(a, g), divexact(b, g))
   end
end

@doc raw"""
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

@doc raw"""
    lcm(fs::AbstractArray{<:T}) where T <: RingElem

Return a least common multiple of the elements in `fs`.
Requires that `fs` is not empty.
"""
function lcm(fs::AbstractArray{<:T}) where T <: RingElem
   length(fs) > 0 || error("Empty collection")
   return reduce(lcm, fs)
end

@doc raw"""
    lcm(f::T, g::T, hs::T...) where T <: RingElem

Return a least common multiple of $f$, $g$ and the elements in `hs`.
"""
function lcm(f::T, g::T, hs::T...) where T <: RingElem
   return lcm(f, lcm(g, hs...))
end

@doc raw"""
    gcdx(f::T, g::T) where T <: RingElem

Return a triple `d, s, t` such that $d = gcd(f, g)$ and $d = sf + tg$, with $s$
loosely reduced modulo $g/d$ and $t$ loosely reduced modulo $f/d$.
"""
function gcdx(a::T, b::T) where T <: RingElem
   check_parent(a, b)
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

@doc raw"""
    gcdinv(f::T, g::T) where T <: RingElem

Return a tuple `d, s` such that $d = gcd(f, g)$ and $s = (f/d)^{-1} \pmod{g/d}$. Note
that $d = 1$ iff $f$ is invertible modulo $g$, in which case $s = f^{-1} \pmod{g}$.
"""
function gcdinv(a::T, b::T) where T <: RingElem
   g, s, t = gcdx(a, b)
   return (g, s)
end

# TODO: Move from CRT from Hecke/src/Misc

function _crt_with_lcm_stub(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement
   diff = r2 - r1
   if iszero(m1)
      check && !is_divisible_by(diff, m2) && error("no crt solution")
      return (r1, m1)
   elseif iszero(m2)
      check && !is_divisible_by(diff, m1) && error("no crt solution")
      return (r2, m2)
   end
   # eliminating one cofactor computation with g, s = gcdinv(m1, m2) should be
   # sufficient, but almost all of Nemo's implementations of gcdinv are
   # non-conforming (i.e. they throw or return a wrong gcd)
   g, s, _ = gcdx(m1, m2)
   if isone(g)
      return (r1 + mulmod(diff, s, m2)*m1, m1*m2)
   elseif !check
      m1og = divexact(m1, g; check=false)
      return (r1 + mulmod(diff, s, m2)*m1og, m1og*m2)
   else
      m2og = divexact(m2, g; check=false)
      diff = divexact(diff, g; check=check)
      return (r1 + mulmod(diff, s, m2og)*m1, m1*m2og)
   end
end

function _crt_stub(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement
    return _crt_with_lcm_stub(r1, m1, r2, m2; check=check)[1]
end

@doc raw"""
    crt(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement

Return an element congruent to $r_1$ modulo $m_1$ and $r_2$ modulo $m_2$.
If `check = true` and no solution exists, an error is thrown.

If `T` is a fixed precision integer type (like `Int`), the result will be
correct if `abs(ri) <= abs(mi)` and `abs(m1 * m2) < typemax(T)`.
"""
function crt(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement
   return _crt_stub(r1, m1, r2, m2; check=check)
end

@doc raw"""
    crt_with_lcm(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement

Return a tuple consisting of an element congruent to $r_1$ modulo $m_1$ and
$r_2$ modulo $m_2$ and the least common multiple of $m_1$ and $m_2$.
If `check = true` and no solution exists, an error is thrown.
"""
function crt_with_lcm(r1::T, m1::T, r2::T, m2::T; check::Bool=true) where T <: RingElement
   return _crt_with_lcm_stub(r1, m1, r2, m2; check=check)
end

function _crt_with_lcm_stub(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement
   n = length(r)
   @assert n == length(m)
   @assert n > 0
   n < 2 && return (r[1], m[1])
   n == 2 && return crt_with_lcm(r[1], m[1], r[2], m[2]; check=check)
   return reduce((a, b) -> crt_with_lcm(a[1], a[2], b[1], b[2]; check=check),
                 ((r[i], m[i]) for i in 1:n))
end

function _crt_stub(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement
    return _crt_with_lcm_stub(r, m; check=check)[1]
end

@doc raw"""
    crt(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement

Return an element congruent to $r_i$ modulo $m_i$ for each $i$.
"""
function crt(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement
   return _crt_stub(r, m; check=check)
end

@doc raw"""
    crt_with_lcm(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement

Return a tuple consisting of an element congruent to $r_i$ modulo $m_i$ for
each $i$ and the least common multiple of the $m_i$.
"""
function crt_with_lcm(r::Vector{T}, m::Vector{T}; check::Bool=true) where T <: RingElement
   return _crt_with_lcm_stub(r, m; check=check)
end

###############################################################################
#
# Functions that can't really be implemented generically
#
###############################################################################

@doc raw"""
    is_zero_divisor(a::T) where T <: RingElement

Return `true` if there exists a nonzero $b$ such that $a b = 0$ and
`false` otherwise.
"""
function is_zero_divisor(a::T) where T <: RingElement
   if !is_domain_type(T)
      throw(NotImplementedError(:is_zero_divisor, a))
   end
   return is_zero(a) && !is_trivial(parent(a))
end

ConformanceTests._implements(::Type{T}, f::typeof(is_zero_divisor)) where {T} = is_domain_type(T) || _implements_directly(T, f)

@doc raw"""
    is_zero_divisor_with_annihilator(a::T) where T <: RingElement

Return `(true, b)` if there exists a nonzero $b$ such that $a b = 0$ and
`(false, junk)` otherwise.
"""
function is_zero_divisor_with_annihilator(a::T) where T <: RingElement
   if !is_domain_type(T)
      if is_zero_divisor(a)
         throw(NotImplementedError(:is_zero_divisor_with_annihilator, a))
      end
      return (false, parent(a)())
   end
   theone = one(parent(a))
   return (is_zero(a) && !is_zero(theone), theone)
end

@doc raw"""
    factor(a::T) where T <: RingElement -> Fac{T}

Return a factorization of $a$ into irreducible elements, as a `Fac{T}`.
The irreducible elements in the factorization are pairwise coprime.
"""
function factor(a)
   throw(NotImplementedError(:factor, a))
end

ConformanceTests._implements(::Type{T}, f::typeof(factor)) where {T} = _implements_directly(T, f)

@doc raw"""
    factor_squarefree(a::T) where T <: RingElement -> Fac{T}

Return a factorization of $a$ into squarefree elements, as a `Fac{T}`.
The squarefree elements in the factorization are pairwise coprime.
"""
function factor_squarefree(a)
   throw(NotImplementedError(:factor_squarefree, a))
end

ConformanceTests._implements(::Type{T}, f::typeof(factor_squarefree)) where {T} = _implements_directly(T, f)

@doc raw"""
    is_irreducible(a::RingElement)

Return `true` if $a$ is irreducible, else return `false`.
Zero and units are by definition never irreducible.
"""
function is_irreducible(a)
   is_zero(a) && return false
   is_unit(a) && return false
   af = factor(a)
   return length(af) == 1 && all(isone, values(af.fac))
end

ConformanceTests._implements(::Type{T}, ::typeof(is_irreducible)) where {T} = _implements(T, is_unit) && _implements(T, factor)

@doc raw"""
    is_squarefree(a::RingElement)

Return `true` if $a$ is squarefree, else return `false`.
An element is squarefree if it it is not divisible by any squares
except the squares of units.
"""
function is_squarefree(a)
   iszero(a) && return false
   is_unit(a) && return true
   af = factor_squarefree(a)
   return all(isone, values(af.fac))
end

ConformanceTests._implements(::Type{T}, ::typeof(is_squarefree)) where {T} = _implements(T, is_unit) && _implements(T, factor_squarefree)
