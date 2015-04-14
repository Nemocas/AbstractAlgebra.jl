###############################################################################
#
#   ZZ.jl : BigInts
#
###############################################################################

# Copyright (c) 2009-2014: Jeff Bezanson, Stefan Karpinski, Viral B. Shah,
# and other contributors:
# 
# https://github.com/JuliaLang/julia/contributors
#
# Copyright (C) 2014, William Hart
# 
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import Base: convert, promote_rule, show, string, parseint, serialize,
             deserialize, base, bin, dec, oct, hex, gcd, gcdx, lcm, div, size,
             zero, one, sign, hash

export ZZ, IntegerRing, parent, show, convert, hash, fac, binom, isprime, fdiv,
       cdiv, tdiv, div, rem, mod, gcd, xgcd, lcm, invmod, powmod, abs, divrem,
       isqrt, popcount, prevpow2, nextpow2, ndigits, dec, bin, oct, hex, base,
       one, zero, divexact, fits, sign, nbits, deepcopy, tdivpow2, fdivpow2,
       cdivpow2, flog, clog, cmpabs, clrbit!, setbit!, combit!, crt, divisible,
       divisor_lenstra, fdivrem, tdivrem, fmodpow2, gcdinv, isprobabprime, 
       issquare, jacobi, remove, root, size, isqrtrem, sqrtmod, trailing_zeros,
       sigma, eulerphi, fib, moebiusmu, primorial, risingfac, canonical_unit,
       needs_parentheses, is_negative, show_minus_one, parseint, addeq!, mul!,
       isunit, words

###############################################################################
#
#   Data type and memory management
#
###############################################################################

type IntegerRing <: Ring
end

type fmpz <: RingElem
    d::Int

    function fmpz()
        z = new()
        ccall((:fmpz_init, :libflint), Void, (Ptr{fmpz},), &z)
        finalizer(z, _fmpz_clear_fn)
        return z
    end

    function fmpz(x::Int)
        z = new()
        ccall((:fmpz_init_set_si, :libflint), Void, (Ptr{fmpz}, Int), &z, x)
        finalizer(z, _fmpz_clear_fn)
        return z
    end

    function fmpz(x::BigInt)
        z = new()
        ccall((:fmpz_init, :libflint), Void, (Ptr{fmpz},), &z)
        ccall((:fmpz_set_mpz, :libflint), Void, (Ptr{fmpz}, Ptr{BigInt}), &z, &x)
        finalizer(z, _fmpz_clear_fn)
        return z
    end

    fmpz(x::fmpz) = x
end

function _fmpz_clear_fn(a::fmpz)
   ccall((:fmpz_clear, :libflint), Void, (Ptr{fmpz},), &a)
end

ZZ = IntegerRing()

parent(a::fmpz) = ZZ

elem_type(::IntegerRing) = fmpz

base_ring(a::IntegerRing) = None

hash(a::fmpz) = hash(BigInt(a))

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function deepcopy(a::fmpz)
   z = ZZ()
   ccall((:fmpz_set, :libflint), Void, (Ptr{fmpz}, Ptr{fmpz}), &z, &a)
   return z
end

one(::IntegerRing) = ZZ(1)

zero(::IntegerRing) = ZZ(0)

sign(a::fmpz) = int(ccall((:fmpz_sgn, :libflint), Cint, (Ptr{fmpz},), &a))

fits(::Type{Int}, a::fmpz) = ccall((:fmpz_fits_si, :libflint), Bool, 
                                   (Ptr{fmpz},), &a)

fits(::Type{Uint}, a::fmpz) = sign(a) < 0 ? false : 
              ccall((:fmpz_abs_fits_ui, :libflint), Bool, (Ptr{fmpz},), &a)

size(a::fmpz) = int(ccall((:fmpz_size, :libflint), Cint, (Ptr{fmpz},), &a))

isunit(a::fmpz) = ccall((:fmpz_is_pm1, :libflint), Bool, (Ptr{fmpz},), &a)

iszero(a::fmpz) = ccall((:fmpz_is_zero, :libflint), Bool, (Ptr{fmpz},), &a)

isone(a::fmpz) = ccall((:fmpz_is_one, :libflint), Bool, (Ptr{fmpz},), &a)

function words(a::fmpz)
   return a == 0 ? 0 : div(ndigits(a, 2) + 8*sizeof(Int) - 1, 8*sizeof(Int))
end

###############################################################################
#
#   Serialisation
#
###############################################################################

function serialize(s, n::fmpz)
    Base.serialize_type(s, fmpz)
    serialize(s, base(62, n))
end

deserialize(s, ::Type{fmpz}) = Base.parseint_nocheck(ZZ, deserialize(s), 62)

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::fmpz) = x < 0 ? ZZ(-1) : ZZ(1)

###############################################################################
#
#   Binary operators and functions
#
###############################################################################

# Metaprogram to define functions +, -, *, gcd, lcm, 
#                                 &, |, $

for (fJ, fC) in ((:+, :add), (:-,:sub), (:*, :mul),
                 (:gcd, :gcd), (:lcm, :lcm),
                 (:&, :and), (:|, :or), (:$, :xor))
    @eval begin
        function ($fJ)(x::fmpz, y::fmpz)
            z = ZZ()
            ccall(($(string(:fmpz_, fC)), :libflint), Void, 
                  (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &z, &x, &y)
            return z
        end
    end
end

# Metaprogram to define functions fdiv, cdiv, tdiv, div, mod

for (fJ, fC) in ((:fdiv, :fdiv_q), (:cdiv, :cdiv_q), (:tdiv, :tdiv_q), 
                 (:div, :tdiv_q), (:mod, :mod))
    @eval begin
        function ($fJ)(x::fmpz, y::fmpz)
            y == 0 && throw(DivideError())
            z = ZZ()
            ccall(($(string(:fmpz_, fC)), :libflint), Void, 
                  (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &z, &x, &y)
            return z
        end
    end
end

function divexact(x::fmpz, y::fmpz)
    y == 0 && throw(DivideError())
    z = ZZ()
    ccall((:fmpz_divexact, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &z, &x, &y)
    z
end

function flog(x::fmpz, c::fmpz)
    c <= 0 && throw(DomainError())
    x <= 0 && throw(DomainError())
    return ccall((:fmpz_flog, :libflint), Int, 
                 (Ptr{fmpz}, Ptr{fmpz}), &x, &c)
end

function clog(x::fmpz, c::fmpz)
    c <= 0 && throw(DomainError())
    x <= 0 && throw(DomainError())
    return ccall((:fmpz_clog, :libflint), Int, 
                 (Ptr{fmpz}, Ptr{fmpz}), &x, &c)
end

function %(x::fmpz, c::fmpz)
    c == 0 && throw(DivideError())
    q = ZZ()
    r = ZZ()
    ccall((:fmpz_tdiv_qr, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &q, &r, &x, &c)
    return r
end

rem(x::fmpz, c::fmpz) = %(x, c)

###############################################################################
#
#   Unsafe operators
#
###############################################################################

function mul!(z::fmpz, x::fmpz, y::fmpz)
   ccall((:fmpz_mul, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &z, &x, &y)
end

function addeq!(z::fmpz, x::fmpz)
   ccall((:fmpz_add, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &z, &z, &x)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function +(x::fmpz, c::Int)
    z = ZZ()
    if c >= 0
       ccall((:fmpz_add_ui, :libflint), Void, 
             (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    else
       ccall((:fmpz_sub_ui, :libflint), Void, 
             (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, -c)
    end
    return z
end

+(c::Int, x::fmpz) = x + c

function -(x::fmpz, c::Int)
    z = ZZ()
    if c >= 0
       ccall((:fmpz_sub_ui, :libflint), Void, 
             (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    else
       ccall((:fmpz_add_ui, :libflint), Void, 
             (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, -c)
    end
    return z
end

function -(c::Int, x::fmpz)
    z = ZZ()
    if c >= 0
       ccall((:fmpz_sub_ui, :libflint), Void, 
             (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    else
       ccall((:fmpz_add_ui, :libflint), Void, 
             (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, -c)
    end
    ccall((:__fmpz_neg, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}), &z, &z)
    return z
end

function *(x::fmpz, c::Int)
    z = ZZ()
    ccall((:fmpz_mul_si, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    return z
end

*(c::Int, x::fmpz) = x * c

###############################################################################
#
#   Shifting
#
###############################################################################

function <<(x::fmpz, c::Int)
    c < 0 && throw(DomainError())
    c == 0 && return x
    z = ZZ()
    ccall((:fmpz_mul_2exp, :libflint), Void,
          (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    return z
end

function >>(x::fmpz, c::Int)
    c < 0 && throw(DomainError())
    c == 0 && return x
    z = ZZ()
    ccall((:fmpz_fdiv_q_2exp, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    return z
end

###############################################################################
#
#   Ad hoc division
#
###############################################################################

function %(x::fmpz, c::Int)
   c < 0 && throw(DomainError())
   c == 0 && throw(DivideError())
   r = ccall((:fmpz_tdiv_ui, :libflint), Int, (Ptr{fmpz}, Int), &x, c)
   return sign(x) < 0 ? -r : r
end

rem(x::fmpz, c::Int) = %(x, c)

function tdivpow2(x::fmpz, c::Int)
    c < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_tdiv_q_2exp, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    return z
end

function fdivpow2(x::fmpz, c::Int)
    c < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_fdiv_q_2exp, :libflint), Void,
          (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    return z
end

function fmodpow2(x::fmpz, c::Int)
    c < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_fdiv_r_2exp, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    return z
end

function cdivpow2(x::fmpz, c::Int)
    c < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_cdiv_q_2exp, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    return z
end

function div(x::fmpz, c::Int)
    c == 0 && throw(DivideError())
    z = ZZ()
    ccall((:fmpz_tdiv_q_si, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    return z
end

function tdiv(x::fmpz, c::Int)
    c == 0 && throw(DivideError())
    z = ZZ()
    ccall((:fmpz_tdiv_q_si, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    return z
end

function fdiv(x::fmpz, c::Int)
    c == 0 && throw(DivideError())
    z = ZZ()
    ccall((:fmpz_fdiv_q_si, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    return z
end

function cdiv(x::fmpz, c::Int)
    c == 0 && throw(DivideError())
    z = ZZ()
    ccall((:fmpz_cdiv_si, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, c)
    return z
end

function mod(x::fmpz, c::Int)
    c < 0 && throw(DomainError())
    c == 0 && throw(DivideError())
    ccall((:fmpz_fdiv_ui, :libflint), Int, (Ptr{fmpz}, Int), &x, c)
end

function divexact(x::fmpz, y::Int)
    y == 0 && throw(DivideError())
    z = ZZ()
    ccall((:fmpz_divexact_si, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, y)
    z
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::fmpz, y::Int)
    if y < 0; throw(DomainError()); end
    if x == 1; return x; end
    if x == -1; return isodd(y) ? x : -x; end
    if y > typemax(Uint); throw(DomainError()); end
    if y == 0; return one(ZZ); end
    if y == 1; return x; end
    return x^uint(y)
end

###############################################################################
#
#   Unary operators and functions, e.g. -ZZ(12), ~ZZ(12)
#
###############################################################################

function -(x::fmpz)
    z = ZZ()
    ccall((:__fmpz_neg, :libflint), Void, (Ptr{fmpz}, Ptr{fmpz}), &z, &x)
    return z
end

function ~(x::fmpz)
    z = ZZ()
    ccall((:fmpz_complement, :libflint), Void, (Ptr{fmpz}, Ptr{fmpz}), &z, &x)
    return z
end

function abs(x::fmpz)
    z = ZZ()
    ccall((:fmpz_abs, :libflint), Void, (Ptr{fmpz}, Ptr{fmpz}), &z, &x)
    return z
end

###############################################################################
#
#   Division with remainder
#
###############################################################################

function divrem(x::fmpz, y::fmpz)
    y == 0 && throw(DivideError())
    z1 = ZZ()
    z2 = ZZ()
    ccall((:fmpz_tdiv_qr, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &z1, &z2, &x, &y)
    z1, z2
end

function tdivrem(x::fmpz, y::fmpz)
    y == 0 && throw(DivideError())
    z1 = ZZ()
    z2 = ZZ()
    ccall((:fmpz_tdiv_qr, :libflint), Void,
          (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &z1, &z2, &x, &y)
    z1, z2
end

function fdivrem(x::fmpz, y::fmpz)
    y == 0 && throw(DivideError())
    z1 = ZZ()
    z2 = ZZ()
    ccall((:fmpz_fdiv_qr, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &z1, &z2, &x, &y)
    z1, z2
end

###############################################################################
#
#   Modular arithmetic
#
###############################################################################

function powmod(x::fmpz, p::fmpz, m::fmpz)
    m <= 0 && throw(DomainError())
    if p < 0
       x = invmod(x, m)
       p = -p
    end
    r = ZZ()
    ccall((:fmpz_powm, :libflint), Void,
          (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}),
          &r, &x, &p, &m)
    return r 
end

function powmod(x::fmpz, p::Int, m::fmpz)
    m <= 0 && throw(DomainError())
    if p < 0
       x = invmod(x, m)
       p = -p
    end
    r = ZZ()
    ccall((:fmpz_powm_ui, :libflint), Void,
          (Ptr{fmpz}, Ptr{fmpz}, Int, Ptr{fmpz}),
          &r, &x, p, &m)
    return r 
end

function invmod(x::fmpz, y::fmpz)
    y <= 0 && throw(DomainError())
    z = ZZ()
    if y == 1
        return ZZ(0)
    end
    if ccall((:fmpz_invmod, :libflint), Cint, 
             (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &z, &x, &y) == 0
       error("Impossible inverse in invmod")
    end
    return z
end

function gcdinv(x::fmpz, y::fmpz)
    y <= 0 && throw(DomainError())
    g = ZZ()
    z = ZZ()
    if y == 1
        return ZZ(0), ZZ(0)
    end
    ccall((:fmpz_gcdinv, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &g, &z, &x, &y)
    return g, z
end

function sqrtmod(x::fmpz, y::fmpz)
    y <= 0 && throw(DomainError())
    z = ZZ()
    if (ccall((:fmpz_sqrtmod, :libflint), Cint, 
              (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &z, &x, &y) == 0)
        error("no square root exists")
    end
    return z
end

function crt(r1::fmpz, m1::fmpz, r2::fmpz, m2::fmpz, signed=false)
   z = ZZ()
   ccall((:fmpz_CRT, :libflint), Void,
          (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Cint),
          &z, &r1, &m1, &r2, &m2, signed)
   return z
end

function crt(r1::fmpz, m1::fmpz, r2::Int, m2::Int, signed=false)
   z = ZZ()
   r2 < 0 && throw(DomainError())
   m2 < 0 && throw(DomainError())
   ccall((:fmpz_CRT_ui, :libflint), Void,
          (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Int, Int, Cint),
          &z, &r1, &m1, r2, m2, signed)
   return z
end

###############################################################################
#
#   Extended GCD
#
###############################################################################

function gcdx(a::fmpz, b::fmpz)
    if b == 0 # shortcut this to ensure consistent results with gcdx(a,b)
        return a < 0 ? (-a, -one(ZZ), zero(ZZ)) : (a, one(ZZ), zero(ZZ))
    end
    g = ZZ()
    s = ZZ()
    t = ZZ()
    ccall((:fmpz_xgcd, :libflint), Void,
        (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}),
        &g, &s, &t, &a, &b)
    g, s, t
end

function gcdinv(a::fmpz, b::fmpz)
   a < 0 && throw(DomainError())
   b < a && throw(DomainError())
   g = ZZ()
   s = ZZ()
   ccall((:fmpz_gcdinv, :libflint), Void,
        (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}),
        &g, &s, &a, &b)
   return g, s
end

###############################################################################
#
#   Comparison
#
###############################################################################

function cmp(x::fmpz, y::fmpz)
    int(ccall((:fmpz_cmp, :libflint), Cint, 
              (Ptr{fmpz}, Ptr{fmpz}), &x, &y))
end

==(x::fmpz, y::fmpz) = cmp(x,y) == 0

<=(x::fmpz, y::fmpz) = cmp(x,y) <= 0

>=(x::fmpz, y::fmpz) = cmp(x,y) >= 0

<(x::fmpz, y::fmpz) = cmp(x,y) < 0

>(x::fmpz, y::fmpz) = cmp(x,y) > 0

function cmpabs(x::fmpz, y::fmpz)
    int(ccall((:fmpz_cmpabs, :libflint), Cint, 
              (Ptr{fmpz}, Ptr{fmpz}), &x, &y))
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

function cmp(x::fmpz, y::Int)
    int(ccall((:fmpz_cmp_si, :libflint), Cint, (Ptr{fmpz}, Int), &x, y))
end

==(x::fmpz, y::Int) = cmp(x,y) == 0

<=(x::fmpz, y::Int) = cmp(x,y) <= 0

>=(x::fmpz, y::Int) = cmp(x,y) >= 0

<(x::fmpz, y::Int) = cmp(x,y) < 0

>(x::fmpz, y::Int) = cmp(x,y) > 0

==(x::Int, y::fmpz) = cmp(y,x) == 0

<=(x::Int, y::fmpz) = cmp(y,x) >= 0

>=(x::Int, y::fmpz) = cmp(y,x) <= 0

<(x::Int, y::fmpz) = cmp(y,x) > 0

>(x::Int, y::fmpz) = cmp(y,x) < 0

###############################################################################
#
#   Bit fiddling
#
###############################################################################

popcount(x::fmpz) = int(ccall((:fmpz_popcnt, :libflint), Culong, 
                              (Ptr{fmpz},), &x))

prevpow2(x::fmpz) = x < 0 ? -prevpow2(-x) :
                            (x <= 2 ? x : one(ZZ) << (ndigits(x, 2) - 1))

nextpow2(x::fmpz) = x < 0 ? -nextpow2(-x) : 
                            (x <= 2 ? x : one(ZZ) << ndigits(x - 1, 2))

trailing_zeros(x::fmpz) = ccall((:fmpz_val2, :libflint), Int, 
                                (Ptr{fmpz},), &x)

###############################################################################
#
#   Bitwise operations (unsafe)
#
###############################################################################

function clrbit!(x::fmpz, c::Int)
    c < 0 && throw(DomainError())
    ccall((:fmpz_clrbit, :libflint), Void, (Ptr{fmpz}, Int), &x, c)
end

function setbit!(x::fmpz, c::Int)
    c < 0 && throw(DomainError())
    ccall((:fmpz_setbit, :libflint), Void, (Ptr{fmpz}, Int), &x, c)
end

function combit!(x::fmpz, c::Int)
    c < 0 && throw(DomainError())
    ccall((:fmpz_combit, :libflint), Void, (Ptr{fmpz}, Int), &x, c)
end

###############################################################################
#
#   Roots
#
###############################################################################

function isqrt(x::fmpz)
    x < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_sqrt, :libflint), Void, (Ptr{fmpz}, Ptr{fmpz}), &z, &x)
    return z
end

function isqrtrem(x::fmpz)
    x < 0 && throw(DomainError())
    s = ZZ()
    r = ZZ()
    ccall((:fmpz_sqrtrem, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &s, &r, &x)
    return s, r
end

function root(x::fmpz, y::Int) 
   x < 0 && iseven(y) && throw(DomainError())
   y <= 0 && throw(DomainError())
   z = ZZ()
   ccall((:fmpz_root, :libflint), Void,
         (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, y)
   return z
end

###############################################################################
#
#   Integer logarithm
#
###############################################################################

function flog(x::fmpz, c::Int)
    c <= 0 && throw(DomainError())
    return ccall((:fmpz_flog_ui, :libflint), Int, 
                 (Ptr{fmpz}, Int), &x, c)
end

function clog(x::fmpz, c::Int)
    c <= 0 && throw(DomainError())
    return ccall((:fmpz_clog_ui, :libflint), Int, 
                 (Ptr{fmpz}, Int), &x, c)
end
    
###############################################################################
#
#   Array arithmetic
#
###############################################################################

function sum(arr::AbstractArray{fmpz})
    n = ZZ(0)
    for i in arr
        ccall((:fmpz_add, :libflint), Void,
            (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}),
            &n, &n, &i)
    end
    return n
end

function prod(arr::AbstractArray{fmpz})
    n = ZZ(1)
    for i in arr
        ccall((:fmpz_mul, :libflint), Void,
            (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}),
            &n, &n, &i)
    end
    return n
end

###############################################################################
#
#   Number theoretic/combinatorial
#
###############################################################################

function divisible(x::fmpz, y::fmpz)
   y == 0 && throw(DivideError())
   bool(ccall((:fmpz_divisible, :libflint), Cint, 
              (Ptr{fmpz}, Ptr{fmpz}), &x, &y))
end

function divisible(x::fmpz, y::Int)
   y == 0 && throw(DivideError())
   bool(ccall((:fmpz_divisible_si, :libflint), Cint, 
              (Ptr{fmpz}, Int), &x, y))
end

issquare(x::fmpz) = bool(ccall((:fmpz_is_square, :libflint), Cint, 
                               (Ptr{fmpz},), &x))

# flint's fmpz_is_prime doesn't work yet
isprime(x::fmpz) = bool(ccall((:fmpz_is_probabprime, :libflint), Cint, 
                              (Ptr{fmpz},), &x))

isprobabprime(x::fmpz) = bool(ccall((:fmpz_is_probabprime, :libflint), Cint, 
                                    (Ptr{fmpz},), &x))

function remove(x::fmpz, y::fmpz) 
   y == 0 && throw(DivideError())
   z = ZZ()
   num = ccall((:fmpz_remove, :libflint), Int, 
               (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &z, &x, &y)
   return num, z
end

function divisor_lenstra(n::fmpz, r::fmpz, m::fmpz)
   r <= 0 && throw(DomainError())
   m <= r && throw(DomainError())
   n <= m && throw(DomainError())
   z = ZZ()
   if !bool(ccall((:fmpz_divisor_in_residue_class_lenstra, :libflint), 
       Cint, (Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}, Ptr{fmpz}), &z, &n, &r, &m))
      z = 0
   end
   return z
end

function fac(x::Int)
    x < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_fac_ui, :libflint), Void, (Ptr{fmpz}, Culong), &z, x)
    return z
end

function risingfac(x::fmpz, y::Int)
    y < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_rfac_ui, :libflint), Void, 
          (Ptr{fmpz}, Ptr{fmpz}, Culong), &z, &x, y)
    return z
end

function risingfac(x::Int, y::Int)
    y < 0 && throw(DomainError())
    z = ZZ()
    if x < 0
       if y <= -x # we don't pass zero
          z = isodd(y) ? -risingfac(-x - y + 1, y) : 
                          risingfac(-x - y + 1, y)
       end
    else
       ccall((:fmpz_rfac_uiui, :libflint), Void, 
             (Ptr{fmpz}, Culong, Culong), &z, x, y)
    end
    return z
end

function primorial(x::Int)
    x < 0 && throw(DomainError()) 
    z = ZZ()
    ccall((:fmpz_primorial, :libflint), Void, 
          (Ptr{fmpz}, Culong), &z, x)
    return z
end

function fib(x::Int)
    x < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_fib_ui, :libflint), Void, 
          (Ptr{fmpz}, Culong), &z, x)
    return z
end

function binom(n::Int, k::Int)
    n < 0 && return ZZ(0)
    k < 0 && return ZZ(0)
    z = ZZ()
    ccall((:fmpz_bin_uiui, :libflint), Void, 
          (Ptr{fmpz}, Culong, Culong), &z, n, k)
    return z
end

function moebiusmu(x::fmpz) 
   x < 0 && throw(DomainError())
   return int(ccall((:fmpz_moebius_mu, :libflint), Cint, 
                    (Ptr{fmpz},), &x))
end

function jacobi(x::fmpz, y::fmpz) 
   y <= x && throw(DomainError())
   x < 0 && throw(DomainError())
   return int(ccall((:fmpz_jacobi, :libflint), Cint, 
                    (Ptr{fmpz}, Ptr{fmpz}), &x, &y))
end

function sigma(x::fmpz, y::Int) 
   y < 0 && throw(DomainError())
   z = ZZ()
   ccall((:fmpz_divisor_sigma, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fmpz}, Int), &z, &x, y)
   return z
end

function eulerphi(x::fmpz) 
   x < 0 && throw(DomainError())
   z = ZZ()
   ccall((:fmpz_euler_phi, :libflint), Void, 
         (Ptr{fmpz}, Ptr{fmpz}), &z, &x)
   return z
end

###############################################################################
#
#   String I/O
#
###############################################################################

string(x::fmpz) = dec(x)

show(io::IO, x::fmpz) = print(io, string(x))

show(io::IO, a::IntegerRing) = print(io, "Integer Ring")

needs_parentheses(x::fmpz) = false

is_negative(x::fmpz) = x < 0

show_minus_one(::Type{fmpz}) = false

###############################################################################
#
#   Number bases/digits
#
###############################################################################

bin(n::fmpz) = base(n, 2)

oct(n::fmpz) = base(n, 8)

dec(n::fmpz) = base(n, 10)

hex(n::fmpz) = base(n, 16)

function base(n::fmpz, b::Integer)
    2 <= b <= 62 || error("invalid base: $b")
    p = ccall((:fmpz_get_str,:libflint), Ptr{Uint8}, 
              (Ptr{Uint8}, Cint, Ptr{fmpz}), C_NULL, b, &n)
    len = int(ccall(:strlen, Csize_t, (Ptr{Uint8},), p))
    ASCIIString(pointer_to_array(p, len, true))
end

function ndigits_internal(x::fmpz, b::Integer = 10)
    # fmpz_sizeinbase might return an answer 1 too big
    n = int(ccall((:fmpz_sizeinbase, :libflint), Culong, 
                  (Ptr{fmpz}, Int32), &x, b))
    abs(x) < ZZ(b)^(n - 1) ? n - 1 : n
end

ndigits(x::fmpz, b::Integer = 10) = x == 0 ? 1 : ndigits_internal(x, b)

nbits(x::fmpz) = x == 0 ? 0 : ndigits(x, 2)

###############################################################################
#
#   Parent object overloads
#
###############################################################################

call(::IntegerRing) = fmpz()

call(::IntegerRing, a::Integer) = fmpz(a)

call(::IntegerRing, a::String) = fmpz(a)

call(::IntegerRing, a::fmpz) = a

###############################################################################
#
#   String parser
#
###############################################################################

function parseint(::Type{fmpz}, s::String, base::Int = 10)
    s = bytestring(s)
    sgn = s[1] == '-' ? -1 : 1
    i = 1 + (sgn == -1)
    z = ZZ()
    err = ccall((:fmpz_set_str, :libflint),
               Int32, (Ptr{fmpz}, Ptr{Uint8}, Int32),
               &z, convert(Ptr{Uint8},SubString(s,i)), base)
    err == 0 || error("Invalid big integer: $(repr(s))")
    return sgn < 0 ? -z : z
end

###############################################################################
#
#   Constructors
#
###############################################################################

fmpz(s::String) = parseint(fmpz, s)

fmpz(z::Integer) = ZZ(BigInt(z))

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

convert(::Type{fmpz}, a::Integer) = ZZ(a)

Base.promote_rule{T <: Integer}(::Type{fmpz}, ::Type{T}) = fmpz

function BigInt(z::fmpz)
   r = BigInt()
   ccall((:fmpz_get_mpz, :libflint), Void, (Ptr{BigInt}, Ptr{fmpz}), &r, &z)
   return r
end

convert(::Type{BigInt}, a::fmpz) = BigInt(a)
