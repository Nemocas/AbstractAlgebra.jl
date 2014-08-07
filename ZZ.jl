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


import Base: convert, promote_rule, show, string, parseint, serialize, deserialize,
             base, bin, dec, oct, hex, gcd, div

export ZZ, fac, binom, isprime, fdiv, cdiv, tdiv, div, rem, mod, gcd, xgcd, lcm, invmod, 
       powmod, abs, divrem, isqrt, popcount, prevpow2, nextpow2, ndigits, dec, bin, oct, 
       hex, base, show, convert, one, zero, divexact, fits, sign, nbits, deepcopy,
       tdivpow2, fdivpow2, cdivpow2, flog, clog, cmpabs, clrbit!, setbit!, combit!,
       crt, divisible, divisor_lenstra, fdivrem, tdivrem, fmodpow2, gcdinv, isprobabprime,
       issquare, jacobi, remove, root, size, isqrtrem, sqrtmod, trailing_zeros, sigma,
       eulerphi, fib, moebiusmu, primorial, risingfac

###########################################################################################
#
#   Helper macros/aliases
#
###########################################################################################

# We define unions of integer types to save defining functions for each type in the union

if Clong == Int32
    typealias ClongMax Union(Int8, Int16, Int32)
    typealias CulongMax Union(Uint8, Uint16, Uint32)
else
    typealias ClongMax Union(Int8, Int16, Int32, Int64)
    typealias CulongMax Union(Uint8, Uint16, Uint32, Uint64)
end

###########################################################################################
#
#   Data type and memory management
#
###########################################################################################

# The ZZ datatype. It is a member of the Ring abstract type

type ZZ <: Ring
    d::Int

    function ZZ()
        b = new(zero(Int))
        ccall((:__fmpz_init, :libflint), Void, (Ptr{ZZ},), &b)
        finalizer(b, _fmpz_clear_fn)
        return b
    end
end

function _fmpz_clear_fn(a::ZZ)
   ccall((:__fmpz_clear, :libflint), Void, (Ptr{ZZ},), &a)
end

###########################################################################################
#
#   Constructors
#
###########################################################################################

ZZ(x::ZZ) = x

ZZ(s::String) = parseint(ZZ, s)

# Turn a string into a ZZ

function Base.parseint(::Type{ZZ}, s::String, base::Int)
    s = bytestring(s)
    sgn = s[1] == '-' ? -1 : 1
    i = 1 + (sgn == -1)
    z = ZZ()
    err = ccall((:fmpz_set_str, :libflint),
               Int32, (Ptr{ZZ}, Ptr{Uint8}, Int32),
               &z, convert(Ptr{Uint8},SubString(s,i)), base)
    err == 0 || error("invalid big integer: $(repr(s))")
    return sgn < 0 ? -z : z
end

function ZZ(x::Union(Clong, Int))
    z = ZZ()
    ccall((:__fmpz_set_si, :libflint), Void, (Ptr{ZZ}, Int), &z, x)
    return z
end

function ZZ(x::Union(Culong, Uint))
    z = ZZ()
    ccall((:__fmpz_set_ui, :libflint), Void, (Ptr{ZZ}, Uint), &z, x)
    return z
end

ZZ(x::Bool) = ZZ(uint(x))

function ZZ(x::Float64)
    !isinteger(x) && throw(InexactError())
    z = ZZ()
    ccall((:fmpz_set_d, :libflint), Void, (Ptr{ZZ}, Cdouble), &z, x)
    return z
end

ZZ(x::Union(Float16, Float32)) = ZZ(float64(x))

# Any other kind of Integer will be taken 32 bits at a time

function ZZ(x::Integer)
    if x < 0
        if typemin(Clong) <= x
            return ZZ(convert(Clong, x))
        end
        b = ZZ(0)
        shift = 0
        while x < -1
            b += ZZ(~uint32(x & 0xffffffff)) << shift
            x >>= 32
            shift += 32
        end
        return -b - 1
    else
        if x <= typemax(Culong)
            return ZZ(convert(Culong, x))
        end
        b = ZZ(0)
        shift = 0
        while x > 0
            b += ZZ(uint32(x & 0xffffffff)) << shift
            x >>>= 32
            shift += 32
        end
        return b
    end
end

function deepcopy(a::ZZ)
   z = ZZ()
   ccall((:fmpz_set, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}), &z, &a)
   return z
end
 
###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

one(::Type{ZZ}) = ZZ(1)

Base.zero(::Type{ZZ}) = ZZ(0)

sign(a::ZZ) = int(ccall((:fmpz_sgn, :libflint), Cint, (Ptr{ZZ},), &a))

fits(::Type{Int}, a::ZZ) = ccall((:fmpz_fits_si, :libflint), Bool, (Ptr{ZZ},), &a)

fits(::Type{Uint}, a::ZZ) = sign(a) < 0 ? false : ccall((:fmpz_abs_fits_ui, :libflint), Bool, (Ptr{ZZ},), &a)

size(a::ZZ) = int(ccall((:fmpz_size, :libflint), Cint, (Ptr{ZZ},), &a))

###########################################################################################
#
#   Conversions to/from other Julia types
#
###########################################################################################

convert(::Type{ZZ}, x::Integer) = ZZ(x)

convert(::Type{ZZ}, x::Float16) = ZZ(x)

convert(::Type{ZZ}, x::FloatingPoint) = ZZ(x)

function convert(::Type{Int64}, x::ZZ)
    lo = int64(convert(Culong, x & typemax(Uint32)))
    hi = int64(convert(Clong, x >> 32))
    hi << 32 | lo
end

convert(::Type{Int32}, n::ZZ) = int32(convert(Clong, n))

convert(::Type{Int16}, n::ZZ) = int16(convert(Clong, n))

convert(::Type{Int8}, n::ZZ) = int8(convert(Clong, n))

function convert(::Type{Clong}, n::ZZ)
    fits = ccall((:fmpz_fits_si, :libflint), Int32, (Ptr{ZZ},), &n) != 0
    if fits
        ccall((:fmpz_get_si, :libflint), Clong, (Ptr{ZZ},), &n)
    else
        throw(InexactError())
    end
end

function convert(::Type{Uint64}, x::ZZ)
    lo = uint64(convert(Culong, x & typemax(Uint32)))
    hi = uint64(convert(Culong, x >> 32))
    hi << 32 | lo
end

convert(::Type{Uint32}, x::ZZ) = uint32(convert(Culong, x))

convert(::Type{Uint16}, x::ZZ) = uint16(convert(Culong, x))

convert(::Type{Uint8}, x::ZZ) = uint8(convert(Culong, x))

function convert(::Type{Culong}, n::ZZ)
    fits = ccall((:fmpz_fits_si, :libflint), Int32, (Ptr{ZZ},), &n) != 0
    if fits
        ccall((:fmpz_get_ui, :libflint), Culong, (Ptr{ZZ},), &n)
    else
        throw(InexactError())
    end
end

if sizeof(Int32) == sizeof(Clong)
    function convert(::Type{Uint128}, x::ZZ)
        uint128(uint(x>>>96))<<96 +
        uint128(uint((x>>>64) & typemax(Uint32)))<<64 +
        uint128(uint((x>>>32) & typemax(Uint32)))<<32 +
        uint128(uint(x & typemax(Uint32)))
    end
end
if sizeof(Int64) == sizeof(Clong)
    function convert(::Type{Uint128}, x::ZZ)
        uint128(uint(x>>>64))<<64 +
        uint128(uint(x & typemax(Uint64)))
    end
end

convert(::Type{Int128}, x::ZZ) = copysign(int128(uint128(abs(x))),x)

function convert(::Type{Float64}, n::ZZ)
    # rounds to zero
    ccall((:fmpz_get_d, :libflint), Float64, (Ptr{ZZ},), &n)
end

convert(::Type{Float32}, n::ZZ) = float32(float64(n))

convert(::Type{Float16}, n::ZZ) = float16(float64(n))

###########################################################################################
#
#   Automatic promotion rules for operators
#
###########################################################################################

promote_rule{T<:Integer}(::Type{ZZ}, ::Type{T}) = ZZ

promote_rule{T<:Real}(::Type{ZZ}, ::Type{T}) = ZZ

###########################################################################################
#
#   Serialisation
#
###########################################################################################

function serialize(s, n::ZZ)
    Base.serialize_type(s, ZZ)
    serialize(s, base(62,n))
end

deserialize(s, ::Type{ZZ}) = Base.parseint_nocheck(ZZ, deserialize(s), 62)

###########################################################################################
#
#   Binary operators and functions
#
###########################################################################################

# Metaprogram to define functions +, -, *, gcd, lcm, 
#                                 &, |, $

for (fJ, fC) in ((:+, :add), (:-,:sub), (:*, :mul),
                 (:gcd, :gcd), (:lcm, :lcm),
                 (:&, :and), (:|, :or), (:$, :xor))
    @eval begin
        function ($fJ)(x::ZZ, y::ZZ)
            z = ZZ()
            ccall(($(string(:fmpz_, fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &x, &y)
            return z
        end
    end
end

# Metaprogram to define functionsfdiv, cdiv, tdiv, div, mod

for (fJ, fC) in ((:fdiv, :fdiv_q), (:cdiv, :cdiv_q), (:tdiv, :tdiv_q), 
                 (:div, :tdiv_q), (:mod, :mod))
    @eval begin
        function ($fJ)(x::ZZ, y::ZZ)
            y == 0 && throw(DivideError())
            z = ZZ()
            ccall(($(string(:fmpz_, fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &x, &y)
            return z
        end
    end
end

function divexact(x::ZZ, y::ZZ)
    y == 0 && throw(DivideError())
    z = ZZ()
    ccall((:fmpz_divexact, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &x, &y)
    z
end

function flog(x::ZZ, c::ZZ)
    c <= 0 && throw(DomainError())
    x <= 0 && throw(DomainError())
    return ccall((:fmpz_flog, :libflint), Int, (Ptr{ZZ}, Ptr{ZZ}), &x, &c)
end

function clog(x::ZZ, c::ZZ)
    c <= 0 && throw(DomainError())
    x <= 0 && throw(DomainError())
    return ccall((:fmpz_clog, :libflint), Int, (Ptr{ZZ}, Ptr{ZZ}), &x, &c)
end

function %(x::ZZ, c::ZZ)
    c == 0 && throw(DivideError())
    q = ZZ()
    r = ZZ()
    ccall((:fmpz_tdiv_qr, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &q, &r, &x, &c)
    return r
end

rem(x::ZZ, c::ZZ) = %(x, c)

###########################################################################################
#
#   Multi operators (e.g. d = +(a, b, c) instead of d = a + b + c)
#
###########################################################################################

for (fJ, fC) in ((:+, :add), (:*, :mul), (:&, :and), (:|, :or), (:$, :xor))
    @eval begin
        function ($fJ)(a::ZZ, b::ZZ, c::ZZ)
            z = ZZ()
            ccall(($(string(:fmpz_, fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &a, &b)
            ccall(($(string(:fmpz_, fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &z, &c)
            return z
        end
        function ($fJ)(a::ZZ, b::ZZ, c::ZZ, d::ZZ)
            z = ZZ()
            ccall(($(string(:fmpz_, fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &a, &b)
            ccall(($(string(:fmpz_, fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &z, &c)
            ccall(($(string(:fmpz_, fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &z, &d)
            return z
        end
        function ($fJ)(a::ZZ, b::ZZ, c::ZZ, d::ZZ, e::ZZ)
            z = ZZ()
            ccall(($(string(:fmpz_, fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &a, &b)
            ccall(($(string(:fmpz_, fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &z, &c)
            ccall(($(string(:fmpz_, fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &z, &d)
            ccall(($(string(:fmpz_, fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &z, &e)
            return z
        end
    end
end

###########################################################################################
#
#   Binary operators and functions where one operand is base type (ad hoc polymorphism)
#
###########################################################################################

function +(x::ZZ, c::Uint)
    z = ZZ()
    ccall((:fmpz_add_ui, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Uint), &z, &x, c)
    return z
end

+(c::Uint, x::ZZ) = x + c

+(c::CulongMax, x::ZZ) = x + convert(Uint, c)

+(x::ZZ, c::CulongMax) = x + convert(Uint, c)

+(x::ZZ, c::ClongMax) = c < 0 ? -(x, convert(Uint, -c)) : x + convert(Uint, c)

+(c::ClongMax, x::ZZ) = c < 0 ? -(x, convert(Uint, -c)) : x + convert(Uint, c)

function -(x::ZZ, c::Uint)
    z = ZZ()
    ccall((:fmpz_sub_ui, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Uint), &z, &x, c)
    return z
end

function -(c::Uint, x::ZZ)
    z = ZZ()
    ccall((:fmpz_sub_ui, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Uint), &z, &x, c)
    ccall((:__fmpz_neg, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}), &z, &z)
    return z
end

-(x::ZZ, c::CulongMax) = -(x, convert(Uint, c))

-(c::CulongMax, x::ZZ) = -(convert(Uint, c), x)

-(x::ZZ, c::ClongMax) = c < 0 ? +(x, convert(Uint, -c)) : -(x, convert(Uint, c))

-(c::ClongMax, x::ZZ) = c < 0 ? -(x + convert(Uint, -c)) : -(convert(Uint, c), x)

function *(x::ZZ, c::Uint)
    z = ZZ()
    ccall((:fmpz_mul_ui, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Uint), &z, &x, c)
    return z
end

*(c::Uint, x::ZZ) = x * c

*(c::CulongMax, x::ZZ) = x * convert(Uint, c)

*(x::ZZ, c::CulongMax) = x * convert(Uint, c)

function *(x::ZZ, c::Int)
    z = ZZ()
    ccall((:fmpz_mul_si, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, c)
    return z
end

*(c::Int, x::ZZ) = x * c

*(x::ZZ, c::ClongMax) = x * convert(Int, c)

*(c::ClongMax, x::ZZ) = x * convert(Int, c)

function <<(x::ZZ, c::Int)
    c < 0 && throw(DomainError())
    c == 0 && return x
    z = ZZ()
    ccall((:fmpz_mul_2exp, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, c)
    return z
end

function >>(x::ZZ, c::Int)
    c < 0 && throw(DomainError())
    c == 0 && return x
    z = ZZ()
    ccall((:fmpz_fdiv_q_2exp, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, c)
    return z
end

function %(x::ZZ, c::Int)
   c < 0 && throw(DomainError())
   c == 0 && throw(DivideError())
   r = ccall((:fmpz_tdiv_ui, :libflint), Int, (Ptr{ZZ}, Int), &x, c)
   return sign(x) < 0 ? -r : r
end

rem(x::ZZ, c::Int) = %(X, C)

function tdivpow2(x::ZZ, c::Int)
    c < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_tdiv_q_2exp, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, c)
    return z
end

function fdivpow2(x::ZZ, c::Int)
    c < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_fdiv_q_2exp, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, c)
    return z
end

function fmodpow2(x::ZZ, c::Int)
    c < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_fdiv_r_2exp, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, c)
    return z
end

function cdivpow2(x::ZZ, c::Int)
    c < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_cdiv_q_2exp, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, c)
    return z
end

function div(x::ZZ, c::Int)
    c == 0 && throw(DivideError())
    z = ZZ()
    ccall((:fmpz_tdiv_q_si, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, c)
    return z
end

function tdiv(x::ZZ, c::Int)
    c == 0 && throw(DivideError())
    z = ZZ()
    ccall((:fmpz_tdiv_q_si, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, c)
    return z
end

function fdiv(x::ZZ, c::Int)
    c == 0 && throw(DivideError())
    z = ZZ()
    ccall((:fmpz_fdiv_q_si, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, c)
    return z
end

function cdiv(x::ZZ, c::Int)
    c == 0 && throw(DivideError())
    z = ZZ()
    ccall((:fmpz_cdiv_si, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, c)
    return z
end

function mod(x::ZZ, c::Int)
    c < 0 && throw(DomainError())
    c == 0 && throw(DivideError())
    ccall((:fmpz_fdiv_ui, :libflint), Int, (Ptr{ZZ}, Int), &x, c)
end

function divexact(x::ZZ, y::Int)
    y == 0 && throw(DivideError())
    z = ZZ()
    ccall((:fmpz_divexact_si, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, y)
    z
end

function flog(x::ZZ, c::Int)
    c <= 0 && throw(DomainError())
    return ccall((:fmpz_flog_ui, :libflint), Int, (Ptr{ZZ}, Int), &x, c)
end

function clog(x::ZZ, c::Int)
    c <= 0 && throw(DomainError())
    return ccall((:fmpz_clog_ui, :libflint), Int, (Ptr{ZZ}, Int), &x, c)
end

function ^(x::ZZ, y::Uint)
    z = ZZ()
    ccall((:fmpz_pow_ui, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Culong), &z, &x, y)
    return z
end

function zz_pow(x::ZZ, y::Int)
    if y < 0; throw(DomainError()); end
    if x == 1; return x; end
    if x == -1; return isodd(y) ? x : -x; end
    if y > typemax(Uint); throw(DomainError()); end
    return x^uint(y)
end

^(x::ZZ, y::Bool) = y ? x : one(x)

^(x::ZZ, y::Int) = zz_pow(x, y)

###########################################################################################
#
#   Unary operators and functions, e.g. -ZZ(12), ~ZZ(12)
#
###########################################################################################

function -(x::ZZ)
    z = ZZ()
    ccall((:__fmpz_neg, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}), &z, &x)
    return z
end

function ~(x::ZZ)
    z = ZZ()
    ccall((:fmpz_complement, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}), &z, &x)
    return z
end

function abs(x::ZZ)
    z = ZZ()
    ccall((:fmpz_abs, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}), &z, &x)
    return z
end

###########################################################################################
#
#   Division with remainder
#
###########################################################################################

function divrem(x::ZZ, y::ZZ)
    y == 0 && throw(DivideError())
    z1 = ZZ()
    z2 = ZZ()
    ccall((:fmpz_tdiv_qr, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z1, &z2, &x, &y)
    z1, z2
end

function tdivrem(x::ZZ, y::ZZ)
    y == 0 && throw(DivideError())
    z1 = ZZ()
    z2 = ZZ()
    ccall((:fmpz_tdiv_qr, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z1, &z2, &x, &y)
    z1, z2
end

function fdivrem(x::ZZ, y::ZZ)
    y == 0 && throw(DivideError())
    z1 = ZZ()
    z2 = ZZ()
    ccall((:fmpz_fdiv_qr, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z1, &z2, &x, &y)
    z1, z2
end

###########################################################################################
#
#   Roots
#
###########################################################################################

function isqrt(x::ZZ)
    x < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_sqrt, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}), &z, &x)
    return z
end

function isqrtrem(x::ZZ)
    x < 0 && throw(DomainError())
    s = ZZ()
    r = ZZ()
    ccall((:fmpz_sqrtrem, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &s, &r, &x)
    return s, r
end

function root(x::ZZ, y::Int) 
   x < 0 && iseven(y) && throw(DomainError())
   y <= 0 && throw(DomainError())
   z = ZZ()
   ccall((:fmpz_root, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, y)
   return z
end

###########################################################################################
#
#   Extended GCD
#
###########################################################################################

function xgcd(a::ZZ, b::ZZ)
    if b == 0 # shortcut this to ensure consistent results with gcdx(a,b)
        return a < 0 ? (-a, -one(ZZ), zero(ZZ)) : (a, one(ZZ), zero(ZZ))
    end
    g = ZZ()
    s = ZZ()
    t = ZZ()
    ccall((:fmpz_xgcd, :libflint), Void,
        (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}),
        &g, &s, &t, &a, &b)
    g, s, t
end

function gcdinv(a::ZZ, b::ZZ)
   a < 0 && throw(DomainError())
   b < a && throw(DomainError())
   g = ZZ()
   s = ZZ()
   ccall((:fmpz_gcdinv, :libflint), Void,
        (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}),
        &g, &s, &a, &b)
   return g, s
end

###########################################################################################
#
#   Comparison
#
###########################################################################################

function cmp(x::ZZ, y::ZZ)
    int(ccall((:fmpz_cmp, :libflint), Cint, (Ptr{ZZ}, Ptr{ZZ}), &x, &y))
end

==(x::ZZ, y::ZZ) = cmp(x,y) == 0

<=(x::ZZ, y::ZZ) = cmp(x,y) <= 0

>=(x::ZZ, y::ZZ) = cmp(x,y) >= 0

<(x::ZZ, y::ZZ) = cmp(x,y) < 0

>(x::ZZ, y::ZZ) = cmp(x,y) > 0

function cmp(x::ZZ, y::Int)
    int(ccall((:fmpz_cmp_si, :libflint), Cint, (Ptr{ZZ}, Int), &x, y))
end

==(x::ZZ, y::Int) = cmp(x,y) == 0

<=(x::ZZ, y::Int) = cmp(x,y) <= 0

>=(x::ZZ, y::Int) = cmp(x,y) >= 0

<(x::ZZ, y::Int) = cmp(x,y) < 0

>(x::ZZ, y::Int) = cmp(x,y) > 0

==(x::Int, y::ZZ) = cmp(y,x) == 0

<=(x::Int, y::ZZ) = cmp(y,x) >= 0

>=(x::Int, y::ZZ) = cmp(y,x) <= 0

<(x::Int, y::ZZ) = cmp(y,x) > 0

>(x::Int, y::ZZ) = cmp(y,x) < 0

function cmpabs(x::ZZ, y::ZZ)
    int(ccall((:fmpz_cmpabs, :libflint), Cint, (Ptr{ZZ}, Ptr{ZZ}), &x, &y))
end


###########################################################################################
#
#   Bit fiddling
#
###########################################################################################

popcount(x::ZZ) = int(ccall((:fmpz_popcnt, :libflint), Culong, (Ptr{ZZ},), &x))

prevpow2(x::ZZ) = x < 0 ? -prevpow2(-x) : (x <= 2 ? x : one(ZZ) << (ndigits(x, 2) - 1))

nextpow2(x::ZZ) = x < 0 ? -nextpow2(-x) : (x <= 2 ? x : one(ZZ) << ndigits(x - 1, 2))

trailing_zeros(x::ZZ) = ccall((:fmpz_val2, :libflint), Int, (Ptr{ZZ},), &x)

###########################################################################################
#
#   Bitwise operations (unsafe)
#
###########################################################################################

function clrbit!(x::ZZ, c::Int)
    c < 0 && throw(DomainError())
    ccall((:fmpz_clrbit, :libflint), Void, (Ptr{ZZ}, Int), &x, c)
end

function setbit!(x::ZZ, c::Int)
    c < 0 && throw(DomainError())
    ccall((:fmpz_setbit, :libflint), Void, (Ptr{ZZ}, Int), &x, c)
end

function combit!(x::ZZ, c::Int)
    c < 0 && throw(DomainError())
    ccall((:fmpz_combit, :libflint), Void, (Ptr{ZZ}, Int), &x, c)
end

###########################################################################################
#
#   Modular arithmetic
#
###########################################################################################

function powmod(x::ZZ, p::ZZ, m::ZZ)
    m <= 0 && throw(DomainError())
    if p < 0
       x = invmod(x, m)
       p = -p
    end
    r = ZZ()
    ccall((:fmpz_powm, :libflint), Void,
          (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}),
          &r, &x, &p, &m)
    return r 
end

function powmod(x::ZZ, p::Int, m::ZZ)
    m <= 0 && throw(DomainError())
    if p < 0
       x = invmod(x, m)
       p = -p
    end
    r = ZZ()
    ccall((:fmpz_powm_ui, :libflint), Void,
          (Ptr{ZZ}, Ptr{ZZ}, Int, Ptr{ZZ}),
          &r, &x, p, &m)
    return r 
end

function invmod(x::ZZ, y::ZZ)
    y <= 0 && throw(DomainError())
    z = ZZ()
    if y == 1
        return ZZ(0)
    end
    if ccall((:fmpz_invmod, :libflint), Cint, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &x, &y) == 0
       error("Impossible inverse in invmod")
    end
    return z
end

function gcdinv(x::ZZ, y::ZZ)
    y <= 0 && throw(DomainError())
    g = ZZ()
    z = ZZ()
    if y == 1
        return ZZ(0), ZZ(0)
    end
    ccall((:fmpz_gcdinv, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &g, &z, &x, &y)
    return g, z
end

function sqrtmod(x::ZZ, y::ZZ)
    y <= 0 && throw(DomainError())
    z = ZZ()
    if (ccall((:fmpz_sqrtmod, :libflint), Cint, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &x, &y) == 0)
        error("no square root exists")
    end
    return z
end

function crt(r1::ZZ, m1::ZZ, r2::ZZ, m2::ZZ, sign::Bool)
   z = ZZ()
   ccall((:fmpz_CRT, :libflint), Void,
          (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Cint),
          &z, &r1, &m1, &r2, &m2, sign)
   return z
end

function crt(r1::ZZ, m1::ZZ, r2::Int, m2::Int, sign::Bool)
   z = ZZ()
   r2 < 0 && throw(DomainError())
   m2 < 0 && throw(DomainError())
   ccall((:fmpz_CRT_ui, :libflint), Void,
          (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Int, Int, Cint),
          &z, &r1, &m1, r2, m2, sign)
   return z
end
    
###########################################################################################
#
#   Array arithmetic
#
###########################################################################################

function sum(arr::AbstractArray{ZZ})
    n = ZZ(0)
    for i in arr
        ccall((:fmpz_add, :libflint), Void,
            (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}),
            &n, &n, &i)
    end
    return n
end

function prod(arr::AbstractArray{ZZ})
    n = ZZ(1)
    for i in arr
        ccall((:fmpz_mul, :libflint), Void,
            (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}),
            &n, &n, &i)
    end
    return n
end

###########################################################################################
#
#   Number theoretic/combinatorial
#
###########################################################################################

function divisible(x::ZZ, y::ZZ)
   y == 0 && throw(DivideError())
   bool(ccall((:fmpz_divisible, :libflint), Cint, (Ptr{ZZ}, Ptr{ZZ}), &x, &y))
end

function divisible(x::ZZ, y::Int)
   y == 0 && throw(DivideError())
   bool(ccall((:fmpz_divisible_si, :libflint), Cint, (Ptr{ZZ}, Int), &x, y))
end

issquare(x::ZZ) = bool(ccall((:fmpz_is_square, :libflint), Cint, (Ptr{ZZ},), &x))

isprime(x::ZZ) = bool(ccall((:fmpz_is_prime, :libflint), Cint, (Ptr{ZZ},), &x))

isprobabprime(x::ZZ) = bool(ccall((:fmpz_is_probabprime, :libflint), Cint, (Ptr{ZZ},), &x))

function remove(x::ZZ, y::ZZ) 
   y == 0 && throw(DivideError())
   z = ZZ()
   num = ccall((:fmpz_remove, :libflint), Int, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &x, &y)
   return num, z
end

function divisor_lenstra(n::ZZ, r::ZZ, m::ZZ)
   r <= 0 && throw(DomainError())
   m <= r && throw(DomainError())
   n <= m && throw(DomainError())
   z = ZZ()
   if !bool(ccall((:fmpz_divisor_in_residue_class_lenstra, :libflint), 
       Cint, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &n, &r, &m))
      z = 0
   end
   return z
end

function fac(x::Int)
    x < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_fac_ui, :libflint), Void, (Ptr{ZZ}, Culong), &z, x)
    return z
end

function risingfac(x::ZZ, y::Int)
    y < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_rfac_ui, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Culong), &z, &x, y)
    return z
end

function risingfac(x::Int, y::Int)
    y < 0 && throw(DomainError())
    z = ZZ()
    if x < 0
       if y <= -x # we don't pass zero
          z = isodd(y) ? -risingfac(-x - y + 1, y) : risingfac(-x - y + 1, y)
       end
    else
       ccall((:fmpz_rfac_uiui, :libflint), Void, (Ptr{ZZ}, Culong, Culong), &z, x, y)
    end
    return z
end

function primorial(x::Int)
    x < 0 && throw(DomainError()) 
    z = ZZ()
    ccall((:fmpz_primorial, :libflint), Void, (Ptr{ZZ}, Culong), &z, x)
    return z
end

function fib(x::Int)
    x < 0 && throw(DomainError())
    z = ZZ()
    ccall((:fmpz_fib_ui, :libflint), Void, (Ptr{ZZ}, Culong), &z, x)
    return z
end

function binom(n::Int, k::Int)
    n < 0 && return ZZ(0)
    k < 0 && return ZZ(0)
    z = ZZ()
    ccall((:fmpz_bin_uiui, :libflint), Void, (Ptr{ZZ}, Culong, Culong), &z, n, k)
    return z
end

function moebiusmu(x::ZZ) 
   x < 0 && throw(DomainError())
   return int(ccall((:fmpz_moebius_mu, :libflint), Cint, (Ptr{ZZ},), &x))
end

function jacobi(x::ZZ, y::ZZ) 
   y <= x && throw(DomainError())
   x < 0 && throw(DomainError())
   return int(ccall((:fmpz_jacobi, :libflint), Cint, (Ptr{ZZ}, Ptr{ZZ}), &x, &y))
end

function sigma(x::ZZ, y::Int) 
   y < 0 && throw(DomainError())
   z = ZZ()
   ccall((:fmpz_divisor_sigma, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Int), &z, &x, y)
   return z
end

function eulerphi(x::ZZ) 
   x < 0 && throw(DomainError())
   z = ZZ()
   ccall((:fmpz_euler_phi, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}), &z, &x)
   return z
end

###########################################################################################
#
#   String I/O
#
###########################################################################################

string(x::ZZ) = dec(x)

show(io::IO, x::ZZ) = print(io, string(x))

###########################################################################################
#
#   Number bases/digits
#
###########################################################################################

bin(n::ZZ) = base(n, 2)

oct(n::ZZ) = base(n, 8)

dec(n::ZZ) = base(n, 10)

hex(n::ZZ) = base(n, 16)

function base(n::ZZ, b::Integer)
    2 <= b <= 62 || error("invalid base: $b")
    p = ccall((:fmpz_get_str,:libflint), Ptr{Uint8}, (Ptr{Uint8}, Cint, Ptr{ZZ}), C_NULL, b, &n)
    len = int(ccall(:strlen, Csize_t, (Ptr{Uint8},), p))
    ASCIIString(pointer_to_array(p,len,true))
end

function ndigits_internal(x::ZZ, b::Integer = 10)
    # fmpz_sizeinbase might return an answer 1 too big
    n = int(ccall((:fmpz_sizeinbase, :libflint), Culong, (Ptr{ZZ}, Int32), &x, b))
    abs(x) < ZZ(b)^(n - 1) ? n - 1 : n
end

ndigits(x::ZZ, b::Integer = 10) = x == 0 ? 1 : ndigits_internal(x, b)

nbits(x::ZZ) = x == 0 ? 0 : ndigits(x, 2)
