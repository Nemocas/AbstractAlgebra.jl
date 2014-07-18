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


import Base: convert, promote_rule, show, string, parseint, serialize, deserialize

export fac, binom, isprime, fdiv, div, rem, mod, gcd, xgcd, lcm, invmod, powmod, abs, 
       divrem, isqrt, popcount, prevpow2, nextpow2, num_digits, dec, bin, oct, hex, base

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

# The ZZ datatype. It is a member of Julia's Integer abstract type

type ZZ <: Integer
    d::Clong
    function ZZ()
        b = new(zero(Clong))
        ccall((:__fmpz_init,:libflint), Void, (Ptr{ZZ},), &b)
        finalizer(b, _fmpz_clear_fn)
        return b
    end
end

function _fmpz_clear_fn(a::ZZ)
   ccall((:__fmpz_clear,:libflint), Void, (Ptr{ZZ},), &a)
end

###########################################################################################
#
#   Constructors
#
###########################################################################################

ZZ(x::ZZ) = x

ZZ(s::String) = parseint(ZZ, s)

# Turn a string into a ZZ

function parseint_nocheck(::Type{ZZ}, s::String, base::Int)
    s = bytestring(s)
    sgn, base, i = Base.parseint_preamble(true, s, base)
    z = ZZ()
    err = ccall((:fmpz_set_str, :libflint),
               Int32, (Ptr{ZZ}, Ptr{Uint8}, Int32),
               &z, convert(Ptr{Uint8},SubString(s,i)), base)
    err == 0 || error("invalid big integer: $(repr(s))")
    return sgn < 0 ? -z : z
end

function ZZ(x::Union(Clong, Int32))
    z = ZZ()
    ccall((:__fmpz_set_si, :libflint), Void, (Ptr{ZZ}, Clong), &z, x)
    return z
end

function ZZ(x::Union(Culong, Uint32))
    z = ZZ()
    ccall((:__fmpz_set_ui, :libflint), Void, (Ptr{ZZ}, Culong), &z, x)
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
            return ZZ(Base.convert(Clong, x))
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
            return ZZ(Base.convert(Culong, x))
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

# Metaprogram to define functions +, -, *, fdiv, div, mod, rem, gcd, lcm, &, |, $

for (fJ, fC) in ((:+, :add), (:-,:sub), (:*, :mul),
                 (:fdiv, :fdiv_q), (:div, :tdiv_q), (:mod, :mod),
                 (:gcd, :gcd), (:lcm, :lcm),
                 (:&, :and), (:|, :or), (:$, :xor))
    @eval begin
        function ($fJ)(x::ZZ, y::ZZ)
            z = ZZ()
            ccall(($(string(:fmpz_,fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &x, &y)
            return z
        end
    end
end

###########################################################################################
#
#   Multi operators (e.g. d = +(a, b, c) instead of d = a + b + c)
#
###########################################################################################

for (fJ, fC) in ((:+, :add), (:*, :mul), (:&, :and), (:|, :or), (:$, :xor))
    @eval begin
        function ($fJ)(a::ZZ, b::ZZ, c::ZZ)
            z = ZZ()
            ccall(($(string(:fmpz_,fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &a, &b)
            ccall(($(string(:fmpz_,fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &z, &c)
            return z
        end
        function ($fJ)(a::ZZ, b::ZZ, c::ZZ, d::ZZ)
            z = ZZ()
            ccall(($(string(:fmpz_,fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &a, &b)
            ccall(($(string(:fmpz_,fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &z, &c)
            ccall(($(string(:fmpz_,fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &z, &d)
            return z
        end
        function ($fJ)(a::ZZ, b::ZZ, c::ZZ, d::ZZ, e::ZZ)
            z = ZZ()
            ccall(($(string(:fmpz_,fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &a, &b)
            ccall(($(string(:fmpz_,fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &z, &c)
            ccall(($(string(:fmpz_,fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &z, &d)
            ccall(($(string(:fmpz_,fC)), :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &z, &e)
            return z
        end
    end
end

###########################################################################################
#
#   Specialisations of binary operators when one operand is base type (ad hoc polymorphism)
#
###########################################################################################

function +(x::ZZ, c::Culong)
    z = ZZ()
    ccall((:fmpz_add_ui, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Culong), &z, &x, c)
    return z
end

+(c::Culong, x::ZZ) = x + c

+(c::CulongMax, x::ZZ) = x + convert(Culong, c)

+(x::ZZ, c::CulongMax) = x + convert(Culong, c)

+(x::ZZ, c::ClongMax) = c < 0 ? -(x, convert(Culong, -c)) : x + convert(Culong, c)

+(c::ClongMax, x::ZZ) = c < 0 ? -(x, convert(Culong, -c)) : x + convert(Culong, c)

function -(x::ZZ, c::Culong)
    z = ZZ()
    ccall((:fmpz_sub_ui, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Culong), &z, &x, c)
    return z
end

function -(c::Culong, x::ZZ)
    z = ZZ()
    ccall((:fmpz_sub_ui, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Culong), &z, &x, c)
    ccall((:__fmpz_neg, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}), &z, &z)
    return z
end

-(x::ZZ, c::CulongMax) = -(x, convert(Culong, c))

-(c::CulongMax, x::ZZ) = -(convert(Culong, c), x)

-(x::ZZ, c::ClongMax) = c < 0 ? +(x, convert(Culong, -c)) : -(x, convert(Culong, c))

-(c::ClongMax, x::ZZ) = c < 0 ? -(x + convert(Culong, -c)) : -(convert(Culong, c), x)

function *(x::ZZ, c::Culong)
    z = ZZ()
    ccall((:fmpz_mul_ui, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Culong), &z, &x, c)
    return z
end

*(c::Culong, x::ZZ) = x * c

*(c::CulongMax, x::ZZ) = x * convert(Culong, c)

*(x::ZZ, c::CulongMax) = x * convert(Culong, c)

function *(x::ZZ, c::Clong)
    z = ZZ()
    ccall((:fmpz_mul_si, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Clong), &z, &x, c)
    return z
end

*(c::Clong, x::ZZ) = x * c

*(x::ZZ, c::ClongMax) = x * convert(Clong, c)

*(c::ClongMax, x::ZZ) = x * convert(Clong, c)

function <<(x::ZZ, c::Int32)
    c < 0 && throw(DomainError())
    c == 0 && return x
    z = ZZ()
    ccall((:fmpz_mul_2exp, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Culong), &z, &x, c)
    return z
end

function >>(x::ZZ, c::Int32)
    c < 0 && throw(DomainError())
    c == 0 && return x
    z = ZZ()
    ccall((:fmpz_fdiv_q_2exp, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Culong), &z, &x, c)
    return z
end


###########################################################################################
#
#   Unary operators, e.g. -ZZ(12), ~ZZ(12)
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

###########################################################################################
#
#   Arithmetic functions
#
###########################################################################################

function abs(x::ZZ)
    z = ZZ()
    ccall((:fmpz_abs, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}), &z, &x)
    return z
end

function divrem(x::ZZ, y::ZZ)
    z1 = ZZ()
    z2 = ZZ()
    ccall((:fmpz_tdiv_qr, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z1, &z2, &x, &y)
    z1, z2
end

function isqrt(x::ZZ)
    z = ZZ()
    ccall((:fmpz_sqrt, :libflint), Void, (Ptr{ZZ}, Ptr{ZZ}), &z, &x)
    return z
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

^(x::Int, y::ZZ) = zz_pow(ZZ(x), y)

function xgcd(a::ZZ, b::ZZ)
    if b == 0 # shortcut this to ensure consistent results with xgcd(a,b)
        return a < 0 ? (-a, -one(ZZ), zero(ZZ)) : (a, one(ZZ), zero(ZZ))
    end
    g = ZZ()
    s = ZZ()
    t = ZZ()
    ccall((:fmpz_xgcd, :libflint), Void,
        (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}),
        &g, &s, &t, &a, &b)
    if t == 0
        # work around a difference in some versions of GMP
        if a == b
            return g, t, s
        elseif abs(a) == abs(b)
            return g, t, -s
        end
    end
    g, s, t
end

###########################################################################################
#
#   Comparison
#
###########################################################################################

function cmp(x::ZZ, y::ZZ)
    ccall((:fmpz_cmp, :libflint), Int32, (Ptr{ZZ}, Ptr{ZZ}), &x, &y)
end

==(x::ZZ, y::ZZ) = cmp(x,y) == 0

<=(x::ZZ, y::ZZ) = cmp(x,y) <= 0

>=(x::ZZ, y::ZZ) = cmp(x,y) >= 0

<(x::ZZ, y::ZZ) = cmp(x,y) < 0

>(x::ZZ, y::ZZ) = cmp(x,y) > 0


###########################################################################################
#
#   Bit fiddling
#
###########################################################################################

popcount(x::ZZ) = int(ccall((:fmpz_popcnt, :libflint), Culong, (Ptr{ZZ},), &x))

prevpow2(x::ZZ) = x < 0 ? -prevpow2(-x) : (x <= 2 ? x : one(ZZ) << (num_digits(x, 2) - 1))

nextpow2(x::ZZ) = x < 0 ? -nextpow2(-x) : (x <= 2 ? x : one(ZZ) << num_digits(x - 1, 2))

###########################################################################################
#
#   Modular arithmetic
#
###########################################################################################

function powmod(x::ZZ, p::ZZ, m::ZZ)
    p < 0 && throw(DomainError())
    r = ZZ()
    ccall((:fmpz_powm, :libflint), Void,
          (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}),
          &r, &x, &p, &m)
    return m < 0 && r > 0 ? r + m : r # choose sign consistent with mod(x^p, m)
end

powmod(x::ZZ, p::Integer, m::ZZ) = powmod(x, ZZ(p), m)

powmod(x::ZZ, p::Integer, m::Integer) = powmod(x, ZZ(p), ZZ(m))

function invmod(x::ZZ, y::ZZ)
    z = ZZ()
    y = abs(y)
    if y == 1
        return ZZ(0)
    end
    if (y == 0 || ccall((:fmpz_invmod, :libflint), Cint, (Ptr{ZZ}, Ptr{ZZ}, Ptr{ZZ}), &z, &x, &y) == 0)
        error("no inverse exists")
    end
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

function fac(x::Uint)
    x < 0 && return ZZ(0)
    z = ZZ()
    ccall((:fmpz_fac_ui, :libflint), Void, (Ptr{ZZ}, Culong), &z, x)
    return z
end

fac(n::Integer) = n < 0 ? throw(DomainError()) : fac(uint(n))

function binom(n::Uint, k::Uint)
    z = ZZ()
    ccall((:fmpz_bin_uiui, :libflint), Void, (Ptr{ZZ}, Culong, Culong), &z, n, k)
    return z
end

binom(n::Integer, k::Integer) = k < 0 ? throw(DomainError()) : binom(uint(n), uint(k))

isprime(x::ZZ) = ccall((:fmpz_is_prime, :libflint), Cint, (Ptr{ZZ},), &x) != 0

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

dec(n::ZZ) = base(n, 8)

hex(n::ZZ) = base(n, 16)

function base(n::ZZ, b::Integer)
    2 <= b <= 62 || error("invalid base: $b")
    p = ccall((:fmpz_get_str,:libflint), Ptr{Uint8}, (Ptr{Uint8}, Cint, Ptr{ZZ}), C_NULL, b, &n)
    len = int(ccall(:strlen, Csize_t, (Ptr{Uint8},), p))
    ASCIIString(pointer_to_array(p,len,true))
end

function num_digits_internal(x::ZZ, b::Integer = 10)
    # fmpz_sizeinbase might return an answer 1 too big
    n = int(ccall((:fmpz_sizeinbase, :libflint), Culong, (Ptr{ZZ}, Int32), &x, b))
    abs(x) < ZZ(b)^(n - 1) ? n - 1 : n
end

num_digits(x::ZZ, b::Integer = 10) = x == 0 ? 1 : num_digits_internal(x,b)

