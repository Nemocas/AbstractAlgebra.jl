###############################################################################
#
#   Poly.jl : Generic polynomials over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{Poly{T}}) where T <: RingElement = PolyRing{T}

elem_type(::Type{PolyRing{T}}) where T <: RingElement = Poly{T}

base_ring(R::PolyRing{T}) where T <: RingElement = R.base_ring::base_ring_type(typeof(R))

parent(a::Poly) = a.parent

var(a::PolyRing) = a.S

is_trivial(a::PolyRing) = a.istrivial

###############################################################################
#
#   Basic manipulation
#
###############################################################################

length(a::Poly) = a.length

function setcoeff!(c::Poly{T}, n::Int, a::S) where {T <: RingElement, S <: RingElement}
   if !iszero(a) || n + 1 <= length(c)
      fit!(c, n + 1)
      c.coeffs[n + 1] = a
      c.length = max(length(c), n + 1)
      # don't normalise
   end
   return c
end

function set_coefficient!(c::Poly{T}, n::Int, a::T) where T <: RingElement
   c = setcoeff!(c, n, a)
   c = set_length!(c, normalise(c, length(c)))
   return c
end

function set_coefficient!(c::Poly{T}, n::Int, a::U) where {T <: RingElement, U <: Integer}
   c = setcoeff!(c, n, base_ring(c)(a))
   c = set_length!(c, normalise(c, length(c)))
   return c
end

function set_coefficient!(c::Poly{T}, n::Int, a::T) where T <: Integer
   c = setcoeff!(c, n, a)
   c = set_length!(c, normalise(c, length(c)))
   return c
end

function normalise(a::Poly, n::Int)
   while n > 0 && iszero(a.coeffs[n])
      n -= 1
   end
   return n
end

coeff(a::Poly, n::Int) = n >= length(a) ? base_ring(a)(0) : a.coeffs[n + 1]

function deepcopy_internal(a::Poly{T}, dict::IdDict) where T <: RingElement
   coeffs = Vector{T}(undef, length(a))
   for i = 1:length(a)
      coeffs[i] = deepcopy_internal(a.coeffs[i], dict)
   end
   return parent(a)(coeffs)
end

Base.copy(f::Poly) = deepcopy(f)

###############################################################################
#
#   Karatsuba multiplication
#
###############################################################################

function use_karamul(a::Poly{BigInt}, b::Poly{BigInt})
   minlen = min(length(a), length(b))
   if minlen == 0
      return false
   end
   if minlen > 175
      return true
   end
   bits = 0
   for i = 1:length(a)
      bits += ndigits(a.coeffs[i], base=2)
   end
   for i = 1:length(b)
      bits += ndigits(b.coeffs[i], base=2)
   end
   return minlen*div(bits, length(a) + length(b)) > 30000
end

function use_karamul(a::Poly{Rational{BigInt}}, b::Poly{Rational{BigInt}})
   minlen = min(length(a), length(b))
   if minlen == 0
      return false
   end
   if minlen > 17
      return true
   end
   bits = 0
   for i = 1:length(a)
      bits += ndigits(numerator(a.coeffs[i]), base=2)
      bits += ndigits(denominator(a.coeffs[i]), base=2)
   end
   for i = 1:length(b)
      bits += ndigits(numerator(b.coeffs[i]), base=2)
      bits += ndigits(denominator(b.coeffs[i]), base=2)
   end
   return minlen^1.7*div(bits, 2*(length(a) + length(b))) > 48500
end

function use_karamul(a::Poly{GFElem{Int}}, b::Poly{GFElem{Int}})
   return min(length(a), length(b)) > 75
end

function use_karamul(a::Poly{GFElem{BigInt}}, b::Poly{GFElem{BigInt}})
   minlen = min(length(a), length(b))
   bits = ndigits(characteristic(parent(a)), base=2)
   return minlen^2*bits > 2000
end

@doc raw"""
    mul_karatsuba(a::Poly{T}, b::Poly{T}, cutoff::Int) where T <: RingElement

Return $a \times b$ using the Karatsuba algorithm recursively for problems of
size roughly greater than `cutoff`.
"""
function mul_karatsuba(a::Poly{T}, b::Poly{T}, cutoff::Int) where T <: RingElement
   alen = length(a)
   blen = length(b)
   (alen < 1 || blen < 1) && return zero(parent(a))
   zlen = alen + blen - 1
   zcoeffs = Vector{T}(undef, zlen)
   AbstractAlgebra.DensePoly.mullow_fast!(zcoeffs, zlen,
                          a.coeffs, alen, b.coeffs, blen, base_ring(a), cutoff)
   z = parent(a)(zcoeffs)
   z = set_length!(z, normalise(z, zlen))
   return z
end

@doc raw"""
    mullow_karatsuba(a::Poly{T}, b::Poly{T}, n::Int, cutoff::Int) where T <: RingElement

Return $a \times b$ truncated to $n$ terms using the Karatsuba algorithm
recursively for problems of size roughly greater than `cutoff`.
"""
function mullow_karatsuba(a::Poly{T}, b::Poly{T}, n::Int, cutoff::Int) where T <: RingElement
   alen = length(a)
   blen = length(b)
   (n < 1 || alen < 1 || blen < 1) && return zero(parent(a))
   zlen = min(alen + blen - 1, n)
   zcoeffs = Vector{T}(undef, zlen)
   AbstractAlgebra.DensePoly.mullow_fast!(zcoeffs, zlen,
                          a.coeffs, alen, b.coeffs, blen, base_ring(a), cutoff)
   z = parent(a)(zcoeffs)
   z = set_length!(z, normalise(z, zlen))
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function set_length!(c::Poly{T}, n::Int) where T <: RingElement
   if n < c.length
      for i = n + 1:c.length
         c.coeffs[i] = zero!(c.coeffs[i])
      end
   end
   c.length = n
   return c
end

function fit!(c::Poly{T}, n::Int) where T <: RingElement
   if length(c.coeffs) < n
      resize!(c.coeffs, n)
      for i = length(c) + 1:n
         c.coeffs[i] = zero(base_ring(c))
      end
   end
   return nothing
end

function zero!(c::Poly{T}) where T <: RingElement
   c = set_length!(c, 0)
   return c
end

function one!(c::Poly{T}) where T <: RingElement
   is_trivial(parent(c)) && return zero!(c)
   fit!(c, 1)
   c = set_length!(c, 1)
   c.coeffs[1] = one(base_ring(c))
   return c
end

function neg!(a::Poly{T}) where T <: RingElement
   for i in 1:length(a)
      a.coeffs[i] = neg!(a.coeffs[i])
   end
   return a
end

function neg!(z::Poly{T}, a::Poly{T}) where T <: RingElement
   fit!(z, length(a))
   z = set_length!(z, length(a))
   for i in 1:length(a)
      z.coeffs[i] = neg!(z.coeffs[i], a.coeffs[i])
   end
   return z
end

function mul!(c::Poly{T}, a::Poly{T}, b::Poly{T}) where T <: RingElement
   lena = length(a)
   lenb = length(b)

   if lena == 0 || lenb == 0
      c = set_length!(c, 0)
   else
      if a === c
         a = deepcopy(a)
      end
      if b === c
         b = deepcopy(b)
      end

      t = base_ring(a)()

      lenc = lena + lenb - 1
      fit!(c, lenc)

      for i = 1:lena
         c.coeffs[i] = mul!(c.coeffs[i], coeff(a, i - 1), coeff(b, 0))
      end

      for i = 2:lenb
         c.coeffs[lena + i - 1] = mul!(c.coeffs[lena + i - 1], coeff(a, lena - 1), coeff(b, i - 1))
      end

      for i = 1:lena - 1
         for j = 2:lenb
            t = mul!(t, coeff(a, i - 1), coeff(b, j - 1))
            c.coeffs[i + j - 1] = add!(c.coeffs[i + j - 1], t)
         end
      end

      c = set_length!(c, normalise(c, lenc))
   end
   return c
end

function add!(c::Poly{T}, a::Poly{T}) where T <: RingElement
   lenc = length(c)
   lena = length(a)
   len = max(lenc, lena)
   fit!(c, len)
   for i = 1:lena
      c.coeffs[i] = add!(c.coeffs[i], coeff(a, i - 1))
   end
   c = set_length!(c, normalise(c, len))
   return c
end

function add!(c::Poly{T}, a::Poly{T}, b::Poly{T}) where T <: RingElement
   lena = length(a)
   lenb = length(b)
   len = max(lena, lenb)
   fit!(c, len)
   i = 1
   while i <= min(lena, lenb)
      c.coeffs[i] = add!(c.coeffs[i], coeff(a, i - 1), coeff(b, i - 1))
      i += 1
   end
   while i <= lena
      # mutating operators must ensure they don't introduce new aliasing
      c = setcoeff!(c, i - 1, deepcopy(coeff(a, i - 1)))
      i += 1
   end
   while i <= lenb
      # mutating operators must ensure they don't introduce new aliasing
      c = setcoeff!(c, i - 1, deepcopy(coeff(b, i - 1)))
      i += 1
   end
   c = set_length!(c, normalise(c, len))
   return c
end

###############################################################################
#
#   Iterators
#
###############################################################################

function Base.iterate(a::MPolyExponentVectors{T}, st::Int = -1) where T <: AbstractAlgebra.PolyRingElem
   st += 1
   if st > degree(a.poly)
       return nothing
   else
       return Int[st], st
   end
end

###############################################################################
#
#   Build context
#
###############################################################################

# used to provide a uniform interface to build ctxs for poly and mpoly
mutable struct PolyBuildCtx{T, S}
   poly::Vector{T}
   parent::S
end

# TODO the MPolyBuildCtx function should be renamed BuildCtx
function MPolyBuildCtx(R::AbstractAlgebra.PolyRing{T}) where T
   S = typeof(R)
   return PolyBuildCtx{T, S}(T[], R)
end

function push_term!(B::PolyBuildCtx{T, S}, c::T, expv::Vector{Int}) where {S, T}
   if length(expv) != 1
      error("length of exponent vector should match the number of variables")
   end
   i = expv[1] + 1
   if iszero(c)
      return B
   end
   while i > length(B.poly)
      push!(B.poly, zero(coefficient_ring(B.parent))::T)
   end
   B.poly[i] += c
   return B
end

function finish(B::PolyBuildCtx{T, S}) where {T, S}
   res = B.parent(B.poly)  # construction of univar from array of coeffs
   B.poly = T[]            # apparently res can have ownership of the B.poly
   return res
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{Poly{T}}, ::Type{Poly{T}}) where T <: RingElement = Poly{T}

function promote_rule(::Type{Poly{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? Poly{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::PolyRing{T})(b::RingElement) where T <: RingElement
   return a(base_ring(a)(b))
end

function (a::PolyRing{T})() where T <: RingElement
   z = Poly{T}()
   z.parent = a
   return z
end

function (a::PolyRing{T})(b::JuliaRingElement) where T <: RingElement
   z = Poly{T}(base_ring(a)(b))
   z.parent = a
   return z
end

function (a::PolyRing{T})(b::T) where T <: RingElement
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = Poly{T}(b)
   z.parent = a
   return z
end

function (a::PolyRing{T})(b::AbstractAlgebra.PolyRingElem{T}) where T <: RingElement
   parent(b) != a && error("Unable to coerce polynomial")
   return b
end

function (a::PolyRing{T})(b::Vector{T}) where T <: RingElement
   R = base_ring(a)
   for i = 1:length(b)
      b[i] = R(b[i])
   end
   z = Poly{T}(b)
   z.parent = a
   return z
end

function (a::PolyRing{T})(b::Vector{S}) where {S <: RingElement, T <: RingElement}
   R = base_ring(a)
   len = length(b)
   entries = Vector{T}(undef, len)
   for i = 1:length(b)
      entries[i] = R(b[i])
   end
   z = Poly{T}(entries)
   z.parent = a
   return z
end

# Functions to remove ambiguities on julia 0.7
function (a::PolyRing{T})(b::T) where {T <: Rational}
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = Poly{T}(b)
   z.parent = a
   return z
end

function (a::PolyRing{T})(b::T) where {T <: AbstractFloat}
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = Poly{T}(b)
   z.parent = a
   return z
end

function (a::PolyRing{T})(b::T) where {T <: Integer}
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = Poly{T}(b)
   z.parent = a
   return z
end
