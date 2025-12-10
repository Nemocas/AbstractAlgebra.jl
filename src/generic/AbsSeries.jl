###############################################################################
#
#   AbsSeries.jl : Generic power series over rings, capped absolute precision
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{AbsSeries{T}}) where T <: RingElement = AbsPowerSeriesRing{T}

elem_type(::Type{AbsPowerSeriesRing{T}}) where T <: RingElement = AbsSeries{T}

@doc raw"""
    abs_series_type(::Type{T}) where T <: RingElement

Return the type of an absolute series whose coefficients have the given type.
"""
abs_series_type(::Type{T}) where T <: RingElement = AbsSeries{T}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

@doc raw"""
    gen(R::AbsPowerSeriesRing{T}) where T <: RingElement

Return the generator of the power series ring, i.e. $x + O(x^n)$ where
$n$ is the precision of the power series ring $R$.
"""
function gen(R::AbsPowerSeriesRing{T}) where T <: RingElement
   S = base_ring(R)
   return R([S(0), S(1)], 2, max_precision(R))
end

number_of_variables(R::AbsPowerSeriesRing) = 1

@doc raw"""
    max_precision(R::AbsPowerSeriesRing)

Return the maximum absolute precision of power series in the given power
series ring.
"""
max_precision(R::AbsPowerSeriesRing) = R.prec_max

function normalise(a::AbsSeries, len::Int)
   while len > 0 && iszero(a.coeffs[len])
      len -= 1
   end
   return len
end

function coeff(a::AbsSeries, n::Int)
   n < 0  && throw(DomainError(n, "n must be >= 0"))
   return n >= length(a) ? zero(base_ring(a)) : a.coeffs[n + 1]
end

polcoeff(a::AbsSeries, n::Int) = coeff(a, n)

function deepcopy_internal(a::AbsSeries{T}, dict::IdDict) where T <: RingElement
   coeffs = Vector{T}(undef, length(a))
   for i = 1:length(a)
      coeffs[i] = deepcopy_internal(coeff(a, i - 1), dict)
   end
   return parent(a)(coeffs, length(a), precision(a))
end

function characteristic(a::AbsPowerSeriesRing{T}) where T <: RingElement
   return characteristic(base_ring(a))
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function mullow_fast_cutoff(a::AbsSeries{BigInt}, b::AbsSeries{BigInt})
   bits = 0
   for i = 1:length(a)
      bits += ndigits(a.coeffs[i], base=2)
   end
   for i = 1:length(b)
      bits += ndigits(b.coeffs[i], base=2)
   end
   bits = div(bits, length(a) + length(b))
   len = 2
   while len*bits <= 30000
      len *= 2
   end
   return len
end

function mullow_fast_cutoff(a::AbsSeries{Rational{BigInt}}, b::AbsSeries{Rational{BigInt}})
   bits = 0
   for i = 1:length(a)
      bits += ndigits(numerator(a.coeffs[i]), base=2)
      bits += ndigits(denominator(a.coeffs[i]), base=2)
   end
   for i = 1:length(b)
      bits += ndigits(numerator(b.coeffs[i]), base=2)
      bits += ndigits(denominator(b.coeffs[i]), base=2)
   end
   bits = div(bits, 2*(length(a) + length(b)))
   len = 2
   while len^1.7*bits <= 48500
      len *= 2
   end
   return len
end

function mullow_fast_cutoff(a::AbsSeries{GFElem{Int}}, b::AbsSeries{GFElem{Int}})
   return 75
end

function mullow_fast_cutoff(a::AbsSeries{GFElem{BigInt}}, b::AbsSeries{GFElem{BigInt}})
   bits = ndigits(characteristic(parent(a)), base=2)
   len = 2
   while len^2*bits <= 2000
      len *= 2
   end
   return len
end

# generic fallback
function mullow_fast_cutoff(a::T, b::T) where {S <: RingElement, T <: AbsSeries{S}}
   return 5
end

function *(a::AbsSeries{T}, b::AbsSeries{T}) where T <: RingElement
   check_parent(a, b)

   lena = length(a)
   lenb = length(b)

   aval = valuation(a)
   bval = valuation(b)

   prec = min(precision(a) + bval, precision(b) + aval)
   prec = min(prec, max_precision(parent(a)))

   lena = min(lena, prec)
   lenb = min(lenb, prec)

   if lena == 0 || lenb == 0
      return parent(a)(Vector{T}(undef, 0), 0, prec)
   end
   lenz = min(lena + lenb - 1, prec)
   d = Vector{T}(undef, lenz)
   cutoff = mullow_fast_cutoff(a, b)
   AbstractAlgebra.DensePoly.mullow_fast!(d, lenz,
                          a.coeffs, lena, b.coeffs, lenb, base_ring(a), cutoff)
   z = parent(a)(d, lenz, prec)
   z = set_length!(z, normalise(z, lenz))
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function truncate!(a::AbsSeries{T}, n::Int) where T <: RingElement
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   if precision(a) <= n
      return a
   end
   a.length = min(n, length(a))
   while length(a) != 0 && is_zero(coeff(a, length(a) - 1))
      a.length -= 1
   end
   a.prec = n
   return a
end

function zero!(c::AbsSeries{T}) where T <: RingElement
   c.length = 0
   c.prec = parent(c).prec_max
   return c
end

function fit!(c::AbsSeries{T}, n::Int) where T <: RingElement
   if length(c.coeffs) < n
      resize!(c.coeffs, n)
      for i = length(c) + 1:n
         c.coeffs[i] = zero(base_ring(c))
      end
   end
   return nothing
end

function setcoeff!(c::AbsSeries{T}, n::Int, a::T) where T <: RingElement
   if (!iszero(a) && precision(c) > n) || n + 1 <= c.length
      fit!(c, n + 1)
      c.coeffs[n + 1] = a
      c.length = max(length(c), n + 1)
      # don't normalise
   end
   return c
end

function mul!(c::AbsSeries{T}, a::AbsSeries{T}, b::AbsSeries{T}) where T <: RingElement
   lena = length(a)
   lenb = length(b)

   aval = valuation(a)
   bval = valuation(b)

   prec = min(precision(a) + bval, precision(b) + aval)
   prec = min(prec, max_precision(parent(c)))

   lena = min(lena, prec)
   lenb = min(lenb, prec)

   if lena == 0 || lenb == 0
      c.length = 0
   else
      lenc = min(lena + lenb - 1, prec)

      if c === a || c === b
         d = T[base_ring(c)() for i in 1:lenc]
      else
         fit!(c, lenc)
         d = c.coeffs
      end

      cutoff = mullow_fast_cutoff(a, b)
      AbstractAlgebra.DensePoly.mullow_fast!(d, lenc,
                          a.coeffs, lena, b.coeffs, lenb, base_ring(a), cutoff)

      c.coeffs = d
      c.length = normalise(c, lenc)
   end
   c.prec = prec
   return c
end

function add!(c::AbsSeries{T}, a::AbsSeries{T}) where T <: RingElement
   lenc = length(c)
   lena = length(a)

   prec = min(precision(a), precision(c))

   lena = min(lena, prec)
   lenc = min(lenc, prec)

   len = max(lenc, lena)
   fit!(c, len)
   for i = 1:lena
      c.coeffs[i] = add!(c.coeffs[i], coeff(a, i - 1))
   end
   c.length = normalise(c, len)
   c.prec = prec
   return c
end

function add!(c::AbsSeries{T}, a::AbsSeries{T}, b::AbsSeries{T}) where T <: RingElement
   if c === a
      return add!(c, b)
   elseif c === b
      return add!(c, a)
   end
   lena = length(a)
   lenb = length(b)
   prec = min(precision(a), precision(b))
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   lenc = max(lena, lenb)
   fit!(c, lenc)
   i = 1
   while i <= min(lena, lenb)
      c.coeffs[i] = coeff(a, i - 1) + coeff(b, i - 1)
      i += 1
   end
   while i <= lena
      c.coeffs[i] = deepcopy(coeff(a, i - 1))
      i += 1
   end
   while i <= lenb
      c.coeffs[i] = deepcopy(coeff(b, i - 1))
      i += 1
   end
   c.length = normalise(c, i - 1)
   c.prec = prec
   return c
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{AbsSeries{T}}, ::Type{AbsSeries{T}}) where T <: RingElement = AbsSeries{T}

function promote_rule(::Type{AbsSeries{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? AbsSeries{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::AbsPowerSeriesRing)(b::RingElement)
   return a(base_ring(a)(b))
end

function (a::AbsPowerSeriesRing{T})() where T <: RingElement
   z = AbsSeries{T}(Vector{T}(undef, 0), 0, a.prec_max)
   z.parent = a
   return z
end

function (a::AbsPowerSeriesRing{T})(b::JuliaRingElement) where T <: RingElement
   bb = base_ring(a)(b)
   if is_zero(bb)
      z = AbsSeries{T}(Vector{T}(undef, 0), 0, a.prec_max)
   else
      z = AbsSeries{T}([bb], 1, a.prec_max)
   end
   z.parent = a
   return z
end

function (a::AbsPowerSeriesRing{T})(b::T) where {T <: RingElem}
   parent(b) != base_ring(a) && error("Unable to coerce to power series")
   if iszero(b)
      z = AbsSeries{T}(Vector{T}(undef, 0), 0, a.prec_max)
   else
      z = AbsSeries{T}([b], 1, a.prec_max)
   end
   z.parent = a
   return z
end

function (a::AbsPowerSeriesRing{T})(b::AbsPowerSeriesRingElem{T}) where T <: RingElement
   parent(b) != a && error("Unable to coerce power series")
   return b
end

function (a::AbsPowerSeriesRing{T})(b::Vector{T}, len::Int, prec::Int) where T <: RingElement
   if length(b) > 0
      parent(b[1]) != base_ring(a) && error("Unable to coerce to power series")
   end
   z = AbsSeries{T}(b, len, prec)
   z.parent = a
   return z
end

function (a::AbsPowerSeriesRing{T})(b::Vector{S}, len::Int, prec::Int) where {S <: RingElement, T <: RingElement}
   R = base_ring(a)
   lenb = length(b)
   entries = Vector{T}(undef, lenb)
   for i = 1:lenb
      entries[i] = R(b[i])
   end
   z = AbsSeries{T}(entries, len, prec)
   z.parent = a
   return z
end

