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

parent_type(::Type{AbsSeries{T}}) where T <: RingElement = AbsSeriesRing{T}

elem_type(::Type{AbsSeriesRing{T}}) where T <: RingElement = AbsSeries{T}

@doc Markdown.doc"""
    abs_series_type(::Type{T}) where T <: RingElement

Return the type of an absolute series whose coefficients have the given type.
"""
abs_series_type(::Type{T}) where T <: RingElement = AbsSeries{T}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

@doc Markdown.doc"""
    gen(R::AbsSeriesRing{T}) where T <: RingElement

Return the generator of the power series ring, i.e. $x + O(x^n)$ where
$n$ is the precision of the power series ring $R$.
"""
function gen(R::AbsSeriesRing{T}) where T <: RingElement
   S = base_ring(R)
   return R([S(0), S(1)], 2, max_precision(R))
end

@doc Markdown.doc"""
    max_precision(R::AbsSeriesRing)

Return the maximum absolute precision of power series in the given power
series ring.
"""
max_precision(R::AbsSeriesRing) = R.prec_max

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
   coeffs = Array{T}(undef, length(a))
   for i = 1:length(a)
      coeffs[i] = deepcopy(coeff(a, i - 1))
   end
   return parent(a)(coeffs, length(a), precision(a))
end

function characteristic(a::AbsSeriesRing{T}) where T <: RingElement
   return characteristic(base_ring(a))
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

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
      t = base_ring(a)()

      for i = 1:min(lena, lenc)
         d[i] = mul!(d[i], coeff(a, i - 1), coeff(b, 0))
      end

      if lenc > lena
         for i = 2:min(lenb, lenc - lena + 1)
            d[lena + i - 1] = mul!(d[lena + i - 1], coeff(a, lena - 1), coeff(b, i - 1))
         end
      end

      for i = 1:lena - 1
         if lenc > i
            for j = 2:min(lenb, lenc - i + 1)
               t = mul!(t, coeff(a, i - 1), coeff(b, j - 1))
               d[i + j - 1] = addeq!(d[i + j - 1], t)
            end
         end
      end

      c.coeffs = d
      c.length = normalise(c, lenc)
   end
   c.prec = prec
   return c
end

function addeq!(c::AbsSeries{T}, a::AbsSeries{T}) where T <: RingElement
   lenc = length(c)
   lena = length(a)

   prec = min(precision(a), precision(c))

   lena = min(lena, prec)
   lenc = min(lenc, prec)

   len = max(lenc, lena)
   fit!(c, len)
   for i = 1:lena
      c.coeffs[i] = addeq!(c.coeffs[i], coeff(a, i - 1))
   end
   c.length = normalise(c, len)
   c.prec = prec
   return c
end

function add!(c::AbsSeries{T}, a::AbsSeries{T}, b::AbsSeries{T}) where T <: RingElement
   if c === a
      return addeq!(c, b)
   elseif c === b
      return addeq!(c, a)
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

function (a::AbsSeriesRing{T} where T <: RingElement)(b::RingElement)
   return a(base_ring(a)(b))
end

function (a::AbsSeriesRing{T})() where T <: RingElement
   z = AbsSeries{T}(Array{T}(undef, 0), 0, a.prec_max)
   z.parent = a
   return z
end

function (a::AbsSeriesRing{T})(b::Union{Integer, Rational, AbstractFloat}) where T <: RingElement
   if b == 0
      z = AbsSeries{T}(Array{T}(undef, 0), 0, a.prec_max)
   else
      z = AbsSeries{T}([base_ring(a)(b)], 1, a.prec_max)
   end
   z.parent = a
   return z
end

function (a::AbsSeriesRing{T})(b::T) where {T <: RingElem}
   parent(b) != base_ring(a) && error("Unable to coerce to power series")
   if iszero(b)
      z = AbsSeries{T}(Array{T}(undef, 0), 0, a.prec_max)
   else
      z = AbsSeries{T}([b], 1, a.prec_max)
   end
   z.parent = a
   return z
end

function (a::AbsSeriesRing{T})(b::AbsSeriesElem{T}) where T <: RingElement
   parent(b) != a && error("Unable to coerce power series")
   return b
end

function (a::AbsSeriesRing{T})(b::Array{T, 1}, len::Int, prec::Int) where T <: RingElement
   if length(b) > 0
      parent(b[1]) != base_ring(a) && error("Unable to coerce to power series")
   end
   z = AbsSeries{T}(b, len, prec)
   z.parent = a
   return z
end

function (a::AbsSeriesRing{T})(b::Array{S, 1}, len::Int, prec::Int) where {S <: RingElement, T <: RingElement}
   R = base_ring(a)
   lenb = length(b)
   entries = Array{T}(undef, lenb)
   for i = 1:lenb
      entries[i] = R(b[i])
   end
   z = AbsSeries{T}(entries, len, prec)
   z.parent = a
   return z
end

