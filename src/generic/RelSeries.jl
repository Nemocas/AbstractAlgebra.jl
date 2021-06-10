###############################################################################
#
#   RelSeries.jl : Generic power series over rings, capped relative precision
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{RelSeries{T}}) where T <: RingElement = RelSeriesRing{T}

elem_type(::Type{RelSeriesRing{T}}) where T <: RingElement = RelSeries{T}

@doc Markdown.doc"""
    rel_series_type(::Type{T}) where T <: RingElement

Return the type of a relative series whose coefficients have the given type.
"""
rel_series_type(::Type{T}) where T <: RingElement = RelSeries{T}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function normalise(a::RelSeries, len::Int)
   while len > 0 && iszero(a.coeffs[len])
      len -= 1
   end
   return len
end

function polcoeff(a::RelSeries, n::Int)
   n < 0  && throw(DomainError(n, "n must be >= 0"))
   return n >= pol_length(a) ? zero(base_ring(a)) : a.coeffs[n + 1]
end

@doc Markdown.doc"""
    gen(R::RelSeriesRing)

Return the generator of the power series ring, i.e. $x + O(x^{n + 1})$ where
$n$ is the maximum precision of the power series ring $R$.
"""
function gen(R::RelSeriesRing)
   S = base_ring(R)
   return R([S(1)], 1, max_precision(R) + 1, 1)
end

function deepcopy_internal(a::RelSeries{T}, dict::IdDict) where T <: RingElement
   coeffs = Array{T}(undef, pol_length(a))
   for i = 1:pol_length(a)
      coeffs[i] = deepcopy(polcoeff(a, i - 1))
   end
   return parent(a)(coeffs, pol_length(a), precision(a), valuation(a))
end

function characteristic(a::RelSeriesRing{T}) where T <: RingElement
   return characteristic(base_ring(a))
end



###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::RelSeries)
   a.length = 0
   a.prec = parent(a).prec_max
   a.val = a.prec
   return a
end

function fit!(c::RelSeries{T}, n::Int) where T <: RingElement
   if length(c.coeffs) < n
      resize!(c.coeffs, n)
   end
   for i = pol_length(c) + 1:n
      c.coeffs[i] = zero(base_ring(c))
   end
   return nothing
end

function setcoeff!(c::RelSeries{T}, n::Int, a::T) where T <: RingElement
   if (a != 0 && precision(c) - valuation(c) > n) || n + 1 <= c.length
      fit!(c, n + 1)
      c.coeffs[n + 1] = a
      c.length = max(pol_length(c), n + 1)
      # don't normalise
   end
   return c
end

function mul!(c::RelSeries{T}, a::RelSeries{T}, b::RelSeries{T}) where T <: RingElement
   lena = pol_length(a)
   lenb = pol_length(b)
   aval = valuation(a)
   bval = valuation(b)
   prec = min(precision(a) - aval, precision(b) - bval)
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   if lena <= 0 || lenb <= 0
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
         d[i] = mul!(d[i], polcoeff(a, i - 1), polcoeff(b, 0))
      end
      if lenc > lena
         for i = 2:min(lenb, lenc - lena + 1)
            d[lena + i - 1] = mul!(d[lena + i - 1], polcoeff(a, lena - 1), polcoeff(b, i - 1))
         end
      end
      for i = 1:lena - 1
         if lenc > i
            for j = 2:min(lenb, lenc - i + 1)
               t = mul!(t, polcoeff(a, i - 1), polcoeff(b, j - 1))
               d[i + j - 1] = addeq!(d[i + j - 1], t)
            end
         end
      end
      c.coeffs = d
      c.length = normalise(c, lenc)
   end
   c.val = a.val + b.val
   c.prec = prec + c.val
   renormalize!(c)
   return c
end

function addeq!(c::RelSeries{T}, a::RelSeries{T}) where T <: RingElement
   lenc = pol_length(c)
   lena = pol_length(a)
   valc = valuation(c)
   vala = valuation(a)
   valr = min(vala, valc)
   precc = precision(c)
   preca = precision(a)
   prec = min(precc, preca)
   mina = min(vala + lena, prec)
   minc = min(valc + lenc, prec)
   lenr = max(mina, minc) - valr
   R = base_ring(c)
   fit!(c, lenr)
   if valc >= vala
      for i = min(lenr, lenc + valc - vala):-1:max(lena, valc - vala) + 1
         t = c.coeffs[i]
         c.coeffs[i] = c.coeffs[i - valc + vala]
         c.coeffs[i - valc + vala] = t
      end
      for i = min(lenr, lena):-1:valc - vala + 1
         c.coeffs[i] = add!(c.coeffs[i], c.coeffs[i - valc + vala], a.coeffs[i])
      end
      for i = 1:min(lenr, lena, valc - vala)
         c.coeffs[i] = deepcopy(a.coeffs[i])
      end
      for i = lena + 1:min(valc - vala, lenr)
         c.coeffs[i] = R()
      end
      for i = lenc + valc - vala + 1:min(lenr, lena)
         c.coeffs[i] = deepcopy(a.coeffs[i])
      end
   else
      for i = lenc + 1:min(vala - valc, lenr)
         c.coeffs[i] = R()
      end
      for i = vala - valc + 1:min(lenc, lenr, lena + vala - valc)
         c.coeffs[i] = addeq!(c.coeffs[i], a.coeffs[i - vala + valc])
      end
      for i = max(lenc, vala - valc) + 1:min(lena + vala - valc, lenr)
         c.coeffs[i] = deepcopy(a.coeffs[i - vala + valc])
      end
   end
   c.length = normalise(c, lenr)
   c.prec = prec
   c.val = valr
   renormalize!(c)
   return c
end

function add!(c::RelSeries{T}, a::RelSeries{T}, b::RelSeries{T}) where T <: RingElement
   if c === a
      return addeq!(c, b)
   elseif c === b
      return addeq!(c, a)
   end
   lena = pol_length(a)
   lenb = pol_length(b)
   valb = valuation(b)
   vala = valuation(a)
   valr = min(vala, valb)
   precb = precision(b)
   preca = precision(a)
   prec = min(precb, preca)
   mina = min(vala + lena, prec)
   minb = min(valb + lenb, prec)
   lenr = max(mina, minb) - valr
   R = base_ring(c)
   fit!(c, lenr)
   c.prec = prec
   c.val = valr
   if vala > valb
      for i = 1:min(lenr, min(lenb, vala - valb))
         c.coeffs[i] = deepcopy(b.coeffs[i])
      end
      for i = lenb + 1:min(lenr, vala - valb)
         c.coeffs[i] = R()
      end
      for i = vala - valb + 1:min(lenr, lenb)
         c.coeffs[i] = add!(c.coeffs[i], polcoeff(a, i - vala + valb - 1), b.coeffs[i])
      end
      for i = max(lenb, vala - valb) + 1:min(lenr, lena + vala - valb)
         c.coeffs[i] = deepcopy(a.coeffs[i - vala + valb])
      end
      for i = lena + vala - valb + 1:min(lenr, lenb)
         c.coeffs[i] = deepcopy(b.coeffs[i])
      end
   else
      for i = 1:min(lenr, min(lena, valb - vala))
         c.coeffs[i] = deepcopy(a.coeffs[i])
      end
      for i = lena + 1:min(lenr, valb - vala)
         c.coeffs[i] = R()
      end
      for i = valb - vala + 1:min(lenr, lena)
         c.coeffs[i] = add!(c.coeffs[i], a.coeffs[i], polcoeff(b, i - valb + vala - 1))
      end
      for i = max(lena, valb - vala) + 1:min(lenr, lenb + valb - vala)
         c.coeffs[i] = deepcopy(b.coeffs[i - valb + vala])
      end
      for i = lenb + valb - vala + 1:min(lenr, lena)
         c.coeffs[i] = deepcopy(a.coeffs[i])
      end
   end
   c = set_length!(c, normalise(c, lenr))
   renormalize!(c)
   return c
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{RelSeries{T}}, ::Type{RelSeries{T}}) where T <: RingElement = RelSeries{T}

function promote_rule(::Type{RelSeries{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? RelSeries{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (R::RelSeriesRing{T})(b::RingElement) where T <: RingElement
   return R(base_ring(R)(b))
end

function (R::RelSeriesRing{T})() where T <: RingElement
   z = RelSeries{T}(Array{T}(undef, 0), 0, R.prec_max, R.prec_max)
   z.parent = R
   return z
end

function (R::RelSeriesRing{T})(b::Union{Integer, Rational, AbstractFloat}) where T <: RingElement
   if b == 0
      z = RelSeries{T}(Array{T}(undef, 0), 0, R.prec_max, R.prec_max)
   else
      z = RelSeries{T}([base_ring(R)(b)], 1, R.prec_max, 0)
   end
   z.parent = R
   return z
end

function (R::RelSeriesRing{T})(b::T) where {T <: RingElem}
   parent(b) != base_ring(R) && error("Unable to coerce to power series")
   if iszero(b)
      z = RelSeries{T}(Array{T}(undef, 0), 0, R.prec_max, R.prec_max)
   else
      z = RelSeries{T}([b], 1, R.prec_max, 0)
   end
   z.parent = R
   return z
end

function (R::RelSeriesRing{T})(b::RelSeriesElem{T}) where T <: RingElement
   parent(b) != R && error("Unable to coerce power series")
   return b
end

function (R::RelSeriesRing{T})(b::Array{T, 1}, len::Int, prec::Int, val::Int) where T <: RingElement
   if length(b) > 0
      parent(b[1]) != base_ring(R) && error("Unable to coerce to power series")
   end
   z = RelSeries{T}(b, len, prec, val)
   z.parent = R
   return z
end

function (R::RelSeriesRing{T})(b::Array{S, 1}, len::Int, prec::Int, val::Int) where {S <: RingElement, T <: RingElement}
   R0 = base_ring(R)
   lenb = length(b)
   entries = Array{T}(undef, lenb)
   for i = 1:lenb
      entries[i] = R0(b[i])
   end
   z = RelSeries{T}(entries, len, prec, val)
   z.parent = R
   return z
end

###############################################################################
#
#   PowerSeriesRing constructor
#
###############################################################################

function PowerSeriesRing(R::AbstractAlgebra.Ring, prec::Int, s::Symbol; cached=true, model=:capped_relative)
   T = elem_type(R)

   if model == :capped_relative
      parent_obj = RelSeriesRing{T}(R, prec, s, cached)
   elseif model == :capped_absolute
      parent_obj = AbsSeriesRing{T}(R, prec, s, cached)
   else
      error("Unknown model")
   end

   return parent_obj, gen(parent_obj)
end