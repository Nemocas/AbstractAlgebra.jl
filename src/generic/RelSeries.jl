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

parent_type(::Type{RelSeries{T}}) where T <: RingElement = RelPowerSeriesRing{T}

elem_type(::Type{RelPowerSeriesRing{T}}) where T <: RingElement = RelSeries{T}

@doc raw"""
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

@doc raw"""
    gen(R::RelPowerSeriesRing)

Return the generator of the power series ring, i.e. $x + O(x^{n + 1})$ where
$n$ is the maximum precision of the power series ring $R$.
"""
function gen(R::RelPowerSeriesRing)
   S = base_ring(R)
   return R([S(1)], 1, max_precision(R) + 1, 1)
end

number_of_variables(R::RelPowerSeriesRing) = 1

function deepcopy_internal(a::RelSeries{T}, dict::IdDict) where T <: RingElement
   coeffs = Vector{T}(undef, pol_length(a))
   for i = 1:pol_length(a)
      coeffs[i] = deepcopy_internal(polcoeff(a, i - 1), dict)
   end
   return parent(a)(coeffs, pol_length(a), precision(a), valuation(a))
end

function characteristic(a::RelPowerSeriesRing{T}) where T <: RingElement
   return characteristic(base_ring(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function mullow_fast_cutoff(a::RelSeries{BigInt}, b::RelSeries{BigInt})
   bits = 0
   for i = 1:pol_length(a)
      bits += ndigits(a.coeffs[i], base=2)
   end
   for i = 1:pol_length(b)
      bits += ndigits(b.coeffs[i], base=2)
   end
   bits = div(bits, pol_length(a) + pol_length(b))
   len = 2
   while len*bits <= 30000
      len *= 2
   end
   return len
end

function mullow_fast_cutoff(a::RelSeries{Rational{BigInt}}, b::RelSeries{Rational{BigInt}})
   bits = 0
   for i = 1:pol_length(a)
      bits += ndigits(numerator(a.coeffs[i]), base=2)
      bits += ndigits(denominator(a.coeffs[i]), base=2)
   end
   for i = 1:pol_length(b)
      bits += ndigits(numerator(b.coeffs[i]), base=2)
      bits += ndigits(denominator(b.coeffs[i]), base=2)
   end
   bits = div(bits, 2*(pol_length(a) + pol_length(b)))
   len = 2
   while len^1.7*bits <= 48500
      len *= 2
   end
   return len
end


function mullow_fast_cutoff(a::RelSeries{GFElem{Int}}, b::RelSeries{GFElem{Int}})
   return 75
end

function mullow_fast_cutoff(a::RelSeries{GFElem{BigInt}}, b::RelSeries{GFElem{BigInt}})
   bits = ndigits(characteristic(parent(a)), base=2)
   len = 2
   while len^2*bits <= 2000
      len *= 2
   end
   return len
end

# generic fallback
function mullow_fast_cutoff(a::T, b::T) where {S <: RingElement, T <: RelSeries{S}}
   return 5
end

function *(a::RelSeries{T}, b::RelSeries{T}) where T <: RingElement
   check_parent(a, b)
   lena = pol_length(a)
   lenb = pol_length(b)
   aval = valuation(a)
   bval = valuation(b)
   zval = aval + bval
   prec = min(precision(a) - aval, precision(b) - bval)
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   if lena == 0 || lenb == 0
      return parent(a)(Vector{T}(undef, 0), 0, prec + zval, zval)
   end
   t = base_ring(a)()
   lenz = min(lena + lenb - 1, prec)
   d = Vector{T}(undef, lenz)
   cutoff = mullow_fast_cutoff(a, b)
   AbstractAlgebra.DensePoly.mullow_fast!(d, lenz,
                          a.coeffs, lena, b.coeffs, lenb, base_ring(a), cutoff)
   z = parent(a)(d, lenz, prec + zval, zval)
   z = set_length!(z, normalise(z, lenz))
   renormalize!(z)
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function truncate!(a::RelSeries{T}, n::Int) where T <: RingElement
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   if precision(a) <= n
      return a
   end
   if n <= valuation(a)
      a = zero!(a)
      a.val = n
   else
      a.length = min(n - valuation(a), pol_length(a))
      while is_zero(polcoeff(a, pol_length(a) - 1))
         a.length -= 1
      end
   end
   a.prec = n
   return a
end

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
      cutoff = mullow_fast_cutoff(a, b)
      AbstractAlgebra.DensePoly.mullow_fast!(d, lenc,
                          a.coeffs, lena, b.coeffs, lenb, base_ring(a), cutoff)
      c.coeffs = d
      c.length = normalise(c, lenc)
   end
   c.val = a.val + b.val
   c.prec = prec + c.val
   renormalize!(c)
   return c
end

function add!(c::RelSeries{T}, a::RelSeries{T}) where T <: RingElement
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
         c.coeffs[i] = add!(c.coeffs[i], a.coeffs[i - vala + valc])
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
      return add!(c, b)
   elseif c === b
      return add!(c, a)
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

function (R::RelPowerSeriesRing{T})(b::RingElement) where T <: RingElement
   return R(base_ring(R)(b))
end

function (R::RelPowerSeriesRing{T})() where T <: RingElement
   z = RelSeries{T}(Vector{T}(undef, 0), 0, R.prec_max, R.prec_max)
   z.parent = R
   return z
end

function (R::RelPowerSeriesRing{T})(b::JuliaRingElement) where T <: RingElement
   bb = base_ring(R)(b)
   if is_zero(bb)
      z = RelSeries{T}(Vector{T}(undef, 0), 0, R.prec_max, R.prec_max)
   else
      z = RelSeries{T}([bb], 1, R.prec_max, 0)
   end
   z.parent = R
   return z
end

function (R::RelPowerSeriesRing{T})(b::T) where {T <: RingElem}
   parent(b) != base_ring(R) && error("Unable to coerce to power series")
   if iszero(b)
      z = RelSeries{T}(Vector{T}(undef, 0), 0, R.prec_max, R.prec_max)
   else
      z = RelSeries{T}([b], 1, R.prec_max, 0)
   end
   z.parent = R
   return z
end

function (R::RelPowerSeriesRing{T})(b::RelPowerSeriesRingElem{T}) where T <: RingElement
   parent(b) != R && error("Unable to coerce power series")
   return b
end

function (R::RelPowerSeriesRing{T})(b::Vector{T}, len::Int, prec::Int, val::Int) where T <: RingElement
   if length(b) > 0
      parent(b[1]) != base_ring(R) && error("Unable to coerce to power series")
   end
   z = RelSeries{T}(b, len, prec, val)
   z.parent = R
   return z
end

function (R::RelPowerSeriesRing{T})(b::Vector{S}, len::Int, prec::Int, val::Int) where {S <: RingElement, T <: RingElement}
   R0 = base_ring(R)
   lenb = length(b)
   entries = Vector{T}(undef, lenb)
   for i = 1:lenb
      entries[i] = R0(b[i])
   end
   z = RelSeries{T}(entries, len, prec, val)
   z.parent = R
   return z
end

###############################################################################
#
#   power_series_ring constructor
#
###############################################################################

function power_series_ring(R::AbstractAlgebra.Ring, prec::Int, s::VarName; cached::Bool=true, model::Symbol=:capped_relative)
   @req !is_trivial(R) "Zero rings are currently not supported as coefficient ring."
   T = elem_type(R)

   if model == :capped_relative
      parent_obj = RelPowerSeriesRing{T}(R, prec, Symbol(s), cached)
   elseif model == :capped_absolute
      parent_obj = AbsPowerSeriesRing{T}(R, prec, Symbol(s), cached)
   else
      error("Unknown model")
   end

   return parent_obj, gen(parent_obj)
end
