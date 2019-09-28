###############################################################################
#
#   RelSeries.jl : Power series over rings, capped relative precision
#
###############################################################################

export PowerSeriesRing, O, valuation, exp, precision, max_precision, set_prec!,
       polcoeff, set_val!, pol_length, renormalize!

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

@doc Markdown.doc"""
    O(a::AbstractAlgebra.RelSeriesElem{T}) where T <: RingElement
> Return $0 + O(x^\mbox{deg}(a))$. Usually this function is called with $x^n$
> as parameter. Then the function returns the power series $0 + O(x^n)$, which
> can be used to set the precision of a power series when constructing it.
"""
function O(a::AbstractAlgebra.RelSeriesElem{T}) where T <: RingElement
   val = pol_length(a) + valuation(a) - 1
   val < 0 && throw(DomainError())
   return parent(a)(Array{T}(undef, 0), 0, val, val)
end

parent_type(::Type{RelSeries{T}}) where T <: RingElement = RelSeriesRing{T}

@doc Markdown.doc"""
    parent(a::AbstractAlgebra.SeriesElem)
> Return the parent of the given power series.
"""
parent(a::AbstractAlgebra.SeriesElem) = a.parent

elem_type(::Type{RelSeriesRing{T}}) where T <: RingElement = RelSeries{T}

@doc Markdown.doc"""
    base_ring(R::SeriesRing{T}) where T <: RingElement
> Return the base ring of the given power series ring.
"""
base_ring(R::SeriesRing{T}) where T <: RingElement = R.base_ring::parent_type(T)

@doc Markdown.doc"""
    base_ring(a::AbstractAlgebra.SeriesElem)
> Return the base ring of the power series ring of the given power series.
"""
base_ring(a::AbstractAlgebra.SeriesElem) = base_ring(parent(a))

function isdomain_type(::Type{T}) where {S <: RingElement, T <: AbstractAlgebra.SeriesElem{S}}
   return isdomain_type(S)
end

isexact_type(a::Type{T}) where T <: AbstractAlgebra.SeriesElem = false

@doc Markdown.doc"""
    var(a::SeriesRing)
> Return the internal name of the generator of the power series ring. Note that
> this is returned as a `Symbol` not a `String`.
"""
var(a::SeriesRing) = a.S

function check_parent(a::AbstractAlgebra.SeriesElem, b::AbstractAlgebra.SeriesElem, throw::Bool = true)
   b = parent(a) != parent(b)
   b && throw && error("Incompatible power series rings in power series operation")
   return !b
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::AbstractAlgebra.RelSeriesElem, h::UInt)
   b = 0xb44d6896204881f3%UInt
   for i in 0:pol_length(a) - 1
      b = xor(b, hash(polcoeff(a, i), h), h)
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

@doc Markdown.doc"""
    pol_length(a::AbstractAlgebra.RelSeriesElem)
> Return the length of the polynomial underlying the given power series. This
> will be zero if the power series has no nonzero terms.
"""
pol_length(a::AbstractAlgebra.RelSeriesElem) = a.length

@doc Markdown.doc"""
    precision(a::AbstractAlgebra.RelSeriesElem)
> Return the precision of the given power series in absolute terms. This will
> be the sum of the valuation and the length of the underlying polynomial.
"""
precision(a::AbstractAlgebra.RelSeriesElem) = a.prec

@doc Markdown.doc"""
    valuation(a::AbstractAlgebra.RelSeriesElem)
> Return the valuation of the given power series, i.e. the degree of the first
> nonzero term (or the precision if it is arithmetically zero).
"""
valuation(a::AbstractAlgebra.RelSeriesElem) = a.val

@doc Markdown.doc"""
    max_precision(R::SeriesRing)
> Return the maximum relative precision of power series in the given power
> series ring.
"""
max_precision(R::SeriesRing) = R.prec_max

function normalise(a::RelSeries, len::Int)
   while len > 0 && iszero(a.coeffs[len])
      len -= 1
   end
   return len
end

function set_length!(a::AbstractAlgebra.SeriesElem, len::Int)
   a.length = len
end

function set_prec!(a::AbstractAlgebra.SeriesElem, prec::Int)
   a.prec = prec
end

function set_val!(a::AbstractAlgebra.RelSeriesElem, val::Int)
   a.val = val
end

function polcoeff(a::RelSeries, n::Int)
   n < 0  && throw(DomainError())
   return n >= pol_length(a) ? zero(base_ring(a)) : a.coeffs[n + 1]
end

function coeff(a::AbstractAlgebra.RelSeriesElem, n::Int)
   if n < valuation(a)
      return base_ring(a)()
   else
      return polcoeff(a, n - valuation(a))
   end
end

@doc Markdown.doc"""
    zero(R::SeriesRing)
> Return $0 + O(x^n)$ where $n$ is the maximum precision of the power series
> ring $R$.
"""
zero(R::SeriesRing) = R(0)

@doc Markdown.doc"""
    one(R::SeriesRing)
> Return $1 + O(x^n)$ where $n$ is the maximum precision of the power series
> ring $R$.
"""
one(R::SeriesRing) = R(1)

@doc Markdown.doc"""
    gen(R::RelSeriesRing)
> Return the generator of the power series ring, i.e. $x + O(x^{n + 1})$ where
> $n$ is the maximum precision of the power series ring $R$.
"""
function gen(R::RelSeriesRing)
   S = base_ring(R)
   return R([S(1)], 1, max_precision(R) + 1, 1)
end

@doc Markdown.doc"""
    iszero(a::AbstractAlgebra.RelSeriesElem)
> Return `true` if the given power series is arithmetically equal to zero to
> its current precision, otherwise return `false`.
"""
iszero(a::AbstractAlgebra.RelSeriesElem) = pol_length(a) == 0

@doc Markdown.doc"""
    isone(a::AbstractAlgebra.RelSeriesElem)
> Return `true` if the given power series is arithmetically equal to one to
> its current precision, otherwise return `false`.
"""
function isone(a::AbstractAlgebra.RelSeriesElem)
   return valuation(a) == 0 && pol_length(a) == 1 && isone(polcoeff(a, 0))
end

@doc Markdown.doc"""
    isgen(a::AbstractAlgebra.RelSeriesElem)
> Return `true` if the given power series is arithmetically equal to the
> generator of its power series ring to its current precision, otherwise return
> `false`.
"""
function isgen(a::AbstractAlgebra.RelSeriesElem)
   return valuation(a) == 1 && pol_length(a) == 1 && isone(polcoeff(a, 0))
end

@doc Markdown.doc"""
    isunit(a::AbstractAlgebra.RelSeriesElem)
> Return `true` if the given power series is arithmetically equal to a unit,
> i.e. is invertible, otherwise return `false`.
"""
isunit(a::AbstractAlgebra.RelSeriesElem) = valuation(a) == 0 && isunit(polcoeff(a, 0))

@doc Markdown.doc"""
    modulus(a::AbstractAlgebra.SeriesElem{T}) where {T <: ResElem}
> Return the modulus of the coefficients of the given power series.
"""
modulus(a::AbstractAlgebra.SeriesElem{T}) where {T <: ResElem} = modulus(base_ring(a))

function deepcopy_internal(a::RelSeries{T}, dict::IdDict) where {T <: RingElement}
   coeffs = Array{T}(undef, pol_length(a))
   for i = 1:pol_length(a)
      coeffs[i] = deepcopy(polcoeff(a, i - 1))
   end
   return parent(a)(coeffs, pol_length(a), precision(a), valuation(a))
end

function renormalize!(z::AbstractAlgebra.RelSeriesElem)
   i = 0
   zlen = pol_length(z)
   zval = valuation(z)
   zprec = precision(z)
   while i < zlen && iszero(polcoeff(z, i))
      i += 1
   end
   set_prec!(z, zprec)
   if i == zlen
      set_length!(z, 0)
      set_val!(z, zprec)
   else
      set_val!(z, zval + i)
      for j = 1:zlen - i
         z = setcoeff!(z, j - 1, polcoeff(z, j + i - 1))
      end
      set_length!(z, zlen - i)
   end
   return nothing
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function show(io::IO, x::AbstractAlgebra.RelSeriesElem)
   len = pol_length(x)
   if len == 0
      print(IOContext(io, :compact => true), zero(base_ring(x)))
   else
      coeff_printed = false
      for i = 0:len - 1
         c = polcoeff(x, i)
         bracket = needs_parentheses(c)
         if !iszero(c)
            if coeff_printed && !displayed_with_minus_in_front(c)
               print(io, "+")
            end
            if i + valuation(x) != 0
               if !isone(c) && (c != -1 || show_minus_one(elem_type(base_ring(x))))
                  if bracket
                     print(io, "(")
                  end
                  print(IOContext(io, :compact => true), c)
                  if bracket
                     print(io, ")")
                  end
                  if i + valuation(x) != 0
                     print(io, "*")
                  end
               end
               if c == -1 && !show_minus_one(elem_type(base_ring(x)))
                  print(io, "-")
               end
               print(io, string(var(parent(x))))
               if i + valuation(x) != 1
                  print(io, "^")
                  print(io, valuation(x) + i)
               end
            else
               print(IOContext(io, :compact => true), c)
            end
            coeff_printed = true
         end
      end
   end
   print(io, "+O(", string(var(parent(x))), "^", precision(x), ")")
end

function show(io::IO, a::SeriesRing)
   print(io, "Univariate power series ring in ", var(a), " over ")
   print(IOContext(io, :compact => true), base_ring(a))
end

needs_parentheses(x::AbstractAlgebra.SeriesElem) = (!iszero(x))

displayed_with_minus_in_front(x::AbstractAlgebra.SeriesElem) = false

show_minus_one(::Type{<:AbstractAlgebra.SeriesElem{T}}) where {T <: RingElement} = show_minus_one(T)

###############################################################################
#
#   Unary operators
#
###############################################################################

@doc Markdown.doc"""
    -(a::AbstractAlgebra.RelSeriesElem)
> Return $-a$.
"""
function -(a::AbstractAlgebra.RelSeriesElem)
   len = pol_length(a)
   z = parent(a)()
   set_prec!(z, precision(a))
   set_val!(z, valuation(a))
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, -polcoeff(a, i - 1))
   end
   return z
end

###############################################################################
#
#   Binary operators
#
###############################################################################

@doc Markdown.doc"""
    +(a::AbstractAlgebra.RelSeriesElem{T}, b::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
> Return $a + b$.
"""
function +(a::AbstractAlgebra.RelSeriesElem{T}, b::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
   check_parent(a, b)
   lena = pol_length(a)
   lenb = pol_length(b)
   vala = valuation(a)
   valb = valuation(b)
   valz = min(vala, valb)
   prec = min(precision(a), precision(b))
   mina = min(vala + lena, prec)
   minb = min(valb + lenb, prec)
   lenz = max(mina, minb) - valz
   R = base_ring(a)
   z = parent(a)()
   fit!(z, lenz)
   set_prec!(z, prec)
   set_val!(z, valz)
   if vala >= valb
      for i = 1:min(lenb, vala - valb)
         z = setcoeff!(z, i - 1, polcoeff(b, i - 1))
      end
      for i = lenb + 1:min(vala - valb, lenz)
         z = setcoeff!(z, i - 1, R())
      end
      for i = vala - valb + 1:lenb
         z = setcoeff!(z, i - 1, polcoeff(a, i - vala + valb - 1) + polcoeff(b, i - 1))
      end
      for i = max(lenb, vala - valb) + 1:lena + vala - valb
         z = setcoeff!(z, i - 1, polcoeff(a, i - vala + valb - 1))
      end
      for i = lena + vala - valb + 1:lenb
         z = setcoeff!(z, i - 1, polcoeff(b, i - 1))
      end
   else
      for i = 1:min(lena, valb - vala)
         z = setcoeff!(z, i - 1, polcoeff(a, i - 1))
      end
      for i = lena + 1:min(valb - vala, lenz)
         z = setcoeff!(z, i - 1, R())
      end
      for i = valb - vala + 1:lena
         z = setcoeff!(z, i - 1, polcoeff(a, i - 1) + polcoeff(b, i - valb + vala - 1))
      end
      for i = max(lena, valb - vala) + 1:lenb + valb - vala
         z = setcoeff!(z, i - 1, polcoeff(b, i - valb + vala - 1))
      end
      for i = lenb + valb - vala + 1:lena
         z = setcoeff!(z, i - 1, polcoeff(a, i - 1))
      end
   end
   set_length!(z, normalise(z, lenz))
   renormalize!(z)
   return z
end

@doc Markdown.doc"""
    -(a::AbstractAlgebra.RelSeriesElem{T}, b::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
> Return $a - b$.
"""
function -(a::AbstractAlgebra.RelSeriesElem{T}, b::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
   check_parent(a, b)
   lena = pol_length(a)
   lenb = pol_length(b)
   vala = valuation(a)
   valb = valuation(b)
   valz = min(vala, valb)
   prec = min(precision(a), precision(b))
   mina = min(vala + lena, prec)
   minb = min(valb + lenb, prec)
   lenz = max(mina, minb) - valz
   R = base_ring(a)
   z = parent(a)()
   fit!(z, lenz)
   set_prec!(z, prec)
   set_val!(z, valz)
   if vala >= valb
      for i = 1:min(lenb, vala - valb)
         z = setcoeff!(z, i - 1, -polcoeff(b, i - 1))
      end
      for i = lenb + 1:min(vala - valb, lenz)
         z = setcoeff!(z, i - 1, R())
      end
      for i = vala - valb + 1:lenb
         z = setcoeff!(z, i - 1, polcoeff(a, i - vala + valb - 1) - polcoeff(b, i - 1))
      end
      for i = max(lenb, vala - valb) + 1:lena + vala - valb
         z = setcoeff!(z, i - 1, polcoeff(a, i - vala + valb - 1))
      end
      for i = lena + vala - valb + 1:lenb
         z = setcoeff!(z, i - 1, -polcoeff(b, i - 1))
      end
   else
      for i = 1:min(lena, valb - vala)
         z = setcoeff!(z, i - 1, polcoeff(a, i - 1))
      end
      for i = lena + 1:min(valb - vala, lenz)
         z = setcoeff!(z, i - 1, R())
      end
      for i = valb - vala + 1:lena
         z = setcoeff!(z, i - 1, polcoeff(a, i - 1) - polcoeff(b, i - valb + vala - 1))
      end
      for i = max(lena, valb - vala) + 1:lenb + valb - vala
         z = setcoeff!(z, i - 1, -polcoeff(b, i - valb + vala - 1))
      end
      for i = lenb + valb - vala + 1:lena
         z = setcoeff!(z, i - 1, polcoeff(a, i - 1))
      end
   end
   set_length!(z, normalise(z, lenz))
   renormalize!(z)
   return z
end

@doc Markdown.doc"""
    *(a::AbstractAlgebra.RelSeriesElem{T}, b::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
> Return $a\times b$.
"""
function *(a::AbstractAlgebra.RelSeriesElem{T}, b::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
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
      return parent(a)(Array{T}(undef, 0), 0, prec + zval, zval)
   end
   t = base_ring(a)()
   lenz = min(lena + lenb - 1, prec)
   d = Array{T}(undef, lenz)
   for i = 1:min(lena, lenz)
      d[i] = polcoeff(a, i - 1)*polcoeff(b, 0)
   end
   if lenz > lena
      for j = 2:min(lenb, lenz - lena + 1)
          d[lena + j - 1] = polcoeff(a, lena - 1)*polcoeff(b, j - 1)
      end
   end
   for i = 1:lena - 1
      if lenz > i
         for j = 2:min(lenb, lenz - i + 1)
            t = mul!(t, polcoeff(a, i - 1), polcoeff(b, j - 1))
            d[i + j - 1] = addeq!(d[i + j - 1], t)
         end
      end
   end
   z = parent(a)(d, lenz, prec + zval, zval)
   set_length!(z, normalise(z, lenz))
   renormalize!(z)
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

@doc Markdown.doc"""
    *(a::T, b::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElem}
> Return $a\times b$.
"""
function *(a::T, b::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElem}
   len = pol_length(b)
   z = parent(b)()
   fit!(z, len)
   set_prec!(z, precision(b))
   set_val!(z, valuation(b))
   for i = 1:len
      z = setcoeff!(z, i - 1, a*polcoeff(b, i - 1))
   end
   set_length!(z, normalise(z, len))
   renormalize!(z)
   return z
end

@doc Markdown.doc"""
    *(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.RelSeriesElem)
> Return $a\times b$.
"""
function *(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.RelSeriesElem)
   len = pol_length(b)
   z = parent(b)()
   fit!(z, len)
   set_prec!(z, precision(b))
   set_val!(z, valuation(b))
   for i = 1:len
      z = setcoeff!(z, i - 1, a*polcoeff(b, i - 1))
   end
   set_length!(z, normalise(z, len))
   renormalize!(z)
   return z
end

@doc Markdown.doc"""
    *(a::AbstractAlgebra.RelSeriesElem{T}, b::T) where {T <: RingElem}
> Return $a\times b$.
"""
*(a::AbstractAlgebra.RelSeriesElem{T}, b::T) where {T <: RingElem} = b*a

@doc Markdown.doc"""
    *(a::AbstractAlgebra.RelSeriesElem, b::Union{Integer, Rational, AbstractFloat})
> Return $a\times b$.
"""
*(a::AbstractAlgebra.RelSeriesElem, b::Union{Integer, Rational, AbstractFloat}) = b*a

###############################################################################
#
#   Shifting
#
###############################################################################

@doc Markdown.doc"""
    shift_left(x::AbstractAlgebra.RelSeriesElem{T}, n::Int) where {T <: RingElement}
> Return the power series $x$ shifted left by $n$ terms, i.e. multiplied by
> $x^n$.
"""
function shift_left(x::AbstractAlgebra.RelSeriesElem{T}, n::Int) where {T <: RingElement}
   n < 0 && throw(DomainError())
   xlen = pol_length(x)
   if xlen == 0
      z = zero(parent(x))
      set_prec!(z, precision(x) + n)
      set_val!(z, valuation(x) + n)
      return z
   end
   z = parent(x)()
   fit!(z, xlen)
   set_prec!(z, precision(x) + n)
   set_val!(z, valuation(x) + n)
   for i = 1:xlen
      z = setcoeff!(z, i - 1, polcoeff(x, i - 1))
   end
   return z
end

@doc Markdown.doc"""
    shift_right(x::AbstractAlgebra.RelSeriesElem{T}, n::Int) where {T <: RingElement}
> Return the power series $x$ shifted right by $n$ terms, i.e. divided by
> $x^n$.
"""
function shift_right(x::AbstractAlgebra.RelSeriesElem{T}, n::Int) where {T <: RingElement}
   n < 0 && throw(DomainError())
   xlen = pol_length(x)
   xval = valuation(x)
   xprec = precision(x)
   z = parent(x)()
   if n >= xlen + xval
      set_prec!(z, max(0, xprec - n))
      set_val!(z, max(0, xprec - n))
   else
      zlen = min(xlen + xval - n, xlen)
      fit!(z, zlen)
      set_prec!(z, max(0, xprec - n))
      set_val!(z, max(0, xval - n))
      for i = 1:zlen
         z = setcoeff!(z, i - 1, polcoeff(x, i + xlen  - zlen - 1))
      end
      renormalize!(z)
   end
   return z
end

###############################################################################
#
#   Truncation
#
###############################################################################

@doc Markdown.doc"""
    truncate(a::AbstractAlgebra.RelSeriesElem{T}, n::Int) where {T <: RingElement}
> Return $a$ truncated to (absolute) precision $n$.
"""
function truncate(a::AbstractAlgebra.RelSeriesElem{T}, n::Int) where {T <: RingElement}
   n < 0 && throw(DomainError())
   alen = pol_length(a)
   aprec = precision(a)
   aval = valuation(a)
   if aprec <= n
      return a
   end
   z = parent(a)()
   set_prec!(z, n)
   if n <= aval
      set_length!(z, 0)
      set_val!(z, n)
   else
      fit!(z, n - aval)
      set_val!(z, aval)
      for i = 1:min(n - aval, alen)
         z = setcoeff!(z, i - 1, polcoeff(a, i - 1))
      end
      set_length!(z, normalise(z, n - aval))
   end
   return z
end

# Intended only for internal use, does not renormalize, assumes n >= 0
# Only efficient if valuation(a) == valuation(b) == 0
function mullow(a::AbstractAlgebra.RelSeriesElem{T}, b::AbstractAlgebra.RelSeriesElem{T}, n::Int) where {T <: RingElement}
   lena = pol_length(a)
   lenb = pol_length(b)
   if lena == 0 || lenb == 0
      return zero(parent(a))
   end
   prec = min(precision(a), precision(b))
   t = base_ring(a)()
   lenz = min(lena + lenb - 1, n)
   d = Array{T}(undef, lenz)
   for i = 1:min(lena, lenz)
      d[i] = coeff(a, i - 1)*coeff(b, 0)
   end
   if lenz > lena
      for j = 2:min(lenb, lenz - lena + 1)
          d[lena + j - 1] = coeff(a, lena - 1)*coeff(b, j - 1)
      end
   end
   for i = 1:lena - 1
      if lenz > i
         for j = 2:min(lenb, lenz - i + 1)
            t = mul!(t, coeff(a, i - 1), coeff(b, j - 1))
            d[i + j - 1] = addeq!(d[i + j - 1], t)
         end
      end
   end
   z = parent(a)(d, lenz, prec, 0)
   set_length!(z, normalise(z, lenz))
   return z
end

###############################################################################
#
#   Powering
#
###############################################################################

@doc Markdown.doc"""
    ^(a::AbstractAlgebra.RelSeriesElem{T}, b::Int) where {T <: RingElement}
> Return $a^b$. We require $b \geq 0$.
"""
function ^(a::AbstractAlgebra.RelSeriesElem{T}, b::Int) where {T <: RingElement}
   b < 0 && throw(DomainError())
   # special case powers of x for constructing power series efficiently
   if pol_length(a) == 0
      z = parent(a)()
      set_prec!(z, b*valuation(a))
      set_val!(z, b*valuation(a))
      return z
   elseif b == 0
      # in fact, the result would be exact 1 if we had exact series
      z = one(parent(a))
      return z
   elseif isgen(a)
      z = parent(a)()
      fit!(z, 1)
      set_prec!(z, b + precision(a) - 1)
      set_val!(z, b)
      z = setcoeff!(z, 0, deepcopy(polcoeff(a, 0)))
      set_length!(z, 1)
      return z
   elseif pol_length(a) == 1
      z = parent(a)(polcoeff(a, 0)^b)
      set_prec!(z, (b - 1)*valuation(a) + precision(a))
      set_val!(z, b*valuation(a))
      return z
   elseif b == 1
      return deepcopy(a)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      val = valuation(a)
      a = shift_right(a, val)
      prec = precision(a)
      z = a
      bit >>= 1
      while bit !=0
         z = mullow(z, z, prec)
         if (UInt(bit) & b) != 0
            z = mullow(z, a, prec)
         end
         bit >>= 1
      end
      set_val!(z, b*val)
      set_prec!(z, b*val + prec)
      renormalize!(z)
      return z
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

@doc Markdown.doc"""
    ==(x::AbstractAlgebra.RelSeriesElem{T}, y::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
> Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
> that power series to different precisions may still be arithmetically
> equal to the minimum of the two precisions.
"""
function ==(x::AbstractAlgebra.RelSeriesElem{T}, y::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
   b = check_parent(x, y, false)
   !b && return false

   xval = valuation(x)
   xprec = precision(x)
   yval = valuation(y)
   yprec = precision(y)
   prec = min(xprec, yprec)
   if prec <= xval && prec <= yval
      return true
   end
   if xval != yval
      return false
   end
   xlen = normalise(x, min(pol_length(x), prec - xval))
   ylen = normalise(y, min(pol_length(y), prec - yval))
   if xlen != ylen
      return false
   end
   for i = 1:xlen
      if polcoeff(x, i - 1) != polcoeff(y, i - 1)
         return false
      end
   end
   return true
end

@doc Markdown.doc"""
    isequal(x::AbstractAlgebra.RelSeriesElem{T}, y::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
> Return `true` if $x == y$ exactly, otherwise return `false`. Only if the
> power series are precisely the same, to the same precision, are they declared
> equal by this function.
"""
function isequal(x::AbstractAlgebra.RelSeriesElem{T}, y::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
   if parent(x) != parent(y)
      return false
   end
   if precision(x) != precision(y) || pol_length(x) != pol_length(y) ||
      valuation(x) != valuation(y)
      return false
   end
   for i = 1:pol_length(x)
      if !isequal(polcoeff(x, i - 1), polcoeff(y, i - 1))
         return false
      end
   end
   return true
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

@doc Markdown.doc"""
    ==(x::AbstractAlgebra.RelSeriesElem{T}, y::T) where {T <: RingElem}
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::AbstractAlgebra.RelSeriesElem{T}, y::T) where {T <: RingElem} = precision(x) == 0 ||
           ((pol_length(x) == 0 && iszero(y)) || (pol_length(x) == 1 &&
             valuation(x) == 0 && polcoeff(x, 0) == y))

@doc Markdown.doc"""
    ==(x::T, y::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElem}
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::T, y::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElem} = y == x

@doc Markdown.doc"""
    ==(x::AbstractAlgebra.RelSeriesElem, y::Union{Integer, Rational, AbstractFloat})
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::AbstractAlgebra.RelSeriesElem, y::Union{Integer, Rational, AbstractFloat}) = precision(x) == 0 ||
                  ((pol_length(x) == 0 && iszero(y)) || (pol_length(x) == 1 &&
                    valuation(x) == 0 && polcoeff(x, 0) == y))

@doc Markdown.doc"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.RelSeriesElem)
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.RelSeriesElem) = y == x

###############################################################################
#
#   Approximation
#
###############################################################################

function Base.isapprox(f::AbstractAlgebra.RelSeriesElem, g::AbstractAlgebra.RelSeriesElem; atol::Real=sqrt(eps()))
   check_parent(f, g)
   nmin = min(precision(f), precision(g))
   i = 1
   while i <= nmin
      if !isapprox(coeff(f, i - 1), coeff(g, i - 1); atol=atol)
         return false
      end
      i += 1
   end
   return true
end

###############################################################################
#
#   Exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact(x::AbstractAlgebra.RelSeriesElem{T}, y::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
> Return $x/y$.
"""
function divexact(x::AbstractAlgebra.RelSeriesElem{T}, y::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
   check_parent(x, y)
   iszero(y) && throw(DivideError())
   v2 = valuation(y)
   if v2 != 0
      v1 = valuation(x)
      if v1 >= v2
         x = shift_right(x, v2)
         y = shift_right(y, v2)
      else
         error("Not an exact division")
      end
   else
      x = deepcopy(x)
   end
   y = truncate(y, precision(x))
   res = parent(x)()
   set_prec!(res, min(precision(x), valuation(x) + precision(y)))
   set_val!(res, valuation(x))
   lc = coeff(y, 0)
   lc == 0 && error("Not an exact division")
   lenr = precision(x) - valuation(x)
   for i = 0:lenr - 1
      flag, q = divides(polcoeff(x, i), lc)
      !flag && error("Not an exact division")
      res = setcoeff!(res, i, q)
      for j = 0:min(precision(y) - 1, lenr - i - 1)
         x = setcoeff!(x, i + j, polcoeff(x, i + j) - polcoeff(y, j)*q)
      end
   end
   set_length!(res, normalise(res, pol_length(res)))
   return res
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact(x::AbstractAlgebra.RelSeriesElem, y::Union{Integer, Rational, AbstractFloat})
> Return $x/y$ where the quotient is expected to be exact.
"""
function divexact(x::AbstractAlgebra.RelSeriesElem, y::Union{Integer, Rational, AbstractFloat})
   y == 0 && throw(DivideError())
   lenx = pol_length(x)
   z = parent(x)()
   fit!(z, lenx)
   set_prec!(z, precision(x))
   set_val!(z, valuation(x))
   for i = 1:lenx
      z = setcoeff!(z, i - 1, divexact(polcoeff(x, i - 1), y))
   end
   return z
end

@doc Markdown.doc"""
    divexact(x::AbstractAlgebra.RelSeriesElem{T}, y::T) where {T <: RingElem}
> Return $x/y$ where the quotient is expected to be exact.
"""
function divexact(x::AbstractAlgebra.RelSeriesElem{T}, y::T) where {T <: RingElem}
   iszero(y) && throw(DivideError())
   lenx = pol_length(x)
   z = parent(x)()
   fit!(z, lenx)
   set_prec!(z, precision(x))
   set_val!(z, valuation(x))
   for i = 1:lenx
      z = setcoeff!(z, i - 1, divexact(polcoeff(x, i - 1), y))
   end
   return z
end

###############################################################################
#
#   Inversion
#
###############################################################################

@doc Markdown.doc"""
   inv(a::AbstractAlgebra.RelSeriesElem)
> Return the inverse of the power series $a$, i.e. $1/a$.
"""
function inv(a::AbstractAlgebra.RelSeriesElem)
   iszero(a) && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
   a1 = polcoeff(a, 0)
   ainv = parent(a)()
   fit!(ainv, precision(a))
   set_prec!(ainv, precision(a))
   set_val!(ainv, 0)
   if precision(a) != 0
      ainv = setcoeff!(ainv, 0, divexact(one(base_ring(a)), a1))
   end
   a1 = -a1
   for n = 2:precision(a)
      s = polcoeff(a, 1)*polcoeff(ainv, n - 2)
      for i = 2:min(n, pol_length(a)) - 1
         s += polcoeff(a, i)*polcoeff(ainv, n - i - 1)
      end
      ainv = setcoeff!(ainv, n - 1, divexact(s, a1))
   end
   set_length!(ainv, normalise(ainv, precision(a)))
   return ainv
end

###############################################################################
#
#   Square root
#
###############################################################################

@doc Markdown.doc"""
   sqrt(a::AbstractAlgebra.RelSeriesElem)
> Return the square root of the power series $a$.
"""
function Base.sqrt(a::AbstractAlgebra.RelSeriesElem)
   aval = valuation(a)
   !iseven(aval) && error("Not a square in sqrt")
   R = base_ring(a)
   !isdomain_type(elem_type(R)) && error("Sqrt not implemented over non-integral domains")
   aval2 = div(aval, 2)
   prec = precision(a) - aval
   if prec == 0
      asqrt = parent(a)()
      set_prec!(asqrt, aval2)
      set_val!(asqrt, aval2)
      return asqrt
   end
   asqrt = parent(a)()
   fit!(asqrt, prec)
   set_prec!(asqrt, prec + aval2)
   set_val!(asqrt, aval2)
   if prec > 0
      g = AbstractAlgebra.sqrt(polcoeff(a, 0))
      asqrt = setcoeff!(asqrt, 0, g)
      g2 = g + g
   end
   p = R()
   for n = 1:prec - 1
      c = R()
      for i = 1:div(n - 1, 2)
         j = n - i
         p = mul!(p, polcoeff(asqrt, i), polcoeff(asqrt, j))
         c = addeq!(c, p)
      end
      c *= 2
      if (n % 2) == 0
         i = div(n, 2)
         p = mul!(p, polcoeff(asqrt, i), polcoeff(asqrt, i))
         c = addeq!(c, p)
      end
      c = polcoeff(a, n) - c
      c = divexact(c, g2)
      asqrt = setcoeff!(asqrt, n, c)
   end
   set_length!(asqrt, normalise(asqrt, prec))
   return asqrt
end

###############################################################################
#
#   Special functions
#
###############################################################################

@doc Markdown.doc"""
    exp(a::AbstractAlgebra.RelSeriesElem)
> Return the exponential of the power series $a$.
"""
function Base.exp(a::AbstractAlgebra.RelSeriesElem)
   if iszero(a)
      z = one(parent(a))
      set_prec!(z, precision(a))
      return z
   end
   vala = valuation(a)
   preca = precision(a)
   z = parent(a)()
   R = base_ring(a)
   fit!(z, preca)
   set_prec!(z, preca)
   set_val!(z, 0)
   c = vala == 0 ? polcoeff(a, 0) : R()
   z = setcoeff!(z, 0, AbstractAlgebra.exp(c))
   len = pol_length(a) + vala
   for k = 1 : preca - 1
      s = R()
      for j = 1 : min(k + 1, len) - 1
         c = j >= vala ? polcoeff(a, j - vala) : R()
         s += j * c * polcoeff(z, k - j)
      end
      !isunit(R(k)) && error("Unable to divide in exp")
      z = setcoeff!(z, k, divexact(s, k))
   end
   set_length!(z, normalise(z, preca))
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::RelSeries)
   a.length = 0
   a.prec = parent(a).prec_max
   return a
end

function fit!(c::RelSeries{T}, n::Int) where {T <: RingElement}
   if length(c.coeffs) < n
      t = c.coeffs
      c.coeffs = Array{T}(undef, n)
      for i = 1:c.length
         c.coeffs[i] = t[i]
      end
      for i = pol_length(c) + 1:n
         c.coeffs[i] = zero(base_ring(c))
      end
   end
   return nothing
end

function setcoeff!(c::RelSeries{T}, n::Int, a::T) where {T <: RingElement}
   if (a != 0 && precision(c) - valuation(c) > n) || n + 1 <= c.length
      fit!(c, n + 1)
      c.coeffs[n + 1] = a
      c.length = max(pol_length(c), n + 1)
      # don't normalise
   end
   return c
end

function mul!(c::RelSeries{T}, a::RelSeries{T}, b::RelSeries{T}) where {T <: RingElement}
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
      t = base_ring(a)()
      lenc = min(lena + lenb - 1, prec)
      fit!(c, lenc)
      for i = 1:min(lena, lenc)
         c.coeffs[i] = mul!(c.coeffs[i], polcoeff(a, i - 1), polcoeff(b, 0))
      end
      if lenc > lena
         for i = 2:min(lenb, lenc - lena + 1)
            c.coeffs[lena + i - 1] = mul!(c.coeffs[lena + i - 1], polcoeff(a, lena - 1), polcoeff(b, i - 1))
         end
      end
      for i = 1:lena - 1
         if lenc > i
            for j = 2:min(lenb, lenc - i + 1)
               t = mul!(t, polcoeff(a, i - 1), polcoeff(b, j - 1))
               c.coeffs[i + j - 1] = addeq!(c.coeffs[i + j - 1], t)
            end
         end
      end
      c.length = normalise(c, lenc)
   end
   c.val = a.val + b.val
   c.prec = prec + c.val
   renormalize!(c)
   return c
end

function addeq!(c::RelSeries{T}, a::RelSeries{T}) where {T <: RingElement}
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
      for i = lenc + valc - vala:-1:max(lena, valc - vala) + 1
         t = c.coeffs[i]
         c.coeffs[i] = c.coeffs[i - valc + vala]
         c.coeffs[i - valc + vala] = t
      end
      for i = lena:-1:valc - vala + 1
         c.coeffs[i] = add!(c.coeffs[i], c.coeffs[i - valc + vala], a.coeffs[i])
      end
      for i = 1:min(lena, valc - vala)
         c.coeffs[i] = deepcopy(a.coeffs[i])
      end
      for i = lena + 1:min(valc - vala, lenr)
         c.coeffs[i] = R()
      end
      for i = lenc + valc - vala + 1:lena
         c.coeffs[i] = deepcopy(a.coeffs[i])
      end
   else
      for i = lenc + 1:min(vala - valc, lenr)
         c.coeffs[i] = R()
      end
      for i = vala - valc + 1:lenc
         c.coeffs[i] = addeq!(c.coeffs[i], a.coeffs[i - vala + valc])
      end
      for i = max(lenc, vala - valc) + 1:lena + vala - valc
         c.coeffs[i] = deepcopy(a.coeffs[i - vala + valc])
      end
   end
   c.length = normalise(c, lenr)
   c.prec = prec
   c.val = valr
   renormalize!(c)
   return c
end

function add!(c::RelSeries{T}, a::RelSeries{T}, b::RelSeries{T}) where {T <: RingElement}
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
      for i = 1:min(lenb, vala - valb)
         c.coeffs[i] = deepcopy(b.coeffs[i])
      end
      for i = lenb + 1:vala - valb
         c.coeffs[i] = R()
      end
      for i = vala - valb + 1:lenb
         c.coeffs[i] = add!(c.coeffs[i], a.coeffs[i - vala + valb], b.coeffs[i])
      end
      for i = max(lenb, vala - valb) + 1:lena + vala - valb
         c.coeffs[i] = deepcopy(a.coeffs[i - vala + valb])
      end
      for i = lena + vala - valb + 1:lenb
         c.coeffs[i] = deepcopy(b.coeffs[i])
      end
   else
      for i = 1:min(lena, valb - vala)
         c.coeffs[i] = deepcopy(a.coeffs[i])
      end
      for i = lena + 1:valb - vala
         c.coeffs[i] = R()
      end
      for i = valb - vala + 1:lena
         c.coeffs[i] = add!(c.coeffs[i], a.coeffs[i], b.coeffs[i - valb + vala])
      end
      for i = max(lena, valb - vala) + 1:lenb + valb - vala
         c.coeffs[i] = deepcopy(b.coeffs[i - valb + vala])
      end
      for i = lenb + valb - vala + 1:lena
         c.coeffs[i] = deepcopy(a.coeffs[i])
      end
   end
   set_length!(c, normalise(c, lenr))
   renormalize!(c)
   return c
end

###############################################################################
#
#   Random elements
#
###############################################################################

function rand(rng::AbstractRNG, S::SeriesRing, val_range::UnitRange{Int}, v...)
   R = base_ring(S)
   f = S()
   x = gen(S)
   for i = 0:S.prec_max - 1
      f += rand(rng, R, v...)*x^i
   end
   return shift_left(f, rand(rng, val_range))
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

function (R::RelSeriesRing{T})(b::RingElement) where {T <: RingElement}
   return R(base_ring(R)(b))
end

function (R::RelSeriesRing{T})() where {T <: RingElement}
   z = RelSeries{T}(Array{T}(undef, 0), 0, R.prec_max, R.prec_max)
   z.parent = R
   return z
end

function (R::RelSeriesRing{T})(b::Union{Integer, Rational, AbstractFloat}) where {T <: RingElement}
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

function (R::RelSeriesRing{T})(b::AbstractAlgebra.RelSeriesElem{T}) where {T <: RingElement}
   parent(b) != R && error("Unable to coerce power series")
   return b
end

function (R::RelSeriesRing{T})(b::Array{T, 1}, len::Int, prec::Int, val::Int) where {T <: RingElement}
   if length(b) > 0
      parent(b[1]) != base_ring(R) && error("Unable to coerce to power series")
   end
   z = RelSeries{T}(b, len, prec, val)
   z.parent = R
   return z
end

###############################################################################
#
#   PowerSeriesRing constructor
#
###############################################################################

@doc Markdown.doc"""
   PowerSeriesRing(R::AbstractAlgebra.Ring, prec::Int, s::AbstractString; cached=true, model=:capped_relative)
> Return a tuple $(S, x)$ consisting of the parent object `S` of a power series
> ring over the given base ring and a generator `x` for the power series ring.
> The maximum precision of power series in the ring is set to `prec`. If the
> model is set to `:capped_relative` this is taken as a maximum relative
> precision, and if it is set to `:capped_absolute` this is take to be a
> maximum absolute precision. The supplied string `s` specifies the way the
> generator of the power series ring will be printed. By default, the parent
> object `S` will be cached so that supplying the same base ring, string and
> precision in future will return the same parent object and generator. If
> caching of the parent object is not required, `cached` can be set to `false`.
"""
function PowerSeriesRing(R::AbstractAlgebra.Ring, prec::Int, s::AbstractString; cached=true, model=:capped_relative)
   S = Symbol(s)
   T = elem_type(R)

   if model == :capped_relative
      parent_obj = RelSeriesRing{T}(R, prec, S, cached)
   elseif model == :capped_absolute
      parent_obj = AbsSeriesRing{T}(R, prec, S, cached)
   else
      error("Unknown model")
   end

   return parent_obj, gen(parent_obj)
end
