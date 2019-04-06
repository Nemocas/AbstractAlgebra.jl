###############################################################################
#
#   AbsSeries.jl : Power series over rings, capped relative precision
#
###############################################################################

export O, valuation, exp, precision, max_precision, set_prec!

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

@doc Markdown.doc"""
    O(a::AbstractAlgebra.AbsSeriesElem{T}) where T <: RingElement
> Returns $0 + O(x^\mbox{deg}(a))$. Usually this function is called with $x^n$
> as parameter. Then the function returns the power series $0 + O(x^n)$, which
> can be used to set the precision of a power series when constructing it.
"""
function O(a::AbstractAlgebra.AbsSeriesElem{T}) where T <: RingElement
   if iszero(a)
      return deepcopy(a)    # 0 + O(x^n)
   end
   prec = length(a) - 1
   prec < 0 && throw(DomainError())
   return parent(a)(Array{T}(undef, 0), 0, prec)
end

parent_type(::Type{AbsSeries{T}}) where {T <: RingElement} = AbsSeriesRing{T}

elem_type(::Type{AbsSeriesRing{T}}) where {T <: RingElement} = AbsSeries{T}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

length(x::AbstractAlgebra.AbsSeriesElem) = x.length

precision(x::AbstractAlgebra.AbsSeriesElem) = x.prec

@doc Markdown.doc"""
    max_precision(R::AbsSeriesRing)
> Return the maximum absolute precision of power series in the given power
> series ring.
"""
max_precision(R::AbsSeriesRing) = R.prec_max

function normalise(a::AbsSeries, len::Int)
   while len > 0 && iszero(a.coeffs[len])
      len -= 1
   end
   return len
end

function coeff(a::AbsSeries, n::Int)
   n < 0  && throw(DomainError())
   return n >= length(a) ? zero(base_ring(a)) : a.coeffs[n + 1]
end

@doc Markdown.doc"""
    gen(R::AbsSeriesRing{T}) where T <: RingElement
> Return the generator of the power series ring, i.e. $x + O(x^n)$ where
> $n$ is the precision of the power series ring $R$.
"""
function gen(R::AbsSeriesRing{T}) where T <: RingElement
   S = base_ring(R)
   return R([S(0), S(1)], 2, max_precision(R))
end

@doc Markdown.doc"""
    iszero(a::SeriesElem)
> Return `true` if the given power series is arithmetically equal to zero to
> its current precision, otherwise return `false`.
"""
iszero(a::SeriesElem) = length(a) == 0

@doc Markdown.doc"""
    isone(a::AbsSeries)
> Return `true` if the given power series is arithmetically equal to one to
> its current precision, otherwise return `false`.
"""
function isone(a::AbsSeries)
   return (length(a) == 1 && isone(coeff(a, 0))) || precision(a) == 0
end

@doc Markdown.doc"""
    isgen(a::AbsSeries)
> Return `true` if the given power series is arithmetically equal to the
> generator of its power series ring to its current precision, otherwise return
> `false`.
"""
function isgen(a::AbsSeries)
   return (valuation(a) == 1 && length(a) == 2 && isone(coeff(a, 1))) ||
           precision(a) == 0
end

@doc Markdown.doc"""
    isunit(a::AbstractAlgebra.AbsSeriesElem)
> Return `true` if the given power series is arithmetically equal to a unit,
> i.e. is invertible, otherwise return `false`.
"""
isunit(a::AbstractAlgebra.AbsSeriesElem) = valuation(a) == 0 && isunit(coeff(a, 0))

@doc Markdown.doc"""
    valuation(a::AbstractAlgebra.AbsSeriesElem)
> Return the valuation of the given power series, i.e. the degree of the first
> nonzero term (or the precision if it is arithmetically zero).
"""
function valuation(a::AbstractAlgebra.AbsSeriesElem)
   for i = 1:length(a)
      if !iszero(coeff(a, i - 1))
         return i - 1
      end
   end
   return precision(a)
end

function deepcopy_internal(a::AbsSeries{T}, dict::IdDict) where {T <: RingElement}
   coeffs = Array{T}(undef, length(a))
   for i = 1:length(a)
      coeffs[i] = deepcopy(coeff(a, i - 1))
   end
   return parent(a)(coeffs, length(a), precision(a))
end

function Base.hash(a::AbstractAlgebra.AbsSeriesElem, h::UInt)
   b = 0xb44d6896204881f3%UInt
   for i in 0:length(a) - 1
      b = xor(b, hash(coeff(a, i), h), h)
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function show(io::IO, x::AbstractAlgebra.AbsSeriesElem)
   len = length(x)

   if len == 0
      print(io, zero(base_ring(x)))
   else
      coeff_printed = false
      for i = 0:len - 1
         c = coeff(x, i)
         if !iszero(c)
            if coeff_printed
               print(io, "+")
            end
            if i != 0
               if !isone(c)
                  print(io, "(")
                  print(IOContext(io, :compact => true), c)
                  print(io, ")")
                  if i != 0
                     print(io, "*")
                  end
               end
               print(io, string(var(parent(x))))
               if i != 1
                  print(io, "^")
                  print(io, i)
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

###############################################################################
#
#   Unary operators
#
###############################################################################

@doc Markdown.doc"""
    -(a::AbstractAlgebra.AbsSeriesElem)
> Return $-a$.
"""
function -(a::AbstractAlgebra.AbsSeriesElem)
   len = length(a)
   z = parent(a)()
   set_prec!(z, precision(a))
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, -coeff(a, i - 1))
   end
   return z
end

###############################################################################
#
#   Binary operators
#
###############################################################################

@doc Markdown.doc"""
    +(a::AbstractAlgebra.AbsSeriesElem{T}, b::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
> Return $a + b$.
"""
function +(a::AbstractAlgebra.AbsSeriesElem{T}, b::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   prec = min(precision(a), precision(b))
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   lenz = max(lena, lenb)
   z = parent(a)()
   fit!(z, lenz)
   set_prec!(z, prec)
   i = 1
   while i <= min(lena, lenb)
      z = setcoeff!(z, i - 1, coeff(a, i - 1) + coeff(b, i - 1))
      i += 1
   end
   while i <= lena
      z = setcoeff!(z, i - 1, coeff(a, i - 1))
      i += 1
   end
   while i <= lenb
      z = setcoeff!(z, i - 1, coeff(b, i - 1))
      i += 1
   end
   set_length!(z, normalise(z, i - 1))
   return z
end

@doc Markdown.doc"""
    -(a::AbstractAlgebra.AbsSeriesElem{T}, b::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
> Return $a - b$.
"""
function -(a::AbstractAlgebra.AbsSeriesElem{T}, b::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   prec = min(precision(a), precision(b))
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   lenz = max(lena, lenb)
   z = parent(a)()
   fit!(z, lenz)
   set_prec!(z, prec)
   i = 1
   while i <= min(lena, lenb)
      z = setcoeff!(z, i - 1, coeff(a, i - 1) - coeff(b, i - 1))
      i += 1
   end
   while i <= lena
      z = setcoeff!(z, i - 1, coeff(a, i - 1))
      i += 1
   end
   while i <= lenb
      z = setcoeff!(z, i - 1, -coeff(b, i - 1))
      i += 1
   end
   set_length!(z, normalise(z, i - 1))
   return z
end

@doc Markdown.doc"""
    *(a::AbstractAlgebra.AbsSeriesElem{T}, b::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
> Return $a\times b$.
"""
function *(a::AbstractAlgebra.AbsSeriesElem{T}, b::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
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
      return parent(a)(Array{T}(undef, 0), 0, prec)
   end
   t = base_ring(a)()
   lenz = min(lena + lenb - 1, prec)
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
   z = parent(a)(d, lenz, prec)
   set_length!(z, normalise(z, lenz))
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

@doc Markdown.doc"""
    *(a::T, b::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElem}
> Return $a\times b$.
"""
function *(a::T, b::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElem}
   len = length(b)
   z = parent(b)()
   fit!(z, len)
   set_prec!(z, precision(b))
   for i = 1:len
      z = setcoeff!(z, i - 1, a*coeff(b, i - 1))
   end
   set_length!(z, normalise(z, len))
   return z
end

@doc Markdown.doc"""
    *(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.AbsSeriesElem)
> Return $a\times b$.
"""
function *(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.AbsSeriesElem)
   len = length(b)
   z = parent(b)()
   fit!(z, len)
   set_prec!(z, precision(b))
   for i = 1:len
      z = setcoeff!(z, i - 1, a*coeff(b, i - 1))
   end
   set_length!(z, normalise(z, len))
   return z
end

@doc Markdown.doc"""
    *(a::AbstractAlgebra.AbsSeriesElem{T}, b::T) where {T <: RingElem}
> Return $a\times b$.
"""
*(a::AbstractAlgebra.AbsSeriesElem{T}, b::T) where {T <: RingElem} = b*a

@doc Markdown.doc"""
    *(a::AbstractAlgebra.AbsSeriesElem, b::Union{Integer, Rational, AbstractFloat})
> Return $a\times b$.
"""
*(a::AbstractAlgebra.AbsSeriesElem, b::Union{Integer, Rational, AbstractFloat}) = b*a

###############################################################################
#
#   Shifting
#
###############################################################################

@doc Markdown.doc"""
    shift_left(x::AbstractAlgebra.AbsSeriesElem{T}, n::Int) where {T <: RingElement}
> Return the power series $x$ shifted left by $n$ terms, i.e. multiplied by
> $x^n$.
"""
function shift_left(x::AbstractAlgebra.AbsSeriesElem{T}, n::Int) where {T <: RingElement}
   n < 0 && throw(DomainError())
   xlen = length(x)
   prec = precision(x) + n
   prec = min(prec, max_precision(parent(x)))
   if xlen == 0
      z = zero(parent(x))
      set_prec!(z, prec)
      return z
   end
   zlen = min(prec, xlen + n)
   z = parent(x)()
   fit!(z, zlen)
   set_prec!(z, prec)
   for i = 1:n
      z = setcoeff!(z, i - 1, zero(base_ring(x)))
   end
   for i = 1:xlen
      z = setcoeff!(z, i + n - 1, coeff(x, i - 1))
   end
   set_length!(z, normalise(z, zlen))
   return z
end

@doc Markdown.doc"""
    shift_right(x::AbstractAlgebra.AbsSeriesElem{T}, n::Int) where {T <: RingElement}
> Return the power series $x$ shifted right by $n$ terms, i.e. divided by
> $x^n$.
"""
function shift_right(x::AbstractAlgebra.AbsSeriesElem{T}, n::Int) where {T <: RingElement}
   n < 0 && throw(DomainError())
   xlen = length(x)
   if n >= xlen
      z = zero(parent(x))
      set_prec!(z, max(0, precision(x) - n))
      return z
   end
   z = parent(x)()
   fit!(z, xlen - n)
   set_prec!(z, precision(x) - n)
   for i = 1:xlen - n
      z = setcoeff!(z, i - 1, coeff(x, i + n - 1))
   end
   return z
end

###############################################################################
#
#   Truncation
#
###############################################################################

@doc Markdown.doc"""
    truncate(a::AbstractAlgebra.AbsSeriesElem{T}, n::Int) where {T <: RingElement}
> Return $a$ truncated to $n$ terms.
"""
function truncate(a::AbstractAlgebra.AbsSeriesElem{T}, n::Int) where {T <: RingElement}
   n < 0 && throw(DomainError())
   len = length(a)
   if precision(a) <= n
      return a
   end
   z = parent(a)()
   fit!(z, n)
   set_prec!(z, n)
   for i = 1:min(n, len)
      z = setcoeff!(z, i - 1, coeff(a, i - 1))
   end
   for i = len + 1:n
      z = setcoeff!(z, i - 1, zero(base_ring(a)))
   end
   set_length!(z, normalise(z, n))
   return z
end

###############################################################################
#
#   Powering
#
###############################################################################

@doc Markdown.doc"""
    ^(a::AbstractAlgebra.AbsSeriesElem{T}, b::Int) where {T <: RingElement}
> Return $a^b$. We require $b \geq 0$.
"""
function ^(a::AbstractAlgebra.AbsSeriesElem{T}, b::Int) where {T <: RingElement}
   b < 0 && throw(DomainError())
   # special case powers of x for constructing power series efficiently
   if b == 0
      z = one(parent(a))
      set_prec!(z, precision(a))
      return z
   elseif precision(a) > 0 && isgen(a) && b > 0
      # arithmetic operators must not introduce new aliasing
      return deepcopy(shift_left(a, b - 1))
   elseif length(a) == 1
      z = parent(a)(coeff(a, 0)^b)
      set_prec!(z, precision(a))
      return z
   elseif b == 1
      return deepcopy(a)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit !=0
         z = z*z
         if (UInt(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

@doc Markdown.doc"""
    ==(x::AbstractAlgebra.AbsSeriesElem{T}, y::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
> Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
> that power series to different precisions may still be arithmetically
> equal to the minimum of the two precisions.
"""
function ==(x::AbstractAlgebra.AbsSeriesElem{T}, y::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
   b = check_parent(x, y, false)
   !b && return false

   prec = min(precision(x), precision(y))
   m1 = min(length(x), length(y))
   m2 = max(length(x), length(y))
   m1 = min(m1, prec)
   m2 = min(m2, prec)
   if length(x) >= m2
      for i = m1 + 1: m2
         if !iszero(coeff(x, i - 1))
            return false
          end
      end
   else
      for i = m1 + 1: m2
         if !iszero(coeff(y, i - 1))
            return false
          end
      end
   end
   for i = 1:m1
      if coeff(x, i - 1) != coeff(y, i - 1)
         return false
      end
   end
   return true
end

@doc Markdown.doc"""
    isequal(x::AbstractAlgebra.AbsSeriesElem{T}, y::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
> Return `true` if $x == y$ exactly, otherwise return `false`. Only if the
> power series are precisely the same, to the same precision, are they declared
> equal by this function.
"""
function isequal(x::AbstractAlgebra.AbsSeriesElem{T}, y::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
   if parent(x) != parent(y)
      return false
   end
   if precision(x) != precision(y) || length(x) != length(y)
      return false
   end
   for i = 1:length(x)
      if !isequal(coeff(x, i - 1), coeff(y, i - 1))
         return false
      end
   end
   return true
end

###############################################################################
#
#   Approximation
#
###############################################################################

function Base.isapprox(f::AbstractAlgebra.AbsSeriesElem, g::AbstractAlgebra.AbsSeriesElem; atol::Real=sqrt(eps()))
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
#   Ad hoc comparison
#
###############################################################################

@doc Markdown.doc"""
    ==(x::AbstractAlgebra.AbsSeriesElem{T}, y::T) where {T <: RingElem}
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::AbstractAlgebra.AbsSeriesElem{T}, y::T) where {T <: RingElem} = precision(x) == 0 ||
      ((length(x) == 0 && iszero(y)) || (length(x) == 1 && coeff(x, 0) == y))

@doc Markdown.doc"""
    ==(x::T, y::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElem}
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::T, y::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElem} = y == x

@doc Markdown.doc"""
    ==(x::AbstractAlgebra.AbsSeriesElem, y::Union{Integer, Rational, AbstractFloat})
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::AbstractAlgebra.AbsSeriesElem, y::Union{Integer, Rational, AbstractFloat}) = precision(x) == 0 || ((length(x) == 0 && iszero(y))
                                       || (length(x) == 1 && coeff(x, 0) == y))

@doc Markdown.doc"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.AbsSeriesElem)
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.AbsSeriesElem) = y == x

###############################################################################
#
#   Exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact(x::AbstractAlgebra.AbsSeriesElem{T}, y::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
> Return $x/y$.
"""
function divexact(x::AbstractAlgebra.AbsSeriesElem{T}, y::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
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
   set_prec!(res, min(precision(x), precision(y) + valuation(x)))
   lc = coeff(y, 0)
   lc == 0 && error("Not an exact division")
   lenr = precision(x)
   for i = valuation(x):lenr - 1
      flag, q = divides(coeff(x, i), lc)
      !flag && error("Not an exact division")
      res = setcoeff!(res, i, q)
      for j = 0:min(precision(y) - 1, lenr - i - 1)
         x = setcoeff!(x, i + j, coeff(x, i + j) - coeff(y, j)*q)
      end
   end
   set_length!(res, normalise(res, length(res)))
   return res
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact(x::AbstractAlgebra.AbsSeriesElem, y::Union{Integer, Rational, AbstractFloat})
> Return $x/y$ where the quotient is expected to be exact.
"""
function divexact(x::AbstractAlgebra.AbsSeriesElem, y::Union{Integer, Rational, AbstractFloat})
   y == 0 && throw(DivideError())
   lenx = length(x)
   z = parent(x)()
   fit!(z, lenx)
   set_prec!(z, precision(x))
   for i = 1:lenx
      z = setcoeff!(z, i - 1, divexact(coeff(x, i - 1), y))
   end
   return z
end

@doc Markdown.doc"""
    divexact(x::AbstractAlgebra.AbsSeriesElem{T}, y::T) where {T <: RingElem}
> Return $x/y$ where the quotient is expected to be exact.
"""
function divexact(x::AbstractAlgebra.AbsSeriesElem{T}, y::T) where {T <: RingElem}
   iszero(y) && throw(DivideError())
   lenx = length(x)
   z = parent(x)()
   fit!(z, lenx)
   set_prec!(z, precision(x))
   for i = 1:lenx
      z = setcoeff!(z, i - 1, divexact(coeff(x, i - 1), y))
   end
   return z
end

###############################################################################
#
#   Inversion
#
###############################################################################

@doc Markdown.doc"""
   inv(a::AbstractAlgebra.AbsSeriesElem)
> Return the inverse of the power series $a$, i.e. $1/a$.
"""
function inv(a::AbstractAlgebra.AbsSeriesElem)
   iszero(a) && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
   a1 = coeff(a, 0)
   ainv = parent(a)()
   fit!(ainv, precision(a))
   set_prec!(ainv, precision(a))
   if precision(a) != 0
      ainv = setcoeff!(ainv, 0, divexact(one(base_ring(a)), a1))
   end
   a1 = -a1
   for n = 2:precision(a)
      s = coeff(a, 1)*coeff(ainv, n - 2)
      for i = 2:min(n, length(a)) - 1
         s += coeff(a, i)*coeff(ainv, n - i - 1)
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
   sqrt(a::AbstractAlgebra.AbsSeriesElem)
> Return the square root of the power series $a$.
"""
function Base.sqrt(a::AbstractAlgebra.AbsSeriesElem)
   # Given a power series f = f0 + f1*x + f2*x^2 + ..., compute the square root
   # g = g0 + g1*x + g2*x^2 + ... using the relations g0^2 = f0, 2g0*g1 = f1
   # 2g0*g2 = f2 - g1^2, 2g0*g3 = f3 - 2g1*g2, 2g0*g4 = f4 - (2g1*g3 + g2^2), etc.
   # where the terms being subtracted are those contributing to the i-th
   # coefficient of the square of g
   aval = valuation(a)
   !iseven(aval) && error("Not a square in sqrt")
   R = base_ring(a)
   !isdomain_type(elem_type(R)) && error("Sqrt not implemented over non-integral domains")
   if iszero(a)
      return deepcopy(a)
   end
   aval2 = div(aval, 2)
   prec = precision(a) - aval2
   asqrt = parent(a)()
   fit!(asqrt, prec)
   set_prec!(asqrt, prec)
   for n = 1:aval2
      asqrt = setcoeff!(asqrt, n - 1, R())
   end
   if prec > aval2
      g = AbstractAlgebra.sqrt(coeff(a, aval))
      setcoeff!(asqrt, aval2, g)
      g2 = g + g
   end
   p = R()
   for n = 1:prec - aval2 - 1
      c = R()
      for i = 1:div(n - 1, 2)
         j = n - i
         p = mul!(p, coeff(asqrt, aval2 + i), coeff(asqrt, aval2 + j))
         c = addeq!(c, p)
      end
      c *= 2
      if (n % 2) == 0
         i = div(n, 2)
         p = mul!(p, coeff(asqrt, aval2 + i), coeff(asqrt, aval2 + i))
         c = addeq!(c, p)
      end
      c = coeff(a, n + aval) - c
      c = divexact(c, g2)
      asqrt = setcoeff!(asqrt, aval2 + n, c)
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
    exp(a::AbstractAlgebra.AbsSeriesElem)
> Return the exponential of the power series $a$.
"""
function Base.exp(a::AbstractAlgebra.AbsSeriesElem)
   if iszero(a)
      z = one(parent(a))
      set_prec!(z, precision(a))
      return z
   end
   z = parent(a)()
   fit!(z, precision(a))
   set_prec!(z, precision(a))
   z = setcoeff!(z, 0, AbstractAlgebra.exp(coeff(a, 0)))
   len = length(a)
   for k = 1 : precision(a) - 1
      s = zero(base_ring(a))
      for j = 1 : min(k + 1, len) - 1
         s += j * coeff(a, j) * coeff(z, k - j)
      end
      !isunit(base_ring(a)(k)) && error("Unable to divide in exp")
      z = setcoeff!(z, k, divexact(s, k))
   end
   set_length!(z, normalise(z, precision(a)))
   return z
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

function fit!(c::AbsSeries{T}, n::Int) where {T <: RingElement}
   if length(c.coeffs) < n
      t = c.coeffs
      c.coeffs = Array{T}(undef, n)
      for i = 1:c.length
         c.coeffs[i] = t[i]
      end
      for i = length(c) + 1:n
         c.coeffs[i] = zero(base_ring(c))
      end
   end
   return nothing
end

function setcoeff!(c::AbsSeries{T}, n::Int, a::T) where {T <: RingElement}
   if (!iszero(a) && precision(c) > n) || n + 1 <= c.length
      fit!(c, n + 1)
      c.coeffs[n + 1] = a
      c.length = max(length(c), n + 1)
      # don't normalise
   end
   return c
end

function mul!(c::AbsSeries{T}, a::AbsSeries{T}, b::AbsSeries{T}) where {T <: RingElement}
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
      t = base_ring(a)()

      lenc = min(lena + lenb - 1, prec)
      fit!(c, lenc)

      for i = 1:min(lena, lenc)
         c.coeffs[i] = mul!(c.coeffs[i], coeff(a, i - 1), coeff(b, 0))
      end

      if lenc > lena
         for i = 2:min(lenb, lenc - lena + 1)
            c.coeffs[lena + i - 1] = mul!(c.coeffs[lena + i - 1], coeff(a, lena - 1), coeff(b, i - 1))
         end
      end

      for i = 1:lena - 1
         if lenc > i
            for j = 2:min(lenb, lenc - i + 1)
               t = mul!(t, coeff(a, i - 1), coeff(b, j - 1))
               c.coeffs[i + j - 1] = addeq!(c.coeffs[i + j - 1], t)
            end
         end
      end

      c.length = normalise(c, lenc)
   end
   c.prec = prec
   return c
end

function addeq!(c::AbsSeries{T}, a::AbsSeries{T}) where {T <: RingElement}
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

function (a::AbsSeriesRing{T} where {T <: RingElement})(b::RingElement)
   return a(base_ring(a)(b))
end

function (a::AbsSeriesRing{T})() where {T <: RingElement}
   z = AbsSeries{T}(Array{T}(undef, 0), 0, a.prec_max)
   z.parent = a
   return z
end

function (a::AbsSeriesRing{T})(b::Union{Integer, Rational, AbstractFloat}) where {T <: RingElement}
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

function (a::AbsSeriesRing{T})(b::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElement}
   parent(b) != a && error("Unable to coerce power series")
   return b
end

function (a::AbsSeriesRing{T})(b::Array{T, 1}, len::Int, prec::Int) where {T <: RingElement}
   if length(b) > 0
      parent(b[1]) != base_ring(a) && error("Unable to coerce to power series")
   end
   z = AbsSeries{T}(b, len, prec)
   z.parent = a
   return z
end

###############################################################################
#
#   PowerSeriesRing constructor
#
###############################################################################

# see RelSeries.jl
