###############################################################################
#
#   AbsSeries.jl : Power series over rings, capped relative precision
#
###############################################################################

export O, valuation, precision, max_precision, set_precision!

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

@doc Markdown.doc"""
    O(a::AbstractAlgebra.AbsSeriesElem{T}) where T <: RingElement

Return $0 + O(x^\mathrm{deg}(a))$. Usually this function is called with $x^n$
as parameter. Then the function returns the power series $0 + O(x^n)$, which
can be used to set the precision of a power series when constructing it.
"""
function O(a::AbstractAlgebra.AbsSeriesElem{T}) where T <: RingElement
   if iszero(a)
      return deepcopy(a)    # 0 + O(x^n)
   end
   prec = length(a) - 1
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

@doc Markdown.doc"""
    gen(R::AbsSeriesRing{T}) where T <: RingElement

Return the generator of the power series ring, i.e. $x + O(x^n)$ where
$n$ is the precision of the power series ring $R$.
"""
function gen(R::AbsSeriesRing{T}) where T <: RingElement
   S = base_ring(R)
   return R([S(0), S(1)], 2, max_precision(R))
end

iszero(a::SeriesElem) = length(a) == 0

function isone(a::AbsSeries)
   return (length(a) == 1 && isone(coeff(a, 0))) || precision(a) == 0
end

@doc Markdown.doc"""
    isgen(a::AbsSeries)

Return `true` if the given power series is arithmetically equal to the
generator of its power series ring to its current precision, otherwise return
`false`.
"""
function isgen(a::AbsSeries)
   return (valuation(a) == 1 && length(a) == 2 && isone(coeff(a, 1))) ||
           precision(a) == 0
end

isunit(a::AbstractAlgebra.AbsSeriesElem) = valuation(a) == 0 && isunit(coeff(a, 0))

@doc Markdown.doc"""
    valuation(a::AbstractAlgebra.AbsSeriesElem)

Return the valuation of the given power series, i.e. the degree of the first
nonzero term (or the precision if it is arithmetically zero).
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

function characteristic(a::AbsSeriesRing{T}) where T <: RingElement
   return characteristic(base_ring(a))
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function AbstractAlgebra.expressify(a::AbstractAlgebra.AbsSeriesElem,
                                    x = var(parent(a)); context = nothing)
    sum = Expr(:call, :+)
    v = valuation(a)
    len = length(a)

    for k in 0:len - 1
        c = coeff(a, k)
        if !iszero(c)
            if k == 0
                xk = 1
            elseif k == 1
                xk = x
            else
                xk = Expr(:call, :^, x, k)
            end
            if isone(c)
                push!(sum.args, Expr(:call, :*, xk))
            else
                push!(sum.args, Expr(:call, :*, expressify(c, context = context), xk))
            end
        end
    end
    push!(sum.args, Expr(:call, :O, Expr(:call, :^, x, precision(a))))
    return sum
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::AbstractAlgebra.AbsSeriesElem)
   len = length(a)
   z = parent(a)()
   z = set_precision!(z, precision(a))
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
   z = set_precision!(z, prec)
   i = 1
   while i <= min(lena, lenb)
      z = setcoeff!(z, i - 1, coeff(a, i - 1) + coeff(b, i - 1))
      i += 1
   end
   while i <= lena
      z = setcoeff!(z, i - 1, deepcopy(coeff(a, i - 1)))
      i += 1
   end
   while i <= lenb
      z = setcoeff!(z, i - 1, deepcopy(coeff(b, i - 1)))
      i += 1
   end
   z = set_length!(z, normalise(z, i - 1))
   return z
end

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
   z = set_precision!(z, prec)
   i = 1
   while i <= min(lena, lenb)
      z = setcoeff!(z, i - 1, coeff(a, i - 1) - coeff(b, i - 1))
      i += 1
   end
   while i <= lena
      z = setcoeff!(z, i - 1, deepcopy(coeff(a, i - 1)))
      i += 1
   end
   while i <= lenb
      z = setcoeff!(z, i - 1, -coeff(b, i - 1))
      i += 1
   end
   z = set_length!(z, normalise(z, i - 1))
   return z
end

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
   z = set_length!(z, normalise(z, lenz))
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::T, b::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElem}
   len = length(b)
   z = parent(b)()
   fit!(z, len)
   z = set_precision!(z, precision(b))
   for i = 1:len
      z = setcoeff!(z, i - 1, a*coeff(b, i - 1))
   end
   z = set_length!(z, normalise(z, len))
   return z
end

function *(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.AbsSeriesElem)
   len = length(b)
   z = parent(b)()
   fit!(z, len)
   z = set_precision!(z, precision(b))
   for i = 1:len
      z = setcoeff!(z, i - 1, a*coeff(b, i - 1))
   end
   z = set_length!(z, normalise(z, len))
   return z
end

*(a::AbstractAlgebra.AbsSeriesElem{T}, b::T) where {T <: RingElem} = b*a

*(a::AbstractAlgebra.AbsSeriesElem, b::Union{Integer, Rational, AbstractFloat}) = b*a

###############################################################################
#
#   Shifting
#
###############################################################################

@doc Markdown.doc"""
    shift_left(x::AbstractAlgebra.AbsSeriesElem{T}, n::Int) where {T <: RingElement}

Return the power series $x$ shifted left by $n$ terms, i.e. multiplied by
$x^n$.
"""
function shift_left(x::AbstractAlgebra.AbsSeriesElem{T}, n::Int) where {T <: RingElement}
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   xlen = length(x)
   prec = precision(x) + n
   prec = min(prec, max_precision(parent(x)))
   if xlen == 0
      z = zero(parent(x))
      z = set_precision!(z, prec)
      return z
   end
   zlen = min(prec, xlen + n)
   z = parent(x)()
   fit!(z, zlen)
   z = set_precision!(z, prec)
   for i = 1:n
      z = setcoeff!(z, i - 1, zero(base_ring(x)))
   end
   for i = 1:xlen
      z = setcoeff!(z, i + n - 1, coeff(x, i - 1))
   end
   z = set_length!(z, normalise(z, zlen))
   return z
end

@doc Markdown.doc"""
    shift_right(x::AbstractAlgebra.AbsSeriesElem{T}, n::Int) where {T <: RingElement}

Return the power series $x$ shifted right by $n$ terms, i.e. divided by
$x^n$.
"""
function shift_right(x::AbstractAlgebra.AbsSeriesElem{T}, n::Int) where {T <: RingElement}
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   xlen = length(x)
   if n >= xlen
      z = zero(parent(x))
      z = set_precision!(z, max(0, precision(x) - n))
      return z
   end
   z = parent(x)()
   fit!(z, xlen - n)
   z = set_precision!(z, precision(x) - n)
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

Return $a$ truncated to $n$ terms.
"""
function truncate(a::AbstractAlgebra.AbsSeriesElem{T}, n::Int) where {T <: RingElement}
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   len = length(a)
   if precision(a) <= n
      return a
   end
   z = parent(a)()
   fit!(z, n)
   z = set_precision!(z, n)
   for i = 1:min(n, len)
      z = setcoeff!(z, i - 1, coeff(a, i - 1))
   end
   for i = len + 1:n
      z = setcoeff!(z, i - 1, zero(base_ring(a)))
   end
   z = set_length!(z, normalise(z, n))
   return z
end

###############################################################################
#
#   Powering
#
###############################################################################

@doc Markdown.doc"""
    ^(a::AbstractAlgebra.AbsSeriesElem{T}, b::Int) where {T <: RingElement}

Return $a^b$. We require $b \geq 0$.
"""
function ^(a::AbstractAlgebra.AbsSeriesElem{T}, b::Int) where {T <: RingElement}
   b < 0 && throw(DomainError(b, "Can't take negative power"))
   # special case powers of x for constructing power series efficiently
   if b == 0
      z = one(parent(a))
      z = set_precision!(z, precision(a))
      return z
   elseif precision(a) > 0 && isgen(a) && b > 0
      # arithmetic operators must not introduce new aliasing
      return deepcopy(shift_left(a, b - 1))
   elseif length(a) == 1
      z = parent(a)(coeff(a, 0)^b)
      z = set_precision!(z, precision(a))
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

Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
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

Return `true` if $x == y$ exactly, otherwise return `false`. Only if the
power series are precisely the same, to the same precision, are they declared
equal by this function.
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

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::AbstractAlgebra.AbsSeriesElem{T}, y::T) where {T <: RingElem} = precision(x) == 0 ||
      ((length(x) == 0 && iszero(y)) || (length(x) == 1 && coeff(x, 0) == y))

@doc Markdown.doc"""
    ==(x::T, y::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElem}

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::T, y::AbstractAlgebra.AbsSeriesElem{T}) where {T <: RingElem} = y == x

@doc Markdown.doc"""
    ==(x::AbstractAlgebra.AbsSeriesElem, y::Union{Integer, Rational, AbstractFloat})

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::AbstractAlgebra.AbsSeriesElem, y::Union{Integer, Rational, AbstractFloat}) = precision(x) == 0 || ((length(x) == 0 && iszero(y))
                                       || (length(x) == 1 && coeff(x, 0) == y))

@doc Markdown.doc"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.AbsSeriesElem)

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.AbsSeriesElem) = y == x

###############################################################################
#
#   Exact division
#
###############################################################################

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
   res = set_precision!(res, min(precision(x), precision(y) + valuation(x)))
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
   res = set_length!(res, normalise(res, length(res)))
   return res
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::AbstractAlgebra.AbsSeriesElem, y::Union{Integer, Rational, AbstractFloat})
   y == 0 && throw(DivideError())
   lenx = length(x)
   z = parent(x)()
   fit!(z, lenx)
   z = set_precision!(z, precision(x))
   for i = 1:lenx
      z = setcoeff!(z, i - 1, divexact(coeff(x, i - 1), y))
   end
   return z
end

function divexact(x::AbstractAlgebra.AbsSeriesElem{T}, y::T) where {T <: RingElem}
   iszero(y) && throw(DivideError())
   lenx = length(x)
   z = parent(x)()
   fit!(z, lenx)
   z = set_precision!(z, precision(x))
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
    Base.inv(a::AbstractAlgebra.AbsSeriesElem)

Return the inverse of the power series $a$, i.e. $1/a$.
"""
function Base.inv(a::AbstractAlgebra.AbsSeriesElem)
   iszero(a) && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
   R = base_ring(a)
   a1 = coeff(a, 0)
   ainv = parent(a)()
   fit!(ainv, precision(a))
   ainv = set_precision!(ainv, precision(a))
   if precision(a) != 0
      ainv = setcoeff!(ainv, 0, divexact(one(R), a1))
   end
   a1 = -a1
   s = R()
   t = R()
   for n = 2:precision(a)
      s = mul!(s, coeff(a, 1), coeff(ainv, n - 2))
      for i = 2:min(n, length(a)) - 1
         s = addmul_delayed_reduction!(s, coeff(a, i), coeff(ainv, n - i - 1), t)
      end
      s = reduce!(s)
      ainv = setcoeff!(ainv, n - 1, divexact(s, a1))
   end
   ainv = set_length!(ainv, normalise(ainv, precision(a)))
   return ainv
end

function Base.inv(a::AbsSeriesElem{T}) where T <: FieldElement
    prec = precision(a)
    @assert valuation(a) == 0
    @assert prec != 0
    R = parent(a)
    x = R(inv(coeff(a, 0)))
    set_precision!(x, 1)
    la = [prec]
    while la[end] > 1
        push!(la, div(la[end] + 1, 2))
    end 
    two = R(2)
    two = set_precision!(two, prec)
    n = length(la) - 1
    y = R()
    minus_a = -a
    while n > 0
        # x -> x*(2 - xa) is the lifting recursion
        x = set_precision!(x, la[n])
        y = set_precision!(y, la[n])
        y = mul!(y, minus_a, x)
        y = addeq!(y, two)
        x = mul!(x, x, y)
        n -= 1 
    end
    return x
end

###############################################################################
#
#   Square root
#
###############################################################################

@doc Markdown.doc"""
    sqrt(a::AbstractAlgebra.AbsSeriesElem)

Return the square root of the power series $a$.
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
   asqrt = set_precision!(asqrt, prec)
   for n = 1:aval2
      asqrt = setcoeff!(asqrt, n - 1, R())
   end
   if prec > aval2
      g = sqrt(coeff(a, aval))
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
   asqrt = set_length!(asqrt, normalise(asqrt, prec))
   return asqrt
end

###############################################################################
#
#  Derivative and Integral
#
###############################################################################

@doc Markdown.doc"""
    derivative(f::AbsSeriesElem{T}) -> AbsSeriesElem

Return the derivative of the power series $f$.
"""
function derivative(f::AbsSeriesElem{T}) where T <: RingElement
   g = parent(f)()
   set_precision!(g, precision(f) - 1)
   len = length(f) - 1
   fit!(g, len)
   for i = 1:len
      g = setcoeff!(g, i - 1, i*coeff(f, i))
   end
   g = set_length!(g, normalise(g, len))
   return g
end

@doc Markdown.doc"""
    integral(f::AbsSeriesElem{T}) -> AbsSeriesElem

Return the integral of the power series $f$.
"""
function integral(f::AbsSeriesElem{T}) where T <: RingElement
   g = parent(f)()
   len = length(f) + 1
   fit!(g, len)
   set_precision!(g, precision(f) + 1)
   for i = 1:len - 1
      c = coeff(f, i - 1)
      if !iszero(c)
         g = setcoeff!(g, i, divexact(c, i))
      end
   end
   g = set_length!(g, normalise(g, len))
   return g
end

###############################################################################
#
#   Special functions
#
###############################################################################

@doc Markdown.doc"""
    exp(a::AbstractAlgebra.AbsSeriesElem)

Return the exponential of the power series $a$.
"""
function Base.exp(a::AbstractAlgebra.AbsSeriesElem{T}) where T <: RingElement
   if iszero(a)
      z = one(parent(a))
      z = set_precision!(z, precision(a))
      return z
   end
   z = parent(a)()
   fit!(z, precision(a))
   z = set_precision!(z, precision(a))
   z = setcoeff!(z, 0, exp(coeff(a, 0)))
   len = length(a)
   C = base_ring(a)()
   d = derivative(a)
   for k = 1 : precision(a) - 1
      s = zero(base_ring(a))
      for j = 1 : min(k + 1, len) - 1
         s = addmul_delayed_reduction!(s, coeff(d, j - 1), coeff(z, k - j), C)
      end
      s = reduce!(s)
      !isunit(base_ring(a)(k)) && error("Unable to divide in exp")
      z = setcoeff!(z, k, divexact(s, k))
   end
   z = set_length!(z, normalise(z, precision(a)))
   return z
end

function Base.exp(a::AbsSeriesElem{T}) where T <: FieldElement
   if iszero(a)
      b = parent(a)(1)
      set_precision!(b, precision(a))
      return b
   end
   R = base_ring(a)
   c = one(R)
   if valuation(a) == 0
      a = deepcopy(a)
      c = exp(coeff(a, 0))
      a = setcoeff!(a, 0, R())
   end
   x = parent(a)([R(1)], 1, min(2, precision(a)))
   prec = precision(a)
   la = [prec]
   while la[end] > 1
      push!(la, div(la[end] + 1, 2))
   end
   one1 = parent(a)([R(1)], 1, 2)
   n = length(la) - 1
   # x -> x*(1 - log(a) + a) is the recursion
   while n > 0
      set_precision!(x, la[n])
      set_precision!(one1, la[n])
      t = -log(x)
      t = addeq!(t, one1)
      t = addeq!(t, a)
      x = mul!(x, x, t)
      n -= 1 
   end
   if !isone(c)
      x *= c
   end
   return x
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
      resize!(c.coeffs, n)
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
