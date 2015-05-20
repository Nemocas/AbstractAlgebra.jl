###############################################################################
#
#   PowerSeries.jl : Power series over rings
#
###############################################################################    

export PowerSeries, PowerSeriesRing, O, valuation, exp, precision,
       max_precision

import Base: exp, precision

###############################################################################
#
#   Data types and memory management
#
###############################################################################

PowerSeriesID = ObjectIdDict()

type PowerSeriesRing{T <: RingElem} <: Ring
   base_ring::Ring
   prec_max::Int
   S::Symbol

   function PowerSeriesRing(R::Ring, prec::Int, s::Symbol)
      return try
         PowerSeriesID[R, prec, s]
      catch
         PowerSeriesID[R, prec, s] = new(R, prec, s)
      end
   end
end

type PowerSeries{T <: RingElem} <: PowerSeriesElem
   coeffs::Array{T, 1}
   length::Int
   prec::Int
   parent::PowerSeriesRing{T}

   PowerSeries(a::Array{T, 1}, length::Int, prec::Int) = new(a, length, prec)   
   PowerSeries(a::PowerSeries{T}) = a
end

function O{T}(a::PowerSeries{T})
   prec = length(a) - 1
   prec < 0 && throw(DomainError())
   z = PowerSeries{T}(Array(T, 0), 0, prec)
   z.parent = parent(a)
   return z
end

parent(a::PowerSeriesElem) = a.parent

elem_type{T <: RingElem}(::PowerSeriesRing{T}) = PowerSeries{T}

base_ring(R::PowerSeriesRing) = R.base_ring

base_ring(a::PowerSeriesElem) = base_ring(parent(a))

var(a::PowerSeriesRing) = a.S

function check_parent(a::PowerSeriesElem, b::PowerSeriesElem)
   parent(a) != parent(b) && 
             error("Incompatible power series rings in power series operation")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################    
   
function hash(a::PowerSeriesElem)
   h = 0xb44d6896204881f3
   for i in 0:length(a) - 1
      h $= hash(coeff(a, i))
      h = (h << 1) | (h >> (sizeof(Int)*8 - 1))
   end
   return h
end

length(x::PowerSeriesElem) = x.length

precision(x::PowerSeriesElem) = x.prec

max_precision(R::PowerSeriesRing) = R.prec_max

function normalise(a::PowerSeries, len::Int)
   while len > 0 && iszero(a.coeffs[len])
      len -= 1
   end

   return len
end

function coeff(a::PowerSeries, n::Int)
   n < 0  && throw(DomainError())
   return n >= length(a) ? zero(base_ring(a)) : a.coeffs[n + 1]
end

zero(R::PowerSeriesRing) = R(0)

one(R::PowerSeriesRing) = R(1)

function gen{T}(R::PowerSeriesRing{T})
   S = base_ring(R)
   z = PowerSeries{T}([S(0), S(1)], 2, max_precision(R) + 1)
   z.parent = R
   return z
end

iszero(a::PowerSeriesElem) = length(a) == 0

function isone(a::PowerSeriesElem)
   return length(a) == 1 && isone(coeff(a, 0))
end

function isgen(a::PowerSeriesElem)
   return length(a) == 2 && iszero(coeff(a, 0)) && isone(coeff(a, 1))
end

isunit(a::PowerSeriesElem) = isunit(coeff(a, 0))

function valuation(a::PowerSeriesElem)
   if length(a) == 0
      return precision(a)
   end
   for i = 1:length(a)
      if coeff(a, i - 1) != 0
         return i - 1
      end
   end
   error("Power series is not normalised")
end

function deepcopy{T <: RingElem}(a::PowerSeries{T})
   coeffs = Array(T, length(a))
   for i = 1:length(a)
      coeffs[i] = deepcopy(coeff(a, i - 1))
   end
   return parent(a)(coeffs)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show{T <: RingElem}(io::IO, x::PowerSeries{T})
   len = length(x)

   if len == 0
      print(io, zero(base_ring(x)))
   else
      coeff_printed = false
      c = coeff(x, 0)
      bracket = needs_parentheses(c)
      if !iszero(c)
         if bracket
            print(io, "(")
         end
         show(io, c)
         if bracket
            print(io, ")")
         end
         coeff_printed = true
      end
      for i = 1:len - 1
         c = coeff(x, i)
         bracket = needs_parentheses(c)
         if !iszero(c)
            if coeff_printed && !is_negative(c)
               print(io, "+")
            end
            if !isone(c) && (c != -1 || show_minus_one(elem_type(base_ring(x))))
               if bracket
                  print(io, "(")
               end
               show(io, c)
               if bracket
                  print(io, ")")
               end
               print(io, "*")
            end
            if c == -1 && !show_minus_one(elem_type(base_ring(x)))
               print(io, "-")
            end
            print(io, string(var(parent(x))))
            if i != 1
               print(io, "^")
               print(io, i)
            end
            coeff_printed = true
         end
      end
   end
   print(io, "+O(", string(var(parent(x))), "^", x.prec, ")")
end

function show{T <: RingElem}(io::IO, a::PowerSeriesRing{T})
   print(io, "Univariate power series ring in ", var(a), " over ")
   show(io, base_ring(a))
end

needs_parentheses(x::PowerSeriesElem) = length(x) > 1

is_negative(x::PowerSeriesElem) = length(x) <= 1 && is_negative(coeff(x, 0))

show_minus_one{T <: RingElem}(::Type{PowerSeries{T}}) = show_minus_one(T)

###############################################################################
#
#   Unary operators
#
###############################################################################

function -{T <: RingElem}(a::PowerSeries{T})
   len = length(a)
   d = Array(T, len)
   for i = 1:len
      d[i] = -coeff(a, i - 1)
   end
   return parent(a)(d, len, a.prec)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +{T <: RingElem}(a::PowerSeries{T}, b::PowerSeries{T})
   check_parent(a, b)
   lena = a.length
   lenb = b.length
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   d = Array(T, lenz)
   i = 1

   while i <= min(lena, lenb)
      d[i] = coeff(a, i - 1) + coeff(b, i - 1)
      i += 1
   end

   while i <= lena
      d[i] = coeff(a, i - 1)
      i += 1
   end

   while i <= lenb
      d[i] = coeff(b, i - 1)
      i += 1
   end

   z = parent(a)(d, i - 1, prec)

   z.length = normalise(z, i - 1)

   return z
end
  
function -{T <: RingElem}(a::PowerSeries{T}, b::PowerSeries{T})
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   
   prec = min(a.prec, b.prec)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   lenz = max(lena, lenb)
   d = Array(T, lenz)
   i = 1

   while i <= min(lena, lenb)
      d[i] = coeff(a, i - 1) - coeff(b, i - 1)
      i += 1
   end

   while i <= lena
      d[i] = coeff(a, i - 1)
      i += 1
   end

   while i <= lenb
      d[i] = -coeff(b, i - 1)
      i += 1
   end

   z = parent(a)(d, i - 1, prec)

   z.length = normalise(z, i - 1)

   return z
end

function *{T <: RingElem}(a::PowerSeries{T}, b::PowerSeries{T})
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   
   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   if lena == 0 || lenb == 0
      return parent(a)(Array(T, 0), 0, prec)
   end

   t = base_ring(a)()

   lenz = min(lena + lenb - 1, prec)

   d = Array(T, lenz)

   for i = 1:min(lena, lenz)
      d[i] = coeff(a, i - 1)*coeff(b, 0)
   end

   if lenz > lena
      for j = 2:min(lenb, lenz - lena + 1)
          d[lena + j - 1] = coeff(a, lena - 1)*coeff(b, j - 1)
      end
   end

   z = parent(a)(d, lenz, prec)

   for i = 1:lena - 1
      if lenz > i
         for j = 2:min(lenb, lenz - i + 1)
            mul!(t, coeff(a, i - 1), coeff(b, j - 1))
            addeq!(z.coeffs[i + j - 1], t)
         end
      end
   end
        
   z.length = normalise(z, lenz)

   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function fit!{T <: RingElem}(c::PowerSeries{T}, n::Int)
   if c.length < n
      t = c.coeffs
      c.coeffs = Array(T, n)
      for i = 1:c.length
         c.coeffs[i] = t[i]
      end
      for i = length(c) + 1:n
         c.coeffs[i] = zero(base_ring(c))
      end
   end
end

function setcoeff!{T <: RingElem}(c::PowerSeries{T}, n::Int, a::T)
   if (a != 0 && c.prec > n) || n + 1 <= c.length
      fit!(c, n + 1)
      c.coeffs[n + 1] = a
      c.length = max(length(c), n + 1)
      # don't normalise
   end
end

function mul!{T <: RingElem}(c::PowerSeries{T}, a::PowerSeries{T}, b::PowerSeries{T})
   lena = length(a)
   lenb = length(b)

   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   if lena == 0 || lenb == 0
      c.length = 0
   else
      t = base_ring(a)()

      lenc = min(lena + lenb - 1, prec)
      fit!(c, lenc)

      for i = 1:min(lena, lenc)
         mul!(c.coeffs[i], coeff(a, i - 1), coeff(b, 0))
      end

      if lenc > lena
         for i = 2:min(lenb, lenc - lena + 1)
            mul!(c.coeffs[lena + i - 1], coeff(a, lena - 1), coeff(b, i - 1))
         end
      end

      for i = 1:lena - 1
         if lenc > i
            for j = 2:min(lenb, lenc - i + 1)
               mul!(t, coeff(a, i - 1), coeff(b, j - 1))
               addeq!(c.coeffs[i + j - 1], t)
            end
         end
      end
        
      c.length = normalise(c, lenc)
   end
   c.prec = prec
end

function addeq!{T <: RingElem}(c::PowerSeries{T}, a::PowerSeries{T})
   lenc = length(c)
   lena = length(a)
   
   prec = min(a.prec, c.prec)
   
   lena = min(lena, prec)
   lenc = min(lenc, prec)

   len = max(lenc, lena)
   fit!(c, len)
   for i = 1:lena
      addeq!(c.coeffs[i], coeff(a, i - 1))
   end
   c.length = normalise(c, len)
   c.prec = prec
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *{T <: RingElem}(a::Int, b::PowerSeries{T})
   len = length(b)
   d = Array(T, len)
   for i = 1:len
      d[i] = a*coeff(b, i - 1)
   end
   z = parent(b)(d, len, b.prec)
   z.length = normalise(z, len)
   return z
end

function *{T <: RingElem}(a::fmpz, b::PowerSeries{T})
   len = length(b)
   d = Array(T, len)
   for i = 1:len
      d[i] = a*coeff(b, i - 1)
   end
   z = parent(b)(d, len, b.prec)
   z.length = normalise(z, len)
   return z
end

*(a::PowerSeriesElem, b::Int) = b*a

*(a::PowerSeriesElem, b::fmpz) = b*a

###############################################################################
#
#   Shifting
#
###############################################################################

function shift_left{T <: RingElem}(x::PowerSeries{T}, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   if xlen == 0
      return parent(x)(Array(T, 0), 0, x.prec + len)
   end
   v = Array(T, xlen + len)
   for i = 1:len
      v[i] = zero(base_ring(x))
   end
   for i = 1:xlen
      v[i + len] = coeff(x, i - 1)
   end
   return parent(x)(v, xlen + len, x.prec + len)
end

function shift_right{T <: RingElem}(x::PowerSeries{T}, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   if len >= xlen
      return parent(x)(Array(T, 0), 0, max(0, x.prec - len))
   end
   v = Array(T, xlen - len)
   for i = 1:xlen - len
      v[i] = coeff(x, i + len - 1)
   end
   return parent(x)(v, xlen - len, x.prec - len)
end

###############################################################################
#
#   Truncation
#
###############################################################################

function truncate{T <: RingElem}(x::PowerSeries{T}, prec::Int)
   prec < 0 && throw(DomainError())
   len = length(x)
   if x.prec <= prec
      return x
   end
   d = Array(T, prec)
   for i = 1:min(prec, len)
      d[i] = coeff(x, i - 1)
   end
   for i = len + 1:prec
      d[i] = zero(base_ring(x))
   end
   z = parent(x)(d, prec, prec)
   z.length = normalise(z, prec)
   return z
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^{T <: RingElem}(a::PowerSeries{T}, b::Int)
   b < 0 && throw(DomainError())
   # special case powers of x for constructing power series efficiently
   if isgen(a)
      d = Array(T, b + 1)
      d[b + 1] = coeff(a, 1)
      for i = 1:b
         d[i] = coeff(a, 0)
      end
      return parent(a)(d, b + 1, a.prec + b - 1)
   elseif length(a) == 0
      return parent(a)(Array(T, 0), 0, a.prec + (b - 1)*valuation(a))
   elseif length(a) == 1
      return parent(a)([a.coeffs[1]^b], 1, a.prec)
   elseif b == 0
      return parent(a)([one(base_ring(a))], 1, a.prec)
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

function =={T <: RingElem}(x::PowerSeries{T}, y::PowerSeries{T})
   check_parent(x, y)
   prec = min(x.prec, y.prec)
   
   m1 = min(length(x), length(y))
   m2 = max(length(x), length(y))
   
   m1 = min(m1, prec)
   m2 = min(m2, prec)
   if length(x) >= m2
      for i = m1 + 1: m2
         if coeff(x, i - 1) != 0
            return false
          end
      end
   else
      for i = m1 + 1: m2
         if coeff(y, i - 1) != 0
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

function isequal{T <: RingElem}(x::PowerSeries{T}, y::PowerSeries{T})
   check_parent(x, y)
   if x.prec != y.prec || length(x) != length(y)
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
#   Ad hoc comparison
#
###############################################################################

==(x::PowerSeriesElem, y::Int) = x.prec == 0 || ((length(x) == 0 && y == 0)
                                       || (length(x) == 1 && coeff(x, 0) == y))

==(x::PowerSeriesElem, y::fmpz) = x.prec == 0 || ((length(x) == 0 && y == 0)
                                       || (length(x) == 1 && coeff(x, 0) == y))

==(x::Int, y::PowerSeriesElem) = y == x

==(x::fmpz, y::PowerSeriesElem) = y == x

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact{T <: RingElem}(x::PowerSeries{T}, y::PowerSeries{T})
   check_parent(x, y)
   y == 0 && throw(DivideError())
   v2 = valuation(y)
   if v2 != 0
      v1 = valuation(x)
      if v1 >= v2
         x = shift_right(x, v2)
         y = shift_right(y, v2)
      end
   end
   y = truncate(y, x.prec)
   return x*inv(y)
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact{T <: RingElem}(x::PowerSeries{T}, y::Int)
   y == 0 && throw(DivideError())
   lenx = length(x)
   d = Array(T, lenx)
   for i = 1:lenx
      d[i] = divexact(coeff(x, i - 1), y)
   end
   return parent(x)(d, lenx, x.prec)
end

function divexact{T <: RingElem}(x::PowerSeries{T}, y::fmpz)
   y == 0 && throw(DivideError())
   lenx = length(x)
   d = Array(T, lenx)
   for i = 1:lenx
      d[i] = divexact(coeff(x, i - 1), y)
   end
   return parent(x)(d, lenx, x.prec)
end

function divexact{T <: RingElem}(x::PowerSeries{T}, y::T)
   y == 0 && throw(DivideError())
   lenx = length(x)
   d = Array(T, lenx)
   for i = 1:lenx
      d[i] = divexact(coeff(x, i - 1), y)
   end
   return parent(x)(d, lenx, x.prec)
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv{T <: RingElem}(a::PowerSeries{T})
   a == 0 && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
   a1 = coeff(a, 0)
   d = Array(T, a.prec)
   if a.prec != 0
      d[1] = divexact(one(base_ring(a)), a1)
   end
   a1 = -a1
   for n = 2:a.prec
      s = coeff(a, 1)*d[n - 1]
      for i = 2:min(n, a.length) - 1
         s += coeff(a, i)*d[n - i]
      end
      d[n] = divexact(s, a1)
   end
   ainv = parent(a)(d, a.prec, a.prec)
   ainv.length = normalise(ainv, a.prec)
   return ainv
end

###############################################################################
#
#   Special functions
#
###############################################################################

function exp{T <: RingElem}(a::PowerSeries{T})
   if a == 0
      return parent(a)([one(base_ring(a))], 1, a.prec)
   end
   d = Array(T, a.prec)
   d[0 + 1] = exp(coeff(a, 0))
   len = length(a)
   for k = 1 : a.prec - 1
      s = zero(base_ring(a))
      for j = 1 : min(k + 1, len) - 1
         s += j * coeff(a, j) * d[k - j + 1]
      end
      !isunit(base_ring(a)(k)) && error("Unable to divide in exp")
      d[k + 1] = divexact(s, k)
   end
   b = parent(a)(d, a.prec, a.prec)
   b.length = normalise(b, a.prec)
   return b
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

function Base.promote_rule{T <: RingElem, V <: Integer}(::Type{PowerSeries{T}}, ::Type{V})
   return PowerSeries{T}
end

function Base.promote_rule{T <: RingElem}(::Type{PowerSeries{T}}, ::Type{T})
   return PowerSeries{T}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function Base.call{T <: RingElem}(a::PowerSeriesRing{T}, b::RingElem)
   return a(base_ring(a)(b))
end

function Base.call{T <: RingElem}(a::PowerSeriesRing{T})
   z = PowerSeries{T}(Array(T, 0), 0, a.prec_max)
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::PowerSeriesRing{T}, b::Integer)
   if b == 0
      z = PowerSeries{T}(Array(T, 0), 0, a.prec_max)
   else
      z = PowerSeries{T}([base_ring(a)(b)], 1, a.prec_max)
   end
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::PowerSeriesRing{T}, b::T)
   parent(b) != base_ring(a) && error("Unable to coerce to power series")
   if b == 0
      z = PowerSeries{T}(Array(T, 0), 0, a.prec_max)
   else
      z = PowerSeries{T}([b], 1, a.prec_max)
   end
   z.parent = a
   return z
end

function Base.call{T <: RingElem}(a::PowerSeriesRing{T}, b::PowerSeries{T})
   parent(b) != a && error("Unable to coerce power series")
   return b
end

function Base.call{T <: RingElem}(a::PowerSeriesRing{T}, b::Array{T, 1}, len::Int, prec::Int)
   if length(b) > 0
      parent(b[1]) != base_ring(a) && error("Unable to coerce to power series")
   end
   z = PowerSeries{T}(b, len, prec)
   z.parent = a
   return z
end

###############################################################################
#
#   PowerSeriesRing constructor
#
###############################################################################

function PowerSeriesRing(R::Ring, prec::Int, s::String)
   S = symbol(s)
   T = elem_type(R)
   parent_obj = PowerSeriesRing{T}(R, prec, S)

   base = base_ring(R)
   R2 = R
   parent_type = PowerSeries{T}
   while base_ring(R2) != None
      R2 = base_ring(R2)
      T2 = elem_type(R2)
      eval(:(Base.promote_rule(::Type{$parent_type}, ::Type{$T2}) = $parent_type))
   end

   return parent_obj, parent_obj([R(0), R(1)], 2, prec)
end
