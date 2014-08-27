###########################################################################################
#
#   PowerSeries.jl : Power series over rings
#
###########################################################################################    

export PowerSeries, PowerSeriesRing, O, valuation, min, max, isless, Precision, initps, exp

import Base: min, max, isless, exp

###########################################################################################
#
#   Precision type
#
###########################################################################################

typealias Precision Union(Nothing, Int)

max(a::Nothing, b::Int) = nothing

max(a::Int, b::Nothing) = nothing

max(a::Nothing, b::Nothing) = nothing

min(a::Nothing, b::Int) = b

min(a::Int, b::Nothing) = a

min(a::Nothing, b::Nothing) = nothing

isless(a::Nothing, b::Int) = false

isless(a::Int, b::Nothing) = true

isless(a::Nothing, b::Nothing) = false

+(a::Nothing, b::Int) = nothing

+(a::Int, b::Nothing) = nothing

+(a::Nothing, b::Nothing) = nothing

-(a::Nothing, b::Int) = nothing

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

type initps end

type PowerSeries{T <: Ring, S} <: Ring
   coeffs::Ptr{Void}
   len::Int
   alloc::Int
   inv::Int
   prec :: Precision
   data :: PolyStruct{T}
   
   PowerSeries(a :: PolyStruct{T}, n :: Precision) = new(C_NULL, 0, 0, 0, n, a)   

   PowerSeries(a :: initps, n :: Precision) = new(C_NULL, 0, 0, 0, n)
   
   PowerSeries() = PowerSeries(PowerSeries{T, S}, Array(T, 0), nothing)
   
   PowerSeries(a::Integer) = a == 0 ? PowerSeries(PowerSeries{T, S}, Array(T, 0), nothing) : PowerSeries(PowerSeries{T, S}, [T(a)], nothing)
   PowerSeries(a::T) = PowerSeries(PowerSeries{T, S}, [a], nothing)
   PowerSeries(a::PowerSeries{T, S}) = a
   PowerSeries{R <: Ring}(a::R) = convert(PowerSeries{T, S}, a)
end

function PowerSeries{T, S}(::Type{PowerSeries{T, S}}, a :: Array{T, 1}, n :: Precision)
   len = length(a)
   d = PolyStruct(a, len)
   z = PowerSeries{T, S}(d, n)
   z.data.length = normalise(z, len)
   return z
end

function O{T, S}(a :: PowerSeries{T, S})
   prec = length(a) - 1
   prec < 0 && throw(DivideError())
   a.prec != nothing && error("Invalid power series monomial in O()")
   return PowerSeries(PowerSeries{T, S}, Array(T, 0), prec)
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################    
   
length{T <: Ring, S}(x::PowerSeries{T, S}) = x.data.length

function normalise{T <: Ring, S}(a::PowerSeries{T, S}, len::Int)
   while len > 0 && a.data.coeffs[len] == 0 # cannot use coeff(a, len - 1) here
      len -= 1
   end

   return len
end

coeff{T <: Ring, S}(a::PowerSeries{T, S}, n::Int) = n < 0 || n >= a.data.length ? T(0) : a.data.coeffs[n + 1]

zero{T <: Ring, S}(::Type{PowerSeries{T, S}}) = PowerSeries{T, S}(0)

one{T <: Ring, S}(::Type{PowerSeries{T, S}}) = PowerSeries{T, S}(1)

gen{T <: Ring, S}(::Type{PowerSeries{T, S}}) = PowerSeries(PowerSeries{T, S}, [T(0), T(1)], nothing)

isgen{T <: Ring, S}(a::PowerSeries{T, S}) = a.prec == nothing && a.data.length == 2 && coeff(a, 0) == 0 && coeff(a, 1) == 1

isunit{T <: Ring, S}(a::PowerSeries{T, S}) = isunit(coeff(a, 0))

function valuation{T <: Ring, S}(a::PowerSeries{T, S})
   if length(a) == 0
      return a.prec
   end
   for i = 1:length(a)
      if coeff(a, i - 1) != 0
         return i - 1
      end
   end
   error("Power series is not normalised")
end

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{T <: Ring, S}(io::IO, x::PowerSeries{T, S})
   len = x.data.length

   if len == 0
      print(io, zero(T))
   else
      coeff_printed = false
      c = x.data.coeffs[1]
      bracket = needs_parentheses(c)
      if c != 0
         if bracket
            print(io, "(")
         end
         show(io, c)
         if bracket
            print(io, ")")
         end
         coeff_printed = true
      end
      for i = 2:len
         c = x.data.coeffs[i]
         bracket = needs_parentheses(c)
         if c != 0
            if coeff_printed && !is_negative(c)
               print(io, "+")
            end
            if c != 1 && (c != -1 || show_minus_one(typeof(c)))
               if bracket
                  print(io, "(")
               end
               show(io, c)
               if bracket
                  print(io, ")")
               end
               print(io, "*")
            end
            if c == -1 && !show_minus_one(typeof(c))
               print(io, "-")
            end
            print(io, string(S))
            if i != 2
               print(io, "^")
               print(io, i - 1)
            end
            coeff_printed = true
         end
      end
      
   end
   if x.prec != nothing
      print(io, "+O(", string(S), "^", x.prec, ")")
   end
end

function show{T <: Ring, S}(io::IO, ::Type{PowerSeries{T, S}})
   print(io, "Univariate power series ring in ", string(S), " over ")
   show(io, T)
end

needs_parentheses{T <: Ring, S}(x::PowerSeries{T, S}) = length(s) > 1

is_negative{T <: Ring, S}(x::PowerSeries{T, S}) = length(x) <= 1 && is_negative(coeff(x, 0))

show_minus_one{T <: Ring, S}(::Type{PowerSeries{T, S}}) = show_minus_one(T)

###########################################################################################
#
#   Unary operators
#
###########################################################################################

function -{T <: Ring, S}(a::PowerSeries{T, S})
   len = a.data.length
   d = Array(T, len)
   for i = 1:len
      d[i] = -a.data.coeffs[i]
   end
   z = PowerSeries(PowerSeries{T, S}, d, a.prec)
   z.data.length = len
   return z
end

###########################################################################################
#
#   Binary operators
#
###########################################################################################

function +{T <: Ring, S}(a::PowerSeries{T, S}, b::PowerSeries{T, S})
   lena = a.data.length
   lenb = b.data.length
         
   prec = min(a.prec, b.prec)
 
   lena = min(lena, prec)
   lenb = min(lenb, prec)

   lenz = max(lena, lenb)
   d = Array(T, lenz)
   i = 1

   while i <= min(lena, lenb)
      d[i] = a.data.coeffs[i] + b.data.coeffs[i]
      i += 1
   end

   while i <= lena
      d[i] = a.data.coeffs[i]
      i += 1
   end

   while i <= lenb
      d[i] = b.data.coeffs[i]
      i += 1
   end

   z = PowerSeries(PowerSeries{T, S}, d, prec)

   z.data.length = normalise(z, i - 1)

   return z
end
  
function -{T <: Ring, S}(a::PowerSeries{T, S}, b::PowerSeries{T, S})
   lena = a.data.length
   lenb = b.data.length
   
   prec = min(a.prec, b.prec)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   lenz = max(lena, lenb)
   d = Array(T, lenz)
   i = 1

   while i <= min(lena, lenb)
      d[i] = a.data.coeffs[i] - b.data.coeffs[i]
      i += 1
   end

   while i <= lena
      d[i] = a.data.coeffs[i]
      i += 1
   end

   while i <= lenb
      d[i] = -b.data.coeffs[i]
      i += 1
   end

   z = PowerSeries(PowerSeries{T, S}, d, prec)

   z.data.length = normalise(z, i - 1)

   return z
end

function *{T <: Ring, S}(a::PowerSeries{T, S}, b::PowerSeries{T, S})
   lena = a.data.length
   lenb = b.data.length
   
   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   if lena == 0 || lenb == 0
      return PowerSeries(PowerSeries{T, S}, Array(T, 0), prec)
   end

   t = T()

   lenz = prec == nothing ? lena + lenb - 1 : min(lena + lenb - 1, prec)

   d = Array(T, lenz)

   for i = 1:min(lena, lenz)
      d[i] = a.data.coeffs[i]*b.data.coeffs[1]
   end

   if lenz > lena
      for j = 2:min(lenb, lenz - lena + 1)
          d[lena + j - 1] = a.data.coeffs[lena]*b.data.coeffs[j]
      end
   end

   z = PowerSeries(PowerSeries{T, S}, d, prec)

   for i = 1:lena - 1
      if lenz > i
         for j = 2:min(lenb, lenz - i + 1)
            mul!(t, a.data.coeffs[i], b.data.coeffs[j])
            addeq!(z.data.coeffs[i + j - 1], t)
         end
      end
   end
        
   z.data.length = normalise(z, lenz)

   return z
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function fit!{T <: Ring, S}(c::PowerSeries{T, S}, n::Int)
   if c.data.length < n
      t = c.data.coeffs
      c.data.coeffs = Array(T, n)
      for i = 1:c.data.length
         c.data.coeffs[i] = t[i]
      end
      for i = c.data.length + 1:n
         c.data.coeffs[i] = zero(T)
      end
   end
end

function setcoeff!{T <: Ring, S}(c::PowerSeries{T, S}, n::Int, a::T)
   if (a != 0 && (c.prec == nothing || c.prec > n)) || n + 1 <= c.data.length
      fit!(c, n + 1)
      c.data.coeffs[n + 1] = a
      c.data.length = max(c.data.length, n + 1)
      # don't normalise
   end
end

function mul!{T <: Ring, S}(c::PowerSeries{T, S}, a::PowerSeries{T, S}, b::PowerSeries{T, S})
   lena = a.data.length
   lenb = b.data.length

   aval = valuation(a)
   bval = valuation(b)

   prec = min(a.prec + bval, b.prec + aval)
   
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   
   if lena == 0 || lenb == 0
      c.data.length = 0
   else
      t = T()

      lenc = prec == nothing ? lena + lenb - 1 : min(lena + lenb - 1, prec)
      fit!(c, lenc)

      for i = 1:min(lena, lenc)
         mul!(c.data.coeffs[i], a.data.coeffs[i], b.data.coeffs[1])
      end

      if lenc > lena
         for i = 2:min(lenb, lenc - lena + 1)
            mul!(c.data.coeffs[lena + i - 1], a.data.coeffs[lena], b.data.coeffs[i])
         end
      end

      for i = 1:lena - 1
         if lenc > i
            for j = 2:min(lenb, lenc - i + 1)
               mul!(t, a.data.coeffs[i], b.data.coeffs[j])
               addeq!(c.data.coeffs[i + j - 1], t)
            end
         end
      end
        
      c.data.length = normalise(c, lenc)
   end
   c.prec = prec
end

function addeq!{T <: Ring, S}(c::PowerSeries{T, S}, a::PowerSeries{T, S})
   lenc = c.data.length
   lena = a.data.length
   
   prec = min(a.prec, c.prec)
   
   lena = min(lena, prec)
   lenc = min(lenc, prec)

   len = max(lenc, lena)
   fit!(c, len)
   for i = 1:lena
      addeq!(c.data.coeffs[i], a.data.coeffs[i])
   end
   c.data.length = normalise(c, len)
   c.prec = prec
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *{T <: Ring, S}(a::Int, b::PowerSeries{T, S})
   len = b.data.length
   d = Array(T, len)
   for i = 1:len
      d[i] = a*coeff(b, i - 1)
   end
   z = PowerSeries(PowerSeries{T, S}, d, b.prec)
   z.data.length = normalise(z, len)
   return z
end

function *{T <: Ring, S}(a::ZZ, b::PowerSeries{T, S})
   len = b.data.length
   d = Array(T, len)
   for i = 1:len
      d[i] = a*coeff(b, i - 1)
   end
   z = PowerSeries(PowerSeries{T, S}, d, b.prec)
   z.data.length = normalise(z, len)
   return z
end

*{T <: Ring, S}(a::PowerSeries{T, S}, b::Int) = b*a

*{T <: Ring, S}(a::PowerSeries{T, S}, b::ZZ) = b*a

###########################################################################################
#
#   Shifting
#
###########################################################################################

function shift_left{T <: Ring, S}(x::PowerSeries{T, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = x.data.length
   v = Array(T, xlen + len)
   for i = 1:len
      v[i] = zero(T)
   end
   for i = 1:xlen
      v[i + len] = coeff(x, i - 1)
   end
   return PowerSeries(PowerSeries{T, S}, v, x.prec + len)
end

function shift_right{T <: Ring, S}(x::PowerSeries{T, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = x.data.length
   if len >= xlen
      return PowerSeries(PowerSeries{T, S}, Array(T, 0), max(0, x.prec - len))
   end
   v = Array(T, xlen - len)
   for i = 1:xlen - len
      v[i] = coeff(x, i + len - 1)
   end
   return PowerSeries(PowerSeries{T, S}, v, x.prec - len)
end

###########################################################################################
#
#   Truncation
#
###########################################################################################

function truncate{T<: Ring, S}(x::PowerSeries{T, S}, prec::Precision)
   prec < 0 && throw(DomainError())
   if x.prec <= prec
      return x
   end
   d = Array(T, prec)
   for i = 1:min(prec, x.data.length)
      d[i] = coeff(x, i - 1)
   end
   for i = x.data.length + 1:prec
      d[i] = zero(T)
   end
   r = PowerSeries(PowerSeries{T, S}, d, prec)
   r.data.length = normalise(r, prec)
   return r
end

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{T <: Ring, S}(a::PowerSeries{T, S}, b::Int)
   b < 0 && throw(DomainError())
   # special case powers of x for constructing power series efficiently
   if a.prec == nothing && length(a) == 2 && coeff(a, 0) == 0 && coeff(a, 1) == 1
      d = Array(T, b + 1)
      d[b + 1] = coeff(a, 1)
      for i = 1:b
         d[i] = coeff(a, 0)
      end
      return PowerSeries(PowerSeries{T, S}, d, nothing)
   elseif length(a) == 0
      return PowerSeries(PowerSeries{T, S}, Array(T, 0), a.prec + (b - 1)*valuation(a))
   elseif length(a) == 1
      return PowerSeries(PowerSeries{T, S}, [a.data.coeffs[1]^b], a.prec)
   elseif b == 0
      return PowerSeries(PowerSeries{T, S}, [T(1)], nothing)
   else
      bit = ~((~uint(0)) >> 1)
      while (int(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit !=0
         z = z*z
         if (int(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end

###########################################################################################
#
#   Comparisons
#
###########################################################################################

=={T <: Ring, S}(x::PowerSeries{T, S}, y::Int) = x.prec == 0 || ((length(x) == 0 && y == 0)
                                        || (length(x) == 1 && coeff(x, 0) == y))

=={T <: Ring, S}(x::PowerSeries{T, S}, y::ZZ) = x.prec == 0 || ((length(x) == 0 && y == 0)
                                        || (length(x) == 1 && coeff(x, 0) == y))

function =={T<: Ring, S}(x::PowerSeries{T, S}, y::PowerSeries{T, S})
   prec = min(x.prec, y.prec)
   
   m1 = min(x.data.length, y.data.length)
   m2 = max(x.data.length, y.data.length)
   
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
         if coeffs(y, i - 1) != 0
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

=={T<: Ring, S}(x::Int, y::PowerSeries{T, S}) = y == x

=={T<: Ring, S}(x::ZZ, y::PowerSeries{T, S}) = y == x

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{T<: Ring, S}(x::PowerSeries{T, S}, y::PowerSeries{T, S})
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

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv{T<: Ring, S}(a::PowerSeries{T, S})
   a == 0 && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
   a1 = coeff(a, 0)
   if a.prec == nothing
      a.data.length != 1 && error("Unable to invert infinite precision power series")
      return PowerSeries(PowerSeries{T, S}, [divexact(T(1), a1)], nothing)
   end
   d = Array(T, a.prec)
   if a.prec != 0
      d[1] = divexact(T(1), a1)
   end
   a1 = -a1
   for n = 2:a.prec
      s = coeff(a, 1)*d[n - 1]
      for i = 2:min(n, a.data.length) - 1
         s += coeff(a, i)*d[n - i]
      end
      d[n] = divexact(s, a1)
   end
   ainv = PowerSeries(PowerSeries{T, S}, d, a.prec)
   ainv.data.length = normalise(ainv, a.prec)
   return ainv
end

###########################################################################################
#
#   Exponential
#
###########################################################################################

function exp{T <: Ring}(a::T)
   a != 0 && error("Exponential of nonzero element")
   return one(T)
end

function exp{T<: Ring, S}(a::PowerSeries{T, S})
   if a == 0
      return PowerSeries(PowerSeries{T, S}, [one(T)], a.prec)
   elseif a.prec == nothing
      error("Unable to compute exponential of infinite precision power series")
   end
   d = Array(T, a.prec)
   d[0+1] = exp(coeff(a, 0))
   len = length(a)
   for k = 1 : a.prec-1
      s = zero(T)
      for j = 1 : min(k+1, len) - 1
         s += j * coeff(a, j) * d[k-j+1]
      end
      d[k+1] = s / k
   end
   b = PowerSeries(PowerSeries{T, S}, d, a.prec)
   return b
end

###########################################################################################
#
#   PowerSeriesRing constructor
#
###########################################################################################

function PowerSeriesRing{T <: Ring}(::Type{T}, s::String)
   S = symbol(s)
   T1 = PowerSeries{T, S}
   T2 = T
   
   # Conversions and promotions

   Base.convert(::Type{T1}, x::T) = PowerSeries(T1, [x], nothing)
   Base.promote_rule(::Type{T1}, ::Type{T}) = T1

   P = T2.parameters
   while length(P) > 0
      T2 = P[1]
      if isa(T2, DataType) && T2 <: Ring
         Base.convert(::Type{T1}, x::T2) = PowerSeries(T1, [convert(T, x)], nothing)
         Base.promote_rule(::Type{T1}, ::Type{T2}) = T1
         P = T2.parameters
      else
         break
      end
   end

   Base.convert(::Type{T1}, x::Integer) = PowerSeries(T1, [convert(T, x)], nothing)
   Base.promote_rule{R <: Integer}(::Type{T1}, ::Type{R}) = T1

   # (Type, gen) 

   return (PowerSeries{T, S}, PowerSeries(PowerSeries{T, S}, [T(0), T(1)], nothing))
end
