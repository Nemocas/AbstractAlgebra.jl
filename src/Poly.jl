###########################################################################################
#
#   Poly.jl : Generic polynomials over rings
#
###########################################################################################

export Poly, PolyRing, PolynomialRing, coeff, isgen, truncate, mullow, reverse, shift_left,
       shift_right, divexact, pseudorem, pseudodivrem, gcd, content, primpart, evaluate,
       compose, derivative, resultant, discriminant, bezout, zero, one, gen, length,
       iszero, normalise, isone, isunit, addeq!, mul!, fit!, setcoeff!

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

PolyID = ObjectIdDict()

type PolyRing{P <: PolyElem, S} <: Ring
   base_ring :: Ring

   function PolyRing(R::Ring)
      try
         PolyID[R, S]
      catch
         PolyID[R, S] = new(R)
      end
   end
end

type Poly{T <: RingElem, S} <: PolyElem
   coeffs::Array{T, 1}
   length::Int
   parent::PolyRing

   Poly(R::Ring) = new(Array(T, 0), 0)
   
   Poly(R::Ring, a::Array{T, 1}) = new(a, length(a))

   Poly(R::Ring, a::T) = a == 0 ? new(Array(T, 0), 0) : new([a], 1)

   Poly(R::Ring, a::Integer) = a == 0 ? new(Array(T, 0), 0) : new([R(a)], 1)
end

elem_type{P <: PolyElem, S}(::PolyRing{P, S}) = P

base(a::PolyRing) = a.base_ring

base(a::PolyElem) = base(parent(a))

function parent(a::PolyElem)
   return a.parent
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################    

function normalise(a::PolyElem, len::Int)
   while len > 0 && iszero(a.coeffs[len]) # cannot use coeff(a, len - 1)
      len -= 1
   end

   return len
end

length(x::PolyElem) = x.length

degree(x::PolyElem) = length(x) - 1

coeff(a::PolyElem, n::Int) = n >= length(a) ? base(a)(0) : a.coeffs[n + 1]

lead(a::PolyElem) = length(a) == 0 ? base(a)(0) : coeff(a, length(a) - 1)

zero(a::PolyRing) = a(0)

one(a::PolyRing) = a(1)

gen(a::PolyRing) = a([zero(base(a)), one(base(a))])

iszero(a::PolyElem) = length(a) == 0

isone(a::PolyElem) = length(a) == 1 && isone(coeff(a, 0))

isgen(a::PolyElem) = length(a) == 2 && iszero(coeff(a, 0)) && isone(coeff(a, 1))

isunit(a::PolyElem) = length(a) == 1 && isunit(coeff(a, 0))

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{T <: RingElem, S}(io::IO, x::Poly{T, S})
   len = length(x)

   if len == 0
      print(io, base(x)(0))
   else
      for i = 1:len - 1
         c = x.coeffs[len - i + 1]
         bracket = needs_parentheses(c)
         if !iszero(c)
            if i != 1 && !is_negative(c)
               print(io, "+")
            end
            if !isone(c) && (c != -1 || show_minus_one(typeof(c)))
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
            if len - i != 1
               print(io, "^")
               print(io, len - i)
            end
         end
      end
      c = x.coeffs[1]
      bracket = needs_parentheses(c)
      if !iszero(c)
         if len != 1 && !is_negative(c)
            print(io, "+")
         end
         if bracket
            print(io, "(")
         end
         show(io, c)
         if bracket
            print(io, ")")
         end
      end
   end
end

function show{P <: PolyElem, S}(io::IO, p::PolyRing{P, S})
   print(io, "Univariate Polynomial Ring in ")
   print(io, string(S))
   print(io, " over ")
   show(io, p.base_ring)
end

needs_parentheses(x::PolyElem) = length(x) > 1

is_negative(x::PolyElem) = length(x) <= 1 && is_negative(coeff(x, 0))

show_minus_one{T <: RingElem, S}(::Type{Poly{T, S}}) = show_minus_one(T)

###########################################################################################
#
#   Unary operations
#
###########################################################################################

function -{T <: RingElem, S}(a::Poly{T, S})
   len = length(a)
   d = Array(T, len)
   for i = 1:len
      d[i] = -a.coeffs[i]
   end
   z = parent(a)(d)
   z.length = len
   return z
end

###########################################################################################
#
#   Binary operations
#
###########################################################################################

function +{T <: RingElem, S}(a::Poly{T, S}, b::Poly{T, S})
   lena = length(a)
   lenb = length(b)
   lenz = max(lena, lenb)
   d = Array(T, lenz)
   i = 1

   while i <= min(lena, lenb)
      d[i] = a.coeffs[i] + b.coeffs[i]
      i += 1
   end

   while i <= lena
      d[i] = a.coeffs[i]
      i += 1
   end

   while i <= lenb
      d[i] = b.coeffs[i]
      i += 1
   end

   z = parent(a)(d)

   z.length = normalise(z, i - 1)

   return z
end

function -{T <: RingElem, S}(a::Poly{T, S}, b::Poly{T, S})
   lena = length(a)
   lenb = length(b)
   lenz = max(lena, lenb)
   d = Array(T, lenz)
   i = 1

   while i <= min(lena, lenb)
      d[i] = a.coeffs[i] - b.coeffs[i]
      i += 1
   end

   while i <= lena
      d[i] = a.coeffs[i]
      i += 1
   end

   while i <= lenb
      d[i] = -b.coeffs[i]
      i += 1
   end

   z = parent(a)(d)

   z.length = normalise(z, i - 1)

   return z
end

function *{T <: RingElem, S}(a::Poly{T, S}, b::Poly{T, S})
   lena = length(a)
   lenb = length(b)

   if lena == 0 || lenb == 0
      return parent(a)()
   end

   t = base(a)()

   lenz = lena + lenb - 1
   d = Array(T, lenz)
   
   for i = 1:lena
      d[i] = a.coeffs[i]*b.coeffs[1]
   end

   for i = 2:lenb
      d[lena + i - 1] = a.coeffs[lena]*b.coeffs[i]
   end
   
   for i = 1:lena - 1
      for j = 2:lenb
         mul!(t, a.coeffs[i], b.coeffs[j])
         addeq!(d[i + j - 1], t)
      end
   end
   
   z = parent(a)(d)
        
   z.length = normalise(z, lenz)

   return z
end

###########################################################################################
#
#   Ad hoc binary operators
#
###########################################################################################

function *{T <: RingElem, S}(a::Int, b::Poly{T, S})
   len = length(b)
   d = Array(T, len)
   for i = 1:len
      d[i] = a*coeff(b, i - 1)
   end
   z = parent(b)(d)
   z.length = normalise(z, len)
   return z
end

function *{T <: RingElem, S}(a::BigInt, b::Poly{T, S})
   len = length(b)
   d = Array(T, len)
   for i = 1:len
      d[i] = a*coeff(b, i - 1)
   end
   z = parent(b)(d)
   z.length = normalise(z, len)
   return z
end

*(a::Poly, b::Int) = b*a

*(a::Poly, b::BigInt) = b*a

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{T <: RingElem, S}(a::Poly{T, S}, b::Int)
   b < 0 && throw(DomainError())
   # special case powers of x for constructing polynomials efficiently
   if isgen(a)
      d = Array(T, b + 1)
      d[b + 1] = a.coeffs[2]
      for i = 1:b
         d[i] = a.coeffs[1]
      end
      z = parent(a)(d)
      z.length = b + 1
      return z
   elseif length(a) == 0
      return zero(parent(a))
   elseif length(a) == 1
      return parent(a)(a.coeffs[1]^b)
   elseif b == 0
      return one(parent(a))
   else
      bit = ~((~uint(0)) >> 1)
      while (uint(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit !=0
         z = z*z
         if (uint(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end

###########################################################################################
#
#   Unsafe functions
#
###########################################################################################

function fit!{T <: RingElem, S}(c::Poly{T, S}, n::Int)
   if length(c) < n
      t = c.coeffs
      c.coeffs = Array(T, n)
      for i = 1:length(c)
         c.coeffs[i] = t[i]
      end
      for i = length(c) + 1:n
         c.coeffs[i] = zero(base(c))
      end
   end
end

function setcoeff!{T <: RingElem, S}(c::Poly{T, S}, n::Int, a::T)
   if a != 0 || n + 1 <= length(c)
      fit!(c, n + 1)
      c.coeffs[n + 1] = a
      c.length = max(length(c), n + 1)
      # don't normalise
   end
end

function mul!(c::Poly, a::Poly, b::Poly)
   lena = length(a)
   lenb = length(b)

   if lena == 0 || lenb == 0
      c.length = 0
   else
      t = base(a)()

      lenc = lena + lenb - 1
      fit!(c, lenc)

      for i = 1:lena
         mul!(c.coeffs[i], a.coeffs[i], b.coeffs[1])
      end

      for i = 2:lenb
         mul!(c.coeffs[lena + i - 1], a.coeffs[lena], b.coeffs[i])
      end

      for i = 1:lena - 1
         for j = 2:lenb
            mul!(t, a.coeffs[i], b.coeffs[j])
            addeq!(c.coeffs[i + j - 1], t)
         end
      end
        
      c.length = normalise(c, lenc)
   end
end

function addeq!(c::Poly, a::Poly)
   lenc = length(c)
   lena = length(a)
   len = max(lenc, lena)
   fit!(c, len)
   for i = 1:lena
      addeq!(c.coeffs[i], a.coeffs[i])
   end
   c.length = normalise(c, len)
end

###########################################################################################
#
#   PolynomialRing constructor
#
###########################################################################################

function PolynomialRing(R::Ring, s::String)
   S = symbol(s)
   C0 = R(0)
   C1 = R(1)
   T2 = elem_type(R)
   T = Poly{T2, S}
   P = PolyRing{T, S}
   par = P(R)
   
   # Conversions and promotions

   # conversion from base type
   eval(:(Base.convert(::Type{$T}, x::$T2) = $par(x)))
   eval(:(Base.promote_rule(::Type{$T}, ::Type{$T2}) = $T))

   # conversion from base type of base_rings, recursively
   R2 = R
   while base(R2) != None
      R2 = base(R2)
      T3 = elem_type(R2)
      eval(:(Base.convert(::Type{$T}, x::$T3) = $par(convert($T2, x))))
      eval(:(Base.promote_rule(::Type{$T}, ::Type{$T3}) = $T))
      eval(:(Base.call(a::$P, x::$T3) = begin z = $par(convert($T2, x)); z.parent = $par; return z; end))
   end

   # conversion from Integer types
   if R != ZZ
      eval(:(Base.convert(::Type{$T}, x::BigInt) = $par(x)))
      eval(:(Base.promote_rule(::Type{$T}, ::Type{BigInt}) = $T))
   end

   eval(:(Base.convert(::Type{$T}, x::Integer) = $par(x)))
   eval(:(Base.promote_rule{Tint <: Integer}(::Type{$T}, ::Type{Tint}) = $T))

   # overload parent call for all of the above
   eval(:(Base.call(a::$P) = begin z = $T($R); z.parent = $par; return z; end)) 
   eval(:(Base.call(a::$P, x::Int) = begin z = $T($R, x); z.parent = $par; return z; end))
   eval(:(Base.call(a::$P, x::Integer) = begin z = $T($R, BigInt(x)); z.parent = $par; return z; end))
   eval(:(Base.call(a::$P, x::BigInt) = begin z = $T($R, x); z.parent = $par; return z; end))
   eval(:(Base.call(a::$P, x::$T2) = begin z = $T($R, x); z.parent = $par; return z; end))
   eval(:(Base.call(a::$P, x::$T) = begin z = x; z.parent = $par; return z; end))
   eval(:(Base.call(a::$P, x::Array{$T2, 1}) = begin z = $T($R, x); z.parent = $par; return z; end))

   return P(R), P(R)([C0, C1])
end