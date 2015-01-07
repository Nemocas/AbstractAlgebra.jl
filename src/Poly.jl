###########################################################################################
#
#   Poly.jl : Generic polynomials over rings
#
###########################################################################################

export Poly, PolyRing, PolynomialRing, coeff, isgen, truncate, mullow, reverse, shift_left,
       shift_right, divexact, pseudorem, pseudodivrem, gcd, content, primpart, evaluate,
       compose, derivative, resultant, discriminant, bezout, zero, one, gen, length,
       iszero, normalise, isone

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

   Poly(R::Ring) = new(Array(T, 0), 0, PolyRing{Poly{T, S}, S}(R))
   
   Poly(R::Ring, a::Array{T, 1}) = new(a, length(a), PolyRing{Poly{T, S}, S}(R))

   function Poly(R::Ring, a::T) 
      par = PolyRing{Poly{T, S}, S}(R)
      return a == 0 ? new(Array(T, 0), 0, par) : new([a], 1, par)
   end

   function Poly(R::Ring, a::Integer)
      par = PolyRing{Poly{T, S}, S}(R)
      return a == 0 ? new(Array(T, 0), 0, par) : new([R(a)], 1, par)
   end
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

gen{T <: RingElem, S}(a::PolyRing{T, S}) = elem_type(a)(base(a), [zero(base(a)), one(base(a))])

iszero(a::PolyElem) = length(a) == 0

isone(a::PolyElem) = length(a) == 1 && isone(coeff(a, 0))

isgen(a::PolyElem) = length(a) == 2 && iszero(coeff(a, 0)) && isone(coeff(a, 1))

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

   # Conversions and promotions

   # conversion from base type
   eval(:(Base.convert(::Type{$T}, x::$T2) = $T($R, x)))
   eval(:(Base.promote_rule(::Type{$T}, ::Type{$T2}) = $T))

   # conversion from base type of base_rings, recursively
   R2 = R
   while base(R2) != None
      R2 = base(R2)
      T3 = elem_type(R2)
      eval(:(Base.convert(::Type{$T}, x::$T3) = $T($R, convert($T2, x))))
      eval(:(Base.promote_rule(::Type{$T}, ::Type{$T3}) = $T))
      eval(:(Base.call(a::$P, x::$T3) = $T($R, convert($T2, x))))
   end

   # conversion from Integer types
   if R != ZZ
      eval(:(Base.convert(::Type{$T}, x::BigInt) = $T(x)))
      eval(:(Base.promote_rule(::Type{$T}, ::Type{BigInt}) = $T))
   end

   eval(:(Base.convert(::Type{$T}, x::Integer) = $T(x)))
   eval(:(Base.promote_rule{Tint <: Integer}(::Type{$T}, ::Type{Tint}) = $T))

   # overload parent call for all of the above
   eval(:(Base.call(a::$P) = $T($R))) 
   eval(:(Base.call(a::$P, x::Int) = $T($R, x)))
   eval(:(Base.call(a::$P, x::Integer) = $T($R, BigInt(x))))
   eval(:(Base.call(a::$P, x::BigInt) = $T($R, x)))
   eval(:(Base.call(a::$P, x::$T2) = $T($R, x)))
   eval(:(Base.call(a::$P, x::$T) = x))

   return P(R), T(R, [C0, C1])
end