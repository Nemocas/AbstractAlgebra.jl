export FractionField

base_ring(a::FracField{T}) where T <: RingElem = a.base_ring::parent_type(T)

base_ring(a::FracElem) = base_ring(parent(a))

parent(a::FracElem) = a.parent

function check_parent(a::FracElem, b::FracElem, throw::Bool = true)
   fl = parent(a) != parent(b)
   fl && throw && error("Incompatible rings in fraction field operation")
   return !fl
end

function //(x::T, y::T) where {T <: RingElem}
   R = parent(x)
   iszero(y) && throw(DivideError())
   g = gcd(x, y)
   z = Generic.Frac{T}(divexact(x, g), divexact(y, g))
   try
      z.parent = Generic.FracDict[R]
   catch
      z.parent = Generic.FractionField(R)
   end
   return z
end

//(x::T, y::FracElem{T}) where {T <: RingElem} = parent(y)(x)//y

//(x::FracElem{T}, y::T) where {T <: RingElem} = x//parent(x)(y)

function Base.hash(a::FracElem, h::UInt)
   b = 0x8a30b0d963237dd5%UInt
   # We canonicalise before hashing
   return xor(b, hash(numerator(a, true), h), hash(denominator(a, true), h), h)
end

# Fall back method for all other fraction types in system
function Base.numerator(a::FracElem, canonicalise::Bool=true)
   return Base.numerator(a) # all other types ignore canonicalise
end

# Fall back method for all other fraction types in system
function Base.denominator(a::FracElem, canonicalise::Bool=true)
   return Base.denominator(a) # all other types ignore canonicalise
end

zero(R::FracField) = R(0)

one(R::FracField) = R(1)

iszero(a::FracElem) = iszero(numerator(a, false))

isone(a::FracElem) = numerator(a, false) == denominator(a, false)

canonical_unit(a::FracElem) = a

function -(a::FracElem)
   return parent(a)(-numerator(a, false), deepcopy(denominator(a, false)))
end

function +(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   d1 = denominator(a, false)
   d2 = denominator(b, false)
   n1 = numerator(a, false)
   n2 = numerator(b, false)
   if d1 == d2
      rnum = n1 + n2
      if isone(d1)
         rden = deepcopy(d1)
      else
         gd = gcd(rnum, d1)
         if isone(gd)
            rden = deepcopy(d1)
         else
            rnum = divexact(rnum, gd)
            rden = divexact(d1, gd)
         end
      end
   elseif isone(d1)
      rnum = n1*d2 + n2
      rden = deepcopy(d2)
   elseif isone(d2)
      rnum = n1 + n2*d1
      rden = deepcopy(d1)
   else
      gd = gcd(d1, d2)
      if isone(gd)
         rnum = n1*d2 + n2*d1
         rden = d1*d2
      else
         q1 = divexact(d1, gd)
         q2 = divexact(d2, gd)
         rnum = q1*n2 + q2*n1
         t = gcd(rnum, gd)
         if isone(t)
            rden = q2*d1
         else
            rnum = divexact(rnum, t)
            gd = divexact(d1, t)
            rden = gd*q2
         end
      end
   end
   return parent(a)(rnum, rden)
end

function -(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   d1 = denominator(a, false)
   d2 = denominator(b, false)
   n1 = numerator(a, false)
   n2 = numerator(b, false)
   if d1 == d2
      rnum = n1 - n2
      if isone(d1)
         rden = deepcopy(d1)
      else
         gd = gcd(rnum, d1)
         if isone(gd)
            rden = deepcopy(d1)
         else
            rnum = divexact(rnum, gd)
            rden = divexact(d1, gd)
         end
      end
   elseif isone(d1)
      rnum = n1*d2 - n2
      rden = deepcopy(d2)
   elseif isone(d2)
      rnum = n1 - n2*d1
      rden = deepcopy(d1)
   else
      gd = gcd(d1, d2)
      if isone(gd)
         rnum = n1*d2 - n2*d1
         rden = d1*d2
      else
         q1 = divexact(d1, gd)
         q2 = divexact(d2, gd)
         rnum = q2*n1 - q1*n2
         t = gcd(rnum, gd)
         if isone(t)
            rden = q2*d1
         else
            rnum = divexact(rnum, t)
            gd = divexact(d1, t)
            rden = gd*q2
         end
      end
   end
   return parent(a)(rnum, rden)
end

function *(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   n1 = numerator(a, false)
   d2 = denominator(b, false)
   n2 = numerator(b, false)
   d1 = denominator(a, false)
   if d1 == d2
      n = n1*n2
      d = d1*d2
   elseif isone(d1)
      gd = gcd(n1, d2)
      if isone(gd)
         n = n1*n2
         d = deepcopy(d2)
      else
         n = divexact(n1, gd)*n2
         d = divexact(d2, gd)
      end
   elseif isone(d2)
      gd = gcd(n2, d1)
      if isone(gd)
         n = n2*n1
         d = deepcopy(d1)
      else
         n = divexact(n2, gd)*n1
         d = divexact(d1, gd)
      end
   else
      g1 = gcd(n1, d2)
      g2 = gcd(n2, d1)
      if !isone(g1)
         n1 = divexact(n1, g1)
         d2 = divexact(d2, g1)
      end
      if !isone(g2)
         n2 = divexact(n2, g2)
         d1 = divexact(d1, g2)
      end
      n = n1*n2
      d = d1*d2
   end
   return parent(a)(n, d)
end

function ==(x::FracElem{T}, y::FracElem{T}) where {T <: RingElem}
   b  = check_parent(x, y, false)
   !b && return false

   return (denominator(x, false) == denominator(y, false) &&
           numerator(x, false) == numerator(y, false)) ||
          (denominator(x, true) == denominator(y, true) &&
           numerator(x, true) == numerator(y, true)) ||
          (numerator(x, false)*denominator(y, false) ==
           denominator(x, false)*numerator(y, false))
end

function FractionField(R::Ring; cached=true)
   return Generic.FractionField(R; cached=cached)
end
