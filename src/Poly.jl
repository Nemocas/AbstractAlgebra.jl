###########################################################################################
#
#   Poly.jl : Generic polynomials over rings
#
###########################################################################################

export Poly, PolynomialRing, coeff, isgen, truncate, mullow, reverse, shift_left,
       shift_right, divexact, pseudorem, pseudodivrem, gcd, content, primpart, evaluate,
       compose, derivative, resultant, discriminant, bezout, zero, one, gen, length,
       iszero, normalise, isone, isunit, addeq!, mul!, fit!, setcoeff!, mulmod, powmod,
       invmod, lcm, divrem, mod, gcdinv, hash, canonical_unit

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

PolyID = ObjectIdDict()

type PolynomialRing{T <: RingElem, S} <: Ring
   base_ring :: Ring

   function PolynomialRing(R::Ring)
      return try
         PolyID[R, S]
      catch
         PolyID[R, S] = new(R)
      end
   end
end

type Poly{T <: RingElem, S} <: PolyElem
   coeffs::Array{T, 1}
   length::Int
   parent::PolynomialRing{T, S}

   Poly(R::Ring) = new(Array(T, 0), 0)
   
   Poly(R::Ring, a::Array{T, 1}) = new(a, length(a))

   Poly(R::Ring, a::T) = a == 0 ? new(Array(T, 0), 0) : new([a], 1)

   Poly(R::Ring, a::Integer) = a == 0 ? new(Array(T, 0), 0) : new([R(a)], 1)
end

elem_type{T <: RingElem, S}(::PolynomialRing{T, S}) = Poly{T, S}

base_ring(a::PolynomialRing) = a.base_ring

base_ring(a::PolyElem) = base_ring(parent(a))

parent(a::PolyElem) = a.parent

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################    

function hash(a::PolyElem)
   h = 0x53dd43cd511044d1
   for i in 0:length(a) - 1
      h $= hash(coeff(a, i))
      h = (h << 1) | (h >> (sizeof(Int)*8 - 1))
   end
   return h
end

function normalise(a::PolyElem, len::Int)
   while len > 0 && iszero(a.coeffs[len]) # cannot use coeff(a, len - 1)
      len -= 1
   end

   return len
end

length(x::PolyElem) = x.length

degree(x::PolyElem) = length(x) - 1

coeff(a::PolyElem, n::Int) = n >= length(a) ? base_ring(a)(0) : a.coeffs[n + 1]

lead(a::PolyElem) = length(a) == 0 ? base_ring(a)(0) : coeff(a, length(a) - 1)

zero(a::PolynomialRing) = a(0)

one(a::PolynomialRing) = a(1)

gen(a::PolynomialRing) = a([zero(base_ring(a)), one(base_ring(a))])

iszero(a::PolyElem) = length(a) == 0

isone(a::PolyElem) = length(a) == 1 && isone(coeff(a, 0))

isgen(a::PolyElem) = length(a) == 2 && iszero(coeff(a, 0)) && isone(coeff(a, 1))

isunit(a::PolyElem) = length(a) == 1 && isunit(coeff(a, 0))

###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit{T <: RingElem, S}(x::Poly{T, S}) = canonical_unit(lead(x))

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{T <: RingElem, S}(io::IO, x::Poly{T, S})
   len = length(x)

   if len == 0
      print(io, base_ring(x)(0))
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

function show{T <: RingElem, S}(io::IO, p::PolynomialRing{T, S})
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

   t = base_ring(a)()

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
#   Comparisons
#
###########################################################################################

function =={T <: RingElem, S}(x::Poly{T, S}, y::Poly{T, S})
   if length(x) != length(y)
      return false
   else
      for i = 1:length(x)
         if x.coeffs[i] != y.coeffs[i]
            return false
         end
      end
   end
   return true
end

function isequal{T <: RingElem, S}(x::Poly{T, S}, y::Poly{T, S})
   if length(x) != length(y)
      return false
   end
   for i = 1:length(x)
      if !isequal(coeff(x, i - 1), coeff(y, i - 1))
         return false
      end
   end
   return true
end

###########################################################################################
#
#   Ad hoc comparisons
#
###########################################################################################

==(x::Poly, y::Integer) = ((length(x) == 0 && y == 0)
                        || (length(x) == 1 && coeff(x, 0) == y))

==(x::Integer, y::Poly) = y == x

###########################################################################################
#
#   Truncation
#
###########################################################################################

function truncate{T <: RingElem, S}(a::Poly{T, S}, n::Int)
   n < 0 && throw(DomainError())
   
   lena = length(a)

   if lena <= n
      return a
   end

   lenz = min(lena, n)
   d = Array(T, lenz)

   for i = 1:lenz
      d[i] = a.coeffs[i]
   end

   z = parent(a)(d)

   z.length = normalise(z, lenz)

   return z
end

function mullow{T <: RingElem, S}(a::Poly{T, S}, b::Poly{T, S}, n::Int)
   lena = length(a)
   lenb = length(b)

   if lena == 0 || lenb == 0
      return zero(parent(a))
   end

   if n < 0
      n = 0
   end

   t = T()

   lenz = min(lena + lenb - 1, n)

   d = Array(T, lenz)

   for i = 1:min(lena, lenz)
      d[i] = a.coeffs[i]*b.coeffs[1]
   end

   if lenz > lena
      for j = 2:min(lenb, lenz - lena + 1)
          d[lena + j - 1] = a.coeffs[lena]*b.coeffs[j]
      end
   end

   z = parent(a)(d)

   for i = 1:lena - 1
      if lenz > i
         for j = 2:min(lenb, lenz - i + 1)
            mul!(t, a.coeffs[i], b.coeffs[j])
            addeq!(z.coeffs[i + j - 1], t)
         end
      end
   end
        
   z.length = normalise(z, lenz)

   return z
end

###########################################################################################
#
#   Reversal
#
###########################################################################################

function reverse{T <: RingElem, S}(x::Poly{T, S}, len::Int)
   len < 0 && throw(DomainError())
   v = Array(T, len)
   for i = 1:len
      v[i] = coeff(x, len - i)
   end
   r = parent(x)(v)
   r.length = normalise(r, len)
   return r
end

function reverse(x::Poly)
   reverse(x, length(x))
end

###########################################################################################
#
#   Shifting
#
###########################################################################################

function shift_left{T <: RingElem, S}(x::Poly{T, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   v = Array(T, xlen + len)
   for i = 1:len
      v[i] = zero(base_ring(x))
   end
   for i = 1:xlen
      v[i + len] = coeff(x, i - 1)
   end
   return parent(x)(v)
end

function shift_right{T <: RingElem, S}(x::Poly{T, S}, len::Int)
   len < 0 && throw(DomainError())
   xlen = length(x)
   if len >= xlen
      return zero(parent(x))
   end
   v = Array(T, xlen - len)
   for i = 1:xlen - len
      v[i] = coeff(x, i + len - 1)
   end
   return parent(x)(v)
end

###########################################################################################
#
#   Modular arithmetic
#
###########################################################################################

function mulmod{T <: Union(Residue, FieldElem), S}(a::Poly{T, S}, b::Poly{T, S}, d::Poly{T, S})
   return mod(a*b, d)
end

function powmod{T <: Union(Residue, FieldElem), S}(a::Poly{T, S}, b::Int, d::Poly{T, S})
   if length(a) == 0
      return zero(parent(a))
   elseif length(a) == 1
      return parent(a)(coeff(a, 0)^b)
   elseif b == 0
      return one(parent(a))
   else
      if b < 0
         a = invmod(a, d)
         b = -b
      end
      bit = ~((~uint(0)) >> 1)
      while (uint(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit !=0
         z = mulmod(z, z, d)
         if (uint(bit) & b) != 0
            z = mulmod(z, a, d)
         end
         bit >>= 1
      end
      return z
   end
end

function invmod{T <: Union(Residue, FieldElem), S}(a::Poly{T, S}, b::Poly{T, S})
   g, z = gcdinv(a, b)
   if g != 1
      error("Impossible inverse in invmod")
   end
   return z
end

###########################################################################################
#
#   Exact division
#
###########################################################################################

function divexact{T <: RingElem, S}(f::Poly{T, S}, g::Poly{T, S})
   g == 0 && throw(DivideError())
   if f == 0
      return zero(parent(f))
   end
   lenq = length(f) - length(g) + 1
   d = Array(T, lenq)
   for i = 1:lenq
      d[i] = zero(base_ring(f))
   end
   q = parent(f)(d)
   x = gen(parent(f))
   leng = length(g)
   while length(f) >= leng
      lenf = length(f)
      q1 = q.coeffs[lenf - leng + 1] = divexact(f.coeffs[lenf], g.coeffs[leng])
      f = f - q1*g*x^(lenf - leng)
   end
   q.length = lenq
   return q
end

###########################################################################################
#
#   Ad hoc exact division
#
###########################################################################################

function divexact{T <: RingElem, S}(a::Poly{T, S}, b::T)
   b == 0 && throw(DivideError())
   d = Array(T, length(a))
   for i = 1:length(a)
      d[i] = divexact(coeff(a, i - 1), b)
   end
   z = parent(a)(d)
   z.length = length(a)
   return z
end

function divexact{T <: RingElem, S}(a::Poly{T, S}, b::Integer)
   b == 0 && throw(DivideError())
   d = Array(T, length(a))
   for i = 1:length(a)
      d[i] = divexact(coeff(a, i - 1), b)
   end
   z = parent(a)(d)
   z.length = length(a)
   return z
end

###########################################################################################
#
#   Euclidean division
#
###########################################################################################

function mod{T <: Union(Residue, FieldElem), S}(f::Poly{T, S}, g::Poly{T, S})
   if length(g) == 0
      raise(DivideError())
   end
   if length(f) >= length(g)
      b = g.coeffs[length(g)]
      g = inv(b)*g
      x = gen(parent(f))
      while length(f) >= length(g)
         f -= coeff(f, length(f) - 1)*g*x^(length(f) - length(g))
      end
   end
   return f
end

function divrem{T <: Union(FieldElem, Residue), S}(f::Poly{T, S}, g::Poly{T, S})
   if length(g) == 0
      raise(DivideError())
   end
   if length(f) < length(g)
      return zero(parent(f)), f
   end
   binv = inv(lead(g))
   g = binv*g
   x = gen(parent(f))
   qlen = length(f) - length(g) + 1
   d = Array(T, qlen)
   for i = 1:qlen
      d[i] = zero(base_ring(f))
   end
   q = parent(f)(d)
   while length(f) >= length(g)
      q1 = coeff(f, length(f) - 1)
      setcoeff!(q, length(f) - length(g), q1*binv)
      f -= q1*g*x^(length(f) - length(g))
   end
   return q, f
end

###########################################################################################
#
#   Pseudodivision
#
###########################################################################################

function pseudorem{T <: RingElem, S}(f::Poly{T, S}, g::Poly{T, S})
   g == 0 && throw(DivideError())
   b = coeff(g, length(g) - 1)
   x = gen(parent(f))
   while length(f) >= length(g)
      f = f*b - coeff(f, length(f) - 1)*g*x^(length(f) - length(g))
   end
   return f
end

function pseudodivrem{T <: RingElem, S}(f::Poly{T, S}, g::Poly{T, S})
   g == 0 && throw(DivideError())
   if length(f) < length(g)
      return zero(parent(f)), f
   end
   lenq = length(f) - length(g) + 1
   v = Array(T, lenq)
   for i = 1:lenq
      v[i] = zero(base_ring(f))
   end
   q = parent(f)(v)
   b = coeff(g, length(g) - 1)
   x = gen(parent(f))
   while length(f) >= length(g)
      for i = length(f) - length(g) + 2:lenq
         setcoeff!(q, i - 1, coeff(q, i - 1) * b)
      end
      setcoeff!(q, length(f) - length(g), coeff(f, length(f) - 1))
      f = f*b - coeff(f, length(f) - 1)*g*x^(length(f) - length(g))
   end
   while lenq > 0 && coeff(q, lenq - 1) == 0
      lenq -= 1
   end
   q.length = lenq
   return q, f
end

###########################################################################################
#
#   Content, primitive part, GCD and LCM
#
###########################################################################################

function gcd{T <: RingElem, S}(a::Poly{T, S}, b::Poly{T, S})
   if length(b) > length(a)
      (a, b) = (b, a)
   end
   if b == 0
      return a
   end
   c = gcd(content(a), content(b))
   a = divexact(a, c)
   b = divexact(b, c)
   g = one(parent(a))
   h = one(parent(a))
   while true
      d = length(a) - length(b)
      r = pseudorem(a, b)
      if r == 0
         break
      end
      if length(r) == 1
         b = one(parent(a))
         break
      end
      (a, b) = (b, divexact(r, g*h^d))
      g = lead(a)
      if d > 1
         h = divexact(g^d, h^(d - 1))
      else
         h = h^(1 - d)*g^d
      end
   end
   return c*primpart(b)
end

function gcd{T <: Union(FieldElem, Residue), S}(a::Poly{T, S}, b::Poly{T, S})
   if length(a) > length(b)
      (a, b) = (b, a)
   end
   if b == 0
      return a
   end
   g = gcd(content(a), content(b))
   a = divexact(a, g)
   b = divexact(b, g)
   while a != 0
      (a, b) = (mod(b, a), a)
   end
   b = g*b
   return inv(lead(b))*b
end

function lcm{T <: RingElem, S}(a::Poly{T, S}, b::Poly{T, S})
   return a*divexact(b, gcd(a, b))
end

function content(a::Poly)
   z = coeff(a, 0)
   for i = 2:length(a)
      z = gcd(z, coeff(a, i - 1))
   end
   return z
end

function primpart(a::Poly)
   d = content(a)
   return divexact(a, d)
end

###########################################################################################
#
#   Evaluation/composition
#
###########################################################################################

function evaluate{T <: RingElem, S}(a::Poly{T, S}, b::T)
   i = length(a)
   if i == 0
       return zero(base_ring(a))
   end
   z = a.coeffs[i]
   while i > 1
      i -= 1
      z = z*b + a.coeffs[i]
   end
   return z
end

function evaluate{T <: RingElem, S}(a::Poly{T, S}, b::Integer)
   return evaluate(a, base_ring(a)(b))
end

function compose(a::Poly, b::Poly)
   i = length(a)
   if i == 0
       return zero(parent(a))
   end
   z = a.coeffs[i]
   while i > 1
      i -= 1
      z = z*b + a.coeffs[i]
   end
   return z
end

###########################################################################################
#
#   Derivative
#
###########################################################################################

function derivative{T <: RingElem, S}(a::Poly{T, S})
   if a == 0
      return zero(parent(a))
   end
   len = length(a)
   d = Array(T, len - 1)
   for i = 1:len - 1
      d[i] = i*a.coeffs[i + 1]
   end
   z = parent(a)(d)
   z.length = normalise(z, len - 1)
   return z
end

###########################################################################################
#
#   Integral
#
###########################################################################################

function integral{T <: Union(FieldElem, Residue), S}(x::Poly{T, S})
   len = length(x)
   v = Array(T, len + 1)
   v[1] = zero(base_ring(x))
   for i = 1:len
      v[i + 1] = divexact(coeff(x, i - 1), T(i))
   end
   p = parent(x)(v)
   len = len + 1
   while len > 0 && coeff(p, len - 1) == 0 # cannot use normalise here
      len -= 1
   end
   return p
end

###########################################################################################
#
#   Resultant
#
###########################################################################################

function resultant{T <: RingElem, S}(a::Poly{T, S}, b::Poly{T, S})
   if length(a) == 0 || length(b) == 0
      return zero(base_ring(a))
   end
   sgn = 1
   if length(a) < length(b)
      a, b = b, a
      if iseven(length(a)) && iseven(length(b))
         sgn = -sgn
      end
   end
   lena = length(a)
   lenb = length(b)
   if lenb == 1
      return b.coeffs[1]^(lena - 1)
   end
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   g = one(base_ring(a))
   h = one(base_ring(a))
   while lenb > 1
      d = lena - lenb
      if iseven(lena) && iseven(lenb)
         sgn = -sgn
      end
      B, A = pseudorem(A, B), B
      lena = lenb
      lenb = length(B)
      if lenb == 0
         return zero(base_ring(a)) 
      end
      s = h^d
      B = divexact(B, g*s)
      g = lead(A)
      h = divexact(h*g^d, s)
   end
   s = divexact(h*lead(B)^(lena - 1), h^(lena - 1))
   res = c1^(lenb - 1)*c2^(lena - 1)*s*sgn
end

function resultant{T <: Union(FieldElem, Residue), S}(a::Poly{T, S}, b::Poly{T, S})
   if length(a) == 0 || length(b) == 0
      return zero(base_ring(a))
   end
   sgn = 1
   if length(a) < length(b)
      a, b = b, a
      if iseven(length(a)) && iseven(length(b))
         sgn = -sgn
      end
   end
   lena = length(a)
   lenb = length(b)
   if lenb == 1
      return b.coeffs[1]^(lena - 1)
   end
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   s = 1
   while lenb > 1
      if iseven(lena) && iseven(lenb)
         sgn = -sgn
      end
      B, A = mod(A, B), B
      s *= lead(A)^(lena - length(B))
      lena = lenb
      lenb = length(B)
      if lenb == 0
         return zero(base_ring(a)) 
      end
   end
   s *= lead(B)^(lena - 1)
   res = c1^(lenb - 1)*c2^(lena - 1)*s*sgn
end

###########################################################################################
#
#   Discriminant
#
###########################################################################################

function discriminant(a::Poly)
   d = deriv(a)
   z = resultant(a, d)
   if length(a) - length(d) == 1
      z = divexact(z, lead(a))
   else
      z = z*lead(a)^(length(a) - length(d) - 2)
   end
   mod4 = (length(a) + 3)%4 # degree mod 4
   return mod4 == 2 || mod4 == 3 ? -z : z
end

###########################################################################################
#
#   Bezout
#
###########################################################################################

function bezout{T <: RingElem, S}(a::Poly{T, S}, b::Poly{T, S})
   if length(a) == 0 || length(b) == 0
      return zero(base_ring(a)), zero(parent(a)), zero(parent(a))
   end
   sgn = 1
   swap = false
   if length(a) < length(b)
      a, b = b, a
      swap = true
      if iseven(length(a)) && iseven(length(b))
         sgn = -sgn
      end
   end
   lena = length(a)
   lenb = length(b)
   if lenb == 1
      s1 = zero(base_ring(a))
      t1 = one(base_ring(a))
      r1 = b.coeffs[1]^(lena - 1)
      if swap
         s1, t1 = t1, s1
      end
      if sgn
         s1, t1, r1 = -s1, -t1, -r1
      end
      return r1, parent(a)(s1), parent(a)(t1)
   end
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   g = one(base_ring(a))
   h = one(base_ring(a))
   u1, u2 = one(parent(a)), zero(parent(a))
   v1, v2 = zero(parent(a)), one(parent(a))
   while lenb > 1
      d = lena - lenb
      if iseven(lena) && iseven(lenb)
         sgn = -sgn
      end
      (Q, B), A = pseudodivrem(A, B), B
      lena = lenb
      lenb = length(B)
      if lenb == 0
         return zero(base_ring(a)), zero(parent(a)), zero(parent(a))
      end
      s = h^d
      B = divexact(B, g*s)
      t = lead(A)^(d + 1)
      u2, u1 = divexact(u1*t - Q*u2, g*s), u2
      v2, v1 = divexact(v1*t - Q*v2, g*s), v2 
      g = lead(A)
      h = divexact(h*g^d, s)
   end
   s = divexact(h*lead(B)^(lena - 1), h^(lena - 1))
   res = c1^(lenb - 1)*c2^(lena - 1)*s*sgn
   if swap
      u2, v2 = v2, u2
   end
   return res, u2, v2
end

function bezout{T <: Union(FieldElem, Residue), S}(a::Poly{T, S}, b::Poly{T, S})
   if length(a) == 0
      return b, zero(parent(a)), one(parent(a))
   end
   if length(b) == 0
      return a, one(parent(a)), zero(parent(a))
   end
   swap = false
   if length(a) < length(b)
      a, b = b, a
      swap = true
   end
   lena = length(a)
   lenb = length(b)
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   u1, u2 = inv(c1), zero(parent(a))
   v1, v2 = zero(parent(a)), inv(c2)
   while lenb > 0
      d = lena - lenb
      (Q, B), A = divrem(A, B), B
      lena = lenb
      lenb = length(B)
      u2, u1 = u1 - Q*u2, u2
      v2, v1 = v1 - Q*v2, v2 
   end
   if swap
      u1, v1 = v1, u1
   end
   d = gcd(c1, c2)
   A, u1, v1 = d*A, d*u1, d*v1
   d = inv(lead(A))
   return d*A, d*u1, d*v1
end

function gcdinv{T <: Union(FieldElem, Residue), S}(a::Poly{T, S}, b::Poly{T, S})
   if length(a) == 0
      if length(b) == 0
         return zero(base_ring(a)), zero(parent(a))
      else
         d = inv(lead(b))
         return b*d, zero(parent(a))
      end
   end
   if length(b) == 0
      d = inv(lead(b))
      return a*d, d
   end
   if length(a) < length(b)
      a, b = b, a
      u1, u2 = zero(parent(a)), one(parent(a))
   else
      u1, u2 = one(parent(a)), zero(parent(a))
   end
   lena = length(a)
   lenb = length(b)
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   u1 *= inv(c1)
   u2 *= inv(c2)
   while lenb > 0
      d = lena - lenb
      (Q, B), A = divrem(A, B), B
      lena = lenb
      lenb = length(B)
      u2, u1 = u1 - Q*u2, u2
   end
   d = gcd(c1, c2)
   A, u1 = d*A, d*u1
   d = inv(lead(A))
   return d*A, d*u1
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
         c.coeffs[i] = zero(base_ring(c))
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

function mul!{T <: RingElem, S}(c::Poly{T, S}, a::Poly{T, S}, b::Poly{T, S})
   lena = length(a)
   lenb = length(b)

   if lena == 0 || lenb == 0
      c.length = 0
   else
      t = base_ring(a)()

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

function addeq!{T <: RingElem, S}(c::Poly{T, S}, a::Poly{T, S})
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
#   Promotion rules
#
###########################################################################################

Base.promote_rule{T <: RingElem, S, V <: Integer}(::Type{Poly{T, S}}, ::Type{V}) = Poly{T, S}

Base.promote_rule{T <: RingElem, S}(::Type{Poly{T, S}}, ::Type{T}) = Poly{T, S}

###########################################################################################
#
#   Parent object call overload
#
###########################################################################################

function Base.call{T <: RingElem, S}(a::PolynomialRing{T, S}, b::RingElem)
   return a(base_ring(a)(b))
end

function Base.call{T <: RingElem, S}(a::PolynomialRing{T, S})
   z = Poly{T, S}(base_ring(a))
   z.parent = a
   return z
end

function Base.call{T <: RingElem, S}(a::PolynomialRing{T, S}, b::Int)
   z = Poly{T, S}(base_ring(a), b)
   z.parent = a
   return z
end

function Base.call{T <: RingElem, S}(a::PolynomialRing{T, S}, b::Integer)
   z = Poly{T, S}(base_ring(a), BigInt(b))
   z.parent = a
   return z
end

function Base.call{T <: RingElem, S}(a::PolynomialRing{T, S}, b::BigInt)
   z = Poly{T, S}(base_ring(a), b)
   z.parent = a
   return z
end

function Base.call{T <: RingElem, S}(a::PolynomialRing{T, S}, b::T)
   parent(b) != base_ring(a) && error("Unable to coerce element to polynomial")
   z = Poly{T, S}(base_ring(a), b)
   z.parent = a
   return z
end

function Base.call{T <: RingElem, S}(a::PolynomialRing{T, S}, b::Poly{T, S})
   parent(b) != a && error("Unable to coerce polynomial")
   z = Poly{T, S}(base_ring(a), b)
   z.parent = a
   return z
end

function Base.call{T <: RingElem, S}(a::PolynomialRing{T, S}, b::Array{T, 1})
   if length(b) > 0
      parent(b[1]) != base_ring(a) && error("Unable to coerce element to polynomial")
   end
   z = Poly{T, S}(base_ring(a), b)
   z.parent = a
   return z
end

###########################################################################################
#
#   PolynomialRing constructor
#
###########################################################################################

function PolynomialRing(R::Ring, s::String)
   S = symbol(s)
   T = elem_type(R)
   parent_obj = PolynomialRing{T, S}(R)

   base = base_ring(R)
   R2 = R
   parent_type = Poly{T, S}
   while base_ring(R2) != None
      R2 = base_ring(R2)
      T2 = elem_type(R2)
      eval(:(Base.promote_rule(::Type{$parent_type}, ::Type{$T2}) = $parent_type))
   end

   return parent_obj, parent_obj([R(0), R(1)])
end