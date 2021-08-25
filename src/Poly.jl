export PolynomialRing, PolyRing

base_ring(R::PolyRing{T}) where T <: RingElement = R.base_ring::parent_type(T)

base_ring(a::PolynomialElem) = base_ring(parent(a))

parent(a::PolynomialElem) = a.parent

function check_parent(a::PolynomialElem, b::PolynomialElem, throw::Bool = true)
   c = parent(a) != parent(b)
   c && throw && error("Incompatible polynomial rings in polynomial operation")
   return !c
end

function Base.hash(a::PolyElem, h::UInt)
   b = 0x53dd43cd511044d1%UInt
   for i in 0:length(a) - 1
      b = xor(b, xor(hash(coeff(a, i), h), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

length(a::PolynomialElem) = a.length

zero(R::PolyRing) = R(0)

one(R::PolyRing) = R(1)

gen(R::PolyRing) = R([zero(base_ring(R)), one(base_ring(R))])

gens(R::PolyRing) = [gen(R)]

iszero(a::PolynomialElem) = length(a) == 0

isone(a::PolynomialElem) = length(a) == 1 && isone(coeff(a, 0))

canonical_unit(x::PolynomialElem) = canonical_unit(leading_coefficient(x))

function leading_coefficient(a::PolynomialElem)
   return length(a) == 0 ? zero(base_ring(a)) : coeff(a, length(a) - 1)
end

function trailing_coefficient(a::PolynomialElem)
   if iszero(a)
      return zero(base_ring(a))
   else
      for i = 1:length(a)
         c = coeff(a, i - 1)
         if !iszero(c)
            return c
         end
      end
      return coeff(a, length(a) - 1)
   end
end

function -(a::PolynomialElem)
   len = length(a)
   z = parent(a)()
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, -coeff(a, i - 1))
   end
   z = set_length!(z, len)
   return z
end

function +(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   lenz = max(lena, lenb)
   z = parent(a)()
   fit!(z, lenz)
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

function -(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   lenz = max(lena, lenb)
   z = parent(a)()
   fit!(z, lenz)
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

function mul_classical(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   lena = length(a)
   lenb = length(b)
   if lena == 0 || lenb == 0
      return parent(a)()
   end
   R = base_ring(a)
   t = R()
   lenz = lena + lenb - 1
   d = Array{T}(undef, lenz)
   for i = 1:lena
      d[i] = mul_red!(R(), coeff(a, i - 1), coeff(b, 0), false)
   end
   for i = 2:lenb
      d[lena + i - 1] = mul_red!(R(), coeff(a, lena - 1), coeff(b, i - 1), false)
   end
   for i = 1:lena - 1
      for j = 2:lenb
         t = mul_red!(t, coeff(a, i - 1), coeff(b, j - 1), false)
         d[i + j - 1] = addeq!(d[i + j - 1], t)
      end
   end
   for i = 1:lenz
      d[i] = reduce!(d[i])
   end
   z = parent(a)(d)
   z = set_length!(z, normalise(z, lenz))
   return z
end

function *(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   check_parent(a, b)
      return mul_classical(a, b)
end

function *(a::T, b::PolyElem{T}) where {T <: RingElem}
   len = length(b)
   z = parent(b)()
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, a*coeff(b, i - 1))
   end
   z = set_length!(z, normalise(z, len))
   return z
end

function *(a::Union{Integer, Rational, AbstractFloat}, b::PolynomialElem)
   len = length(b)
   z = parent(b)()
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, a*coeff(b, i - 1))
   end
   z = set_length!(z, normalise(z, len))
   return z
end

*(a::PolyElem{T}, b::T) where {T <: RingElem} = b*a

*(a::PolynomialElem, b::Union{Integer, Rational, AbstractFloat}) = b*a

function ^(a::PolyElem{T}, b::Int) where T <: RingElement
   b < 0 && throw(DomainError(b, "exponent must be >= 0"))
   # special case powers of x for constructing polynomials efficiently
   R = parent(a)
   if isgen(a)
      z = R()
      fit!(z, b + 1)
      z = setcoeff!(z, b, deepcopy(coeff(a, 1)))
      for i = 1:b
         z = setcoeff!(z, i - 1, deepcopy(coeff(a, 0)))
      end
      z = set_length!(z, b + 1)
      return z
   elseif b == 0
      return one(R)
   elseif length(a) == 0
      return zero(R)
   elseif length(a) == 1
      return R(coeff(a, 0)^b)
   elseif b == 1
      return deepcopy(a)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit != 0
         z = z*z
         if (UInt(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end

function ==(x::PolyElem{T}, y::PolyElem{T}) where T <: RingElement
   b = check_parent(x, y, false)
   !b && return false
   if length(x) != length(y)
      return false
   else
      for i = 1:length(x)
         if coeff(x, i - 1) != coeff(y, i - 1)
            return false
         end
      end
   end
   return true
end

==(x::PolyElem{T}, y::T) where T <: RingElem = ((length(x) == 0 && iszero(y))
                        || (length(x) == 1 && coeff(x, 0) == y))

==(x::PolynomialElem, y::Union{Integer, Rational, AbstractFloat}) = ((length(x) == 0 && iszero(base_ring(x)(y)))
                        || (length(x) == 1 && coeff(x, 0) == y))

==(x::T, y::PolyElem{T}) where T <: RingElem = y == x

==(x::Union{Integer, Rational, AbstractFloat}, y::PolyElem) = y == x

function shift_left(f::PolynomialElem, n::Int)
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   if n == 0
      return f
   end
   flen = length(f)
   r = parent(f)()
   fit!(r, flen + n)
   for i = 1:n
      r = setcoeff!(r, i - 1, zero(base_ring(f)))
   end
   for i = 1:flen
      r = setcoeff!(r, i + n - 1, coeff(f, i - 1))
   end
   return r
end

function shift_right(f::PolynomialElem, n::Int)
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   flen = length(f)
   if n >= flen
      return zero(parent(f))
   end
   if n == 0
      return f
   end
   r = parent(f)()
   fit!(r, flen - n)
   for i = 1:flen - n
      r = setcoeff!(r, i - 1, coeff(f, i + n - 1))
   end
   return r
end

function divexact(f::PolyElem{T}, g::PolyElem{T}; check::Bool=true) where T <: RingElement
   check_parent(f, g)
   iszero(g) && throw(DivideError())
   if iszero(f)
      return zero(parent(f))
   end
   lenq = length(f) - length(g) + 1
   d = Array{T}(undef, lenq)
   for i = 1:lenq
      d[i] = zero(base_ring(f))
   end
   x = gen(parent(f))
   leng = length(g)
   while length(f) >= leng
      lenf = length(f)
      q1 = d[lenf - leng + 1] = divexact(coeff(f, lenf - 1), coeff(g, leng - 1); check=check)
      f = f - shift_left(q1*g, lenf - leng)
      if length(f) == lenf # inexact case
         f = set_length!(f, normalise(f, lenf - 1))
      end
   end
   check && length(f) != 0 && throw(ArgumentError("not an exact division"))
   q = parent(f)(d)
   q = set_length!(q, lenq)
   return q
end

function divexact(a::PolyElem{T}, b::T; check::Bool=true) where {T <: RingElem}
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact(coeff(a, i - 1), b; check=check))
   end
   z = set_length!(z, length(a))
   return z
end

function divexact(a::PolyElem, b::Union{Integer, Rational, AbstractFloat}; check::Bool=true)
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact(coeff(a, i - 1), b; check=check))
   end
   z = set_length!(z, length(a))
   return z
end

function mod(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
   check_parent(f, g)
   if length(g) == 0
      throw(DivideError())
   end
   if length(f) >= length(g)
      f = deepcopy(f)
      b = leading_coefficient(g)
      g = inv(b)*g
      c = base_ring(f)()
      while length(f) >= length(g)
         l = -leading_coefficient(f)
         for i = 1:length(g) - 1
            c = mul!(c, coeff(g, i - 1), l)
            u = coeff(f, i + length(f) - length(g) - 1)
            u = addeq!(u, c)
            f = setcoeff!(f, i + length(f) - length(g) - 1, u)
         end
         f = set_length!(f, normalise(f, length(f) - 1))
      end
   end
   return f
end

function rem(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
  return mod(f, g)
end

function Base.divrem(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
   check_parent(f, g)
   if length(g) == 0
      throw(DivideError())
   end
   if length(f) < length(g)
      return zero(parent(f)), f
   end
   f = deepcopy(f)
   binv = inv(leading_coefficient(g))
   g = divexact(g, leading_coefficient(g))
   qlen = length(f) - length(g) + 1
   q = parent(f)()
   fit!(q, qlen)
   c = base_ring(f)()
   while length(f) >= length(g)
      q1 = leading_coefficient(f)
      l = -q1
      q = setcoeff!(q, length(f) - length(g), q1*binv)
      for i = 1:length(g) - 1
         c = mul!(c, coeff(g, i - 1), l)
         u = coeff(f, i + length(f) - length(g) - 1)
         u = addeq!(u, c)
         f = setcoeff!(f, i + length(f) - length(g) - 1, u)
      end
      f = set_length!(f, normalise(f, length(f) - 1))
   end
   return q, f
end

function Base.div(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
   q, r = divrem(f, g)
   return q
end

function Base.div(f::PolyElem{T}, g::T) where T <: Union{FieldElem, AbstractFloat, Rational}
   return div(f, parent(f)(g))
end

function pseudorem(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
  check_parent(f, g)
  iszero(g) && throw(DivideError())
  if length(f) < length(g)
     return f
  end
  k = length(f) - length(g) + 1
  b = coeff(g, length(g) - 1)
  x = gen(parent(f))
  while length(f) >= length(g)
     f = f*b - shift_left(coeff(f, length(f) - 1)*g, length(f) - length(g))
     k -= 1
  end
  return f*b^k
end

function pseudodivrem(f::PolyElem{T}, g::PolyElem{T}) where T <: RingElement
  check_parent(f, g)
  iszero(g) && throw(DivideError())
  if length(f) < length(g)
     return zero(parent(f)), f
  end
  lenq = length(f) - length(g) + 1
  k = lenq
  q = parent(f)()
  fit!(q, lenq)
  b = coeff(g, length(g) - 1)
  x = gen(parent(f))
  while length(f) >= length(g)
     for i = length(f) - length(g) + 2:lenq
        q = setcoeff!(q, i - 1, coeff(q, i - 1) * b)
     end
     q = setcoeff!(q, length(f) - length(g), coeff(f, length(f) - 1))
     f = f*b - shift_left(coeff(f, length(f) - 1)*g, length(f) - length(g))
     k -= 1
  end
  while lenq > 0 && iszero(coeff(q, lenq - 1))
     lenq -= 1
  end
  q = set_length!(q, lenq)
  s = b^k
  return q*s, f*s
end

function isterm(a::PolynomialElem)
   if iszero(a)
      return false
   end
   for i = 1:length(a) - 1
      if !iszero(coeff(a, i - 1))
         return false
      end
   end
   return true
end

isterm_recursive(a::T) where T <: RingElement = true

function isterm_recursive(a::PolynomialElem)
   if !isterm_recursive(leading_coefficient(a))
      return false
   end
   for i = 1:length(a) - 1
      if !iszero(coeff(a, i - 1))
         return false
      end
   end
   return true
end

function ismonomial(a::PolynomialElem)
   if !isone(leading_coefficient(a))
      return false
   end
   for i = 1:length(a) - 1
      if !iszero(coeff(a, i - 1))
         return false
      end
   end
   return true
end

ismonomial_recursive(a::T) where T <: RingElement = isone(a)

function ismonomial_recursive(a::PolynomialElem)
   if !ismonomial_recursive(leading_coefficient(a))
      return false
   end
   for i = 1:length(a) - 1
      if !iszero(coeff(a, i - 1))
         return false
      end
   end
   return true
end

function term_gcd(a::T, b::T) where T <: RingElement
   return gcd(a, b)
end

function term_content(a::T) where T <: RingElement
   return a
end

function term_gcd(a::PolyElem{T}, b::PolyElem{T}) where T <: RingElement
   d = min(degree(a), degree(b))
   x = gen(parent(a))
   return term_gcd(coeff(a, degree(a)), coeff(b, degree(b)))*x^d
end

function term_content(a::PolyElem{T}) where T <: RingElement
   for i = 1:length(a)
      c = coeff(a, i - 1)
      if !iszero(c)
         g = term_content(c)
         for j = i + 1:length(a)
            c = coeff(a, j - 1)
            if !iszero(c)
               g = term_gcd(g, term_content(c))
            end
         end
         x = gen(parent(a))
         return g*x^(i - 1)
      end
   end
   return parent(a)()
end

function gcd(a::PolyElem{T}, b::PolyElem{T}, ignore_content::Bool = false) where T <: RingElement
   check_parent(a, b)
   if length(b) > length(a)
      (a, b) = (b, a)
   end
   if iszero(b)
      if iszero(a)
         return a
      else
         return divexact(a, canonical_unit(leading_coefficient(a)))
      end
   end
   if isone(b)
      return b
   end
   if !ignore_content
      c1 = content(a)
      c2 = content(b)
      a = divexact(a, c1)
      b = divexact(b, c2)
      c = gcd(c1, c2)
   end
   lead_monomial = isterm_recursive(leading_coefficient(a)) ||
                   isterm_recursive(leading_coefficient(b))
   trail_monomial = isterm_recursive(trailing_coefficient(a)) ||
                    isterm_recursive(trailing_coefficient(b))
   lead_a = leading_coefficient(a)
   lead_b = leading_coefficient(b)
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
      g = leading_coefficient(a)
      if d > 1
         h = divexact(g^d, h^(d - 1))
      else
         h = h^(1 - d)*g^d
      end
   end
   if !ignore_content
      if !isterm_recursive(leading_coefficient(b)) &&
         !isterm_recursive(trailing_coefficient(b))
         if lead_monomial # lead term monomial, so content contains rest
            d = divexact(leading_coefficient(b),
	                 term_content(leading_coefficient(b)))
            b = divexact(b, d)
         elseif trail_monomial # trail term is monomial, so ditto
            d = divexact(trailing_coefficient(b),
			 term_content(trailing_coefficient(b)))
            b = divexact(b, d)
         else
            glead = gcd(lead_a, lead_b)
            if isterm_recursive(glead)
               d = divexact(leading_coefficient(b),
			     term_content(leading_coefficient(b)))
               b = divexact(b, d)
            else # last ditched attempt to find easy content
               h = gcd(leading_coefficient(b), glead)
               h = divexact(h, term_content(h))
               flag, q = divides(b, h)
               if flag
                  b = q
               end
            end
         end
      end
      b = divexact(b, content(b))
   end
   b = c*b
   return divexact(b, canonical_unit(leading_coefficient(b)))
end

function gcd(a::PolyElem{T}, b::PolyElem{T}) where {T <: FieldElement}
   check_parent(a, b)
   if length(a) > length(b)
      (a, b) = (b, a)
   end
   if iszero(b)
      if iszero(a)
         return(a)
      else
         d = leading_coefficient(a)
         return divexact(a, d)
      end
   end
   g = gcd(content(a), content(b))
   a = divexact(a, g)
   b = divexact(b, g)
   while !iszero(a)
      (a, b) = (mod(b, a), a)
   end
   b = g*b
   d = leading_coefficient(b)
   return divexact(b, d)
end

function content(a::PolyElem)
   z = base_ring(a)() # normalise first coefficient
   for i = 1:length(a)
      z = gcd(z, coeff(a, i - 1))
   end
   return z
end

function primpart(a::PolyElem)
   d = content(a)
   if iszero(d)
      return zero(parent(a))
   else
      return divexact(a, d)
   end
end

function gcdx(a::PolyElem{T}, b::PolyElem{T}) where {T <: FieldElement}
   check_parent(a, b)
   !isexact_type(T) && error("gcdx requires exact Bezout domain")
   if length(a) == 0
      if length(b) == 0
         return zero(parent(a)), zero(parent(a)), zero(parent(a))
      else
         d = leading_coefficient(b)
         return divexact(b, d), zero(parent(a)), divexact(one(parent(a)), d)
      end
   end
   if length(b) == 0
      d = leading_coefficient(a)
      return divexact(a, d), divexact(one(parent(a)), d), zero(parent(a))
   end
   swap = false
   if length(a) < length(b)
      a, b = b, a
      swap = true
   end
   c1 = content(a)
   c2 = content(b)
   A = divexact(a, c1)
   B = divexact(b, c2)
   u1, u2 = parent(a)(inv(c1)), zero(parent(a))
   v1, v2 = zero(parent(a)), parent(a)(inv(c2))
   while length(B) > 0
      (Q, B), A = divrem(A, B), B
      u2, u1 = u1 - Q*u2, u2
      v2, v1 = v1 - Q*v2, v2
   end
   if swap
      u1, v1 = v1, u1
   end
   d = gcd(c1, c2)
   A, u1, v1 = d*A, d*u1, d*v1
   d = leading_coefficient(A)
   return divexact(A, d), divexact(u1, d), divexact(v1, d)
end

function gcdinv(a::PolyElem{T}, b::PolyElem{T}) where {T <: FieldElement}
   check_parent(a, b)
   R = base_ring(a)
   if length(a) == 0
      if length(b) == 0
         return zero(parent(a)), zero(parent(a))
      else
         d = leading_coefficient(b)
         return divexact(b, d), zero(parent(a))
      end
   end
   if length(b) == 0
      d = leading_coefficient(a)
      return divexact(a, d), parent(a)(inv(d))
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
   d = leading_coefficient(A)
   return divexact(A, d), divexact(u1, d)
end

function PolynomialRing(R::Ring, s::Symbol; cached::Bool = true)
   return Generic.PolynomialRing(R, s; cached=cached)
end

function PolynomialRing(R::Ring, s::AbstractString; cached::Bool = true)
   return PolynomialRing(R, Symbol(s); cached=cached)
end

