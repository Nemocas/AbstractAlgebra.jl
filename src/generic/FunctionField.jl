###############################################################################
#
#   FunctionField.jl : Generic univariate function fields (algebraic extension
#                      of rational function field)
#
###############################################################################

export FunctionField, Frac

###############################################################################
#
#   Rational arithmetic : equivalent of _fmpq in Flint for k(x)
#
###############################################################################

function _rat_add(n1::PolyElem{T}, d1::PolyElem{T},
           n2::PolyElem{T}, d2::PolyElem{T}) where T <: FieldElement
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
   return rnum, rden
end

function _rat_sub(n1::PolyElem{T}, d1::PolyElem{T},
                  n2::PolyElem{T}, d2::PolyElem{T}) where T <: FieldElement
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
   return rnum, rden
end

function _rat_mul(n1::PolyElem{T}, d1::PolyElem{T},
                  n2::PolyElem{T}, d2::PolyElem{T}) where T <: FieldElement
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
   return n, d
end

###############################################################################
#
#   Rational polynomial arithmetic : equivalent of _fmpq_poly in Flint for k(x)
#
###############################################################################

function _rat_poly_add(poly1::Poly{S}, den1::S,
                        poly2::Poly{S}, den2::S) where
                        {T <: FieldElement, S <: PolyElem{T}}
   R = base_ring(poly1)
   
   if den1 == den2
      rpoly = poly1 + poly2

      if isone(den1)
         rden = deepcopy(den1)
      else
         d = content(rpoly)

         if !isone(d)
            d = gcd(d, den1)
         end

         if isone(d)
            rden = deepcopy(den1)
         else
            rpoly = divexact(rpoly, d)
            rden = divexact(den1, d)
         end
      end
   
      return rpoly, rden
   end

   d = gcd(den1, den2)

   if isone(d)
      rpoly = poly1*den2
      t = poly2*den1
      rpoly = addeq!(rpoly, t)
      rden = den1*den2
   else
      den11 = divexact(den1, d)
      den22 = divexact(den2, d)

      rpoly = poly1*den22
      t = poly2*den11
      rpoly = addeq!(rpoly, t)

      if iszero(rpoly)
         rden = one(R)
      else
         c = content(rpoly)

         if !isone(c)
            c = gcd(c, d)
         end

         if isone(c)
            rden = den1*den22
         else
            rpoly = divexact(rpoly, c)
            den11 = divexact(den1, c)
            rden = den11*den22
         end
      end
   end
   return rpoly, rden
end

function _rat_poly_sub(poly1::Poly{S}, den1::S,
                        poly2::Poly{S}, den2::S) where
                        {T <: FieldElement, S <: PolyElem{T}}
   R = base_ring(poly1)
   
   if den1 == den2
      rpoly = poly1 - poly2

      if isone(den1)
         rden = deepcopy(den1)
      else
         d = content(rpoly)

         if !isone(d)
            d = gcd(d, den1)
         end

         if isone(d)
            rden = deepcopy(den1)
         else
            rpoly = divexact(rpoly, d)
            rden = divexact(den1, d)
         end
      end
   
      return rpoly, rden
   end

   d = gcd(den1, den2)

   if isone(d)
      rpoly = poly1*den2 - poly2*den1
      rden = den1*den2
   else
      den11 = divexact(den1, d)
      den22 = divexact(den2, d)

      rpoly = poly1*den22 - poly2*den11

      if iszero(rpoly)
         rden = one(R)
      else
         c = content(rpoly)

         if !isone(c)
            c = gcd(c, d)
         end

         if isone(c)
            rden = den1*den22
         else
            rpoly = divexact(rpoly, c)
            den11 = divexact(den1, c)
            rden = den11*den22
         end
      end
   end
   return rpoly, rden
end

function _rat_poly_mul(poly1::Poly{S}, den1::S,
                       poly2::Poly{S}, den2::S) where
                       {T <: FieldElement, S <: PolyElem{T}}
   R = base_ring(poly1)

   if !isone(den2)
      gcd1 = content(poly1)
      gcd1 = gcd(gcd1, den2)
   else
      gcd1 = one(R)
   end
   
   if !isone(den1)
      gcd2 = content(poly2)
      gcd2 = gcd(gcd2, den1)
   else
      gcd2 = one(R)
   end

   rpoly = poly1*poly2
   rden = den1*den2

   if !isone(gcd1) || !isone(gcd2)
      g = gcd1*gcd2
      rpoly = divexact(rpoly, g)
      rden = divexact(rden, g)
   end
   return rpoly, rden
end

function _rat_poly_canonicalise(poly::Poly{S}, den::S) where
                       {T <: FieldElement, S <: PolyElem{T}}
   R = base_ring(poly)

   if isone(den)
      return poly, den
   end

   if isone(-den)
      return -poly, one(R)
   end

   if length(poly) == 0
      return poly, one(R)
   end

   g = content(poly)
   g = gcd(g, den)

   if !isone(g)
      poly = divexact(poly, g)
      den = divexact(den, g)
   end

   return poly, den
end

function _rat_poly_rem(poly1::Poly{S}, den1::S,
                       poly2::Poly{S}, den2::S) where
                       {T <: FieldElement, S <: PolyElem{T}}
   R = base_ring(poly1)

   len1 = length(poly1)
   len2 = length(poly2)

   if len1 < len2
      rpoly = deepcopy(poly1)
      rden = deepcopy(den1)
      return rpoly, rden
   end
   
   if len2 == 1
      rpoly = zero(parent(poly1))
      rden = one(R)
      return rpoly, rden
   end

   lenq = len1 - len2 + 1

   rpoly = pseudorem(poly1, poly2)

   lead = leading_coefficient(poly2)

   if isone(lead)
      rden = deepcopy(den1)
   elseif isone(-lead)
      rden = -den1
      rpoly = -rpoly
   else
      rden = den1*lead^lenq
   end

   return _rat_poly_canonicalise(rpoly, rden)
end

# convert a polynomial over a rational function field to
# a numerator and denominator
function _rat_poly(p::Poly{Rat{T}}, var=parent(p).S; cached::Bool=true) where T <: FieldElement
   K = base_ring(p)
   R = base_ring(fraction_field(K))
   S = elem_type(R)

   par = PolyRing{S}(R, var, cached)

   len = length(p)

   if len == 0
      rpol = Poly{S}(R())
      rpol.parent = par
      return rpol, R()
   end

   d = one(R)
   for i = 1:len
      d = lcm(d, denominator(coeff(p, i - 1)))
   end

   V = Vector{S}(undef, len)
   for i = 1:len
      c = coeff(p, i - 1)
      den_i = denominator(c)
      if den_i == d
         V[i] = deepcopy(numerator(c))
      else
         V[i] = divexact(d, den_i)*numerator(c)
      end
   end

   rpol = Poly{S}(V)
   rpol.parent = par

   return rpol, d
end

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{FunctionFieldElem{T}}) where T <: FieldElement = FunctionField{T}

elem_type(::Type{FunctionField{T}}) where T <: FieldElement = FunctionFieldElem{T}

function base_ring(a::FunctionField{T}) where T <: FieldElement
   return a.base_ring::RationalFunctionField{T}
end

base_ring(a::FunctionFieldElem) = base_ring(parent(a))

parent(a::FunctionFieldElem) = a.parent

function isexact_type(a::Type{T}) where {S <: FieldElement, T <: FunctionFieldElem{S}}
   return isexact_type(S)
end

var(a::FunctionField) = a.S

function check_parent(a::FunctionFieldElem{T}, b::FunctionFieldElem{T}, throw::Bool = true) where T <: FieldElement
   fl = parent(a) != parent(b)
   fl && throw && error("Incompatible function fields in function field operation")
   return !fl
end

characteristic(a::FunctionField) = characteristic(base_ring(a))

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.numerator(a::FunctionFieldElem, canonicalise::Bool=true)
   if canonicalise
      u = canonical_unit(a.den)
      return divexact(a.num, u)
   else
      return a.num
   end
end

function Base.denominator(a::FunctionFieldElem, canonicalise::Bool=true)
   if canonicalise
      u = canonical_unit(a.den)
      return divexact(a.den, u)
   else
      return a.den
   end
end

degree(R::FunctionField) = degree(R.num)

zero(R::FunctionField) = R()

one(R::FunctionField) = R(1)

function gen(R::FunctionField{T}) where T <: FieldElement
   return FunctionFieldElem{T}(R,
                              deepcopy(R.powers[2]), deepcopy(R.powers_den[2]))
end

iszero(a::FunctionFieldElem) = iszero(a.num)

isone(a::FunctionFieldElem) = isone(a.num) && isone(a.den)

isunit(a::FunctionFieldElem) = !iszero(a)

isgen(a::FunctionFieldElem) = isgen(a.num) && isone(a.den)

function coeff(a::FunctionFieldElem, n::Int)
   R = base_ring(a)
   n = coeff(a.num, n)
   d = a.den
   return R(n//d)
end

function num_coeff(a::FunctionFieldElem, n::Int)
   return coeff(a.num, n)
end

function deepcopy_internal(a::FunctionFieldElem, dict::IdDict)
   R = parent(a)
   return R(deepcopy_internal(a.num, dict), deepcopy_internal(a.den, dict))
end

function Base.hash(a::FunctionFieldElem, h::UInt)
   b = 0x52fd76bf2694aa02%UInt
   b = xor(hash(a.den, h), xor(hash(a.num, h), b))
   return b
end

function _rat_poly(a::FunctionFieldElem)
   return numerator(a), denominator(a)
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function AbstractAlgebra.expressify(a::FunctionFieldElem; context = nothing)
   n = numerator(a, true)
   d = denominator(a, true)
   if isone(d)
       return expressify(n)
   else
       return Expr(:call, ://, expressify(n), expressify(d))
   end
end

function show(io::IO, ::MIME"text/plain", a::FunctionFieldElem)
   print(io, AbstractAlgebra.obj_to_string(a, context = io))
 end
 
 function show(io::IO, a::FunctionFieldElem)
   print(io, AbstractAlgebra.obj_to_string(a, context = io))
 end

function show(io::IO, a::FunctionField)
   print(IOContext(io, :compact => true), "Function Field over ",
         base_ring(base_ring(a)), " with defining polynomial ", a.num)
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::FunctionFieldElem)
   R = parent(a)
   return R(-numerator(a), denominator(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::FunctionFieldElem{T}, b::FunctionFieldElem{T}) where T <: FieldElement
   check_parent(a, b)
   R = parent(a)
   n1, d1 = _rat_poly(a)
   n2, d2 = _rat_poly(b)
   return R(_rat_poly_add(n1, d1, n2, d2)...)
end

function -(a::FunctionFieldElem{T}, b::FunctionFieldElem{T}) where T <: FieldElement
   check_parent(a, b)
   R = parent(a)
   n1, d1 = _rat_poly(a)
   n2, d2 = _rat_poly(b)
   return R(_rat_poly_sub(n1, d1, n2, d2)...)
end

function *(a::FunctionFieldElem{T}, b::FunctionFieldElem{T}) where T <: FieldElement
   check_parent(a, b)
   R = parent(a)
   n1, d1 = _rat_poly(a)
   n2, d2 = _rat_poly(b)
   z = R(_rat_poly_mul(n1, d1, n2, d2)...)
   return reduce!(z)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::FunctionFieldElem, b::Union{Integer, Rational})
   R = parent(a)
   num = a.num*b
   return R(_rat_poly_canonicalise(num, a.den)...)
end

*(a::Union{Integer, Rational}, b::FunctionFieldElem) = b*a

function *(a::FunctionFieldElem{T}, b::T) where T <: FieldElem
   R = parent(a)
   num = a.num*b
   return R(_rat_poly_canonicalise(num, a.den)...)
end

*(a::T, b::FunctionFieldElem{T}) where T <: FieldElem = b*a

function *(a::FunctionFieldElem{T}, b::Rat{T}) where T <: FieldElement
   parent(b) != base_ring(a) && error("Could not coerce element")
   R = parent(a)
   num = a.num*numerator(b, false)
   den = a.den*denominator(b, false)
   return R(_rat_poly_canonicalise(num, den)...)
end

*(a::Rat{T}, b::FunctionFieldElem{T}) where T <: FieldElement = b*a

*(a::FunctionFieldElem, b::RingElem) = a*base_ring(a)(b)

*(a::RingElem, b::FunctionFieldElem) = b*a

function +(a::FunctionFieldElem{T}, b::Rat{T}) where T <: FieldElement
   parent(b) != base_ring(a) && error("Unable to coerce element")
   return a + parent(a)(b)
end

+(a::Rat{T}, b::FunctionFieldElem{T}) where T <: FieldElement = b + a

+(a::FunctionFieldElem, b::Union{Integer, Rational}) = a + base_ring(a)(b)

+(a::Union{Integer, Rational}, b::FunctionFieldElem) = b + a

+(a::FunctionFieldElem, b::RingElem) = a + base_ring(a)(b)

+(a::RingElem, b::FunctionFieldElem) = b + a

function -(a::FunctionFieldElem{T}, b::Rat{T}) where T <: FieldElement
   parent(b) != base_ring(a) && error("Unable to coerce element")
   return a - parent(a)(b)
end

function -(a::Rat{T}, b::FunctionFieldElem{T}) where T <: FieldElement
   parent(a) != base_ring(b) && error("Unable to coerce element")
   return parent(b)(a) - b
end

-(a::FunctionFieldElem, b::Union{Integer, Rational}) = a - base_ring(a)(b)

-(a::Union{Integer, Rational}, b::FunctionFieldElem) = base_ring(b)(a) - b

-(a::FunctionFieldElem, b::RingElem) = a - base_ring(a)(b)

-(a::RingElem, b::FunctionFieldElem) = base_ring(b)(a) - b

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::FunctionFieldElem{T}, b::Int) where T <: FieldElement
   b < 0 && error("Not implemented")
   R = parent(a)
   if isgen(a) && b < 2*length(R.num) - 3 # special case powers of generator
      return R(deepcopy(R.powers[b + 1]), deepcopy(R.powers_den[b + 1]))
   elseif b == 0
      return one(R)
   elseif iszero(a)
      return zero(R)
   elseif length(a.num) == 1
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

###############################################################################
#
#   Unsafe operators
#
###############################################################################

function zero!(a::FunctionFieldElem)
   a.num = zero!(a.num)
   R = parent(a.den)
   a.den = one(R)
   return a
end

function setcoeff!(a::FunctionFieldElem{T}, n::Int, c::Rat{T}) where T <: FieldElement
   base_ring(a) != parent(c) && error("Unable to coerce element")
   n < 0 || n > degree(parent(a)) && error("Degree not in range")
   cnum = numerator(c.d, false)
   cden = denominator(c.d, false)
   anum = a.num
   aden = a.den
   g = gcd(cden, aden)
   if g != aden
      u = divexact(aden, g)
      cnum *= u
   end
   if g != cden
      u = divexact(cden, g)
      anum *= u
      aden *= u
   end
   anum = setcoeff!(anum, n, cnum)
   a.num, a.den = _rat_poly_canonicalise(anum, aden)
   return a
end

function setcoeff!(a::FunctionFieldElem{T}, n::Int, c::PolyElem{T}) where T <: FieldElement
   parent(c) != parent(a.den) && error("Unable to coerce element")
   n < 0 || n > degree(parent(a)) && error("Degree not in range")
   if !isone(a.den)
      c *= a.den
   end
   anum = setcoeff!(a.num, n, c)
   a.num, a.den = _rat_poly_canonicalise(anum, a.den)
   return a
end

function setcoeff!(a::FunctionFieldElem, n::Int, c::RingElement)
   return setcoeff!(a, n, base_ring(a)(c))
end

function reduce!(a::FunctionFieldElem)
   R = parent(a)
   len = length(R.num)
   num = numerator(a)
   den = denominator(a)
   z = truncate(num, len - 1)
   zden = den
   for i = len:length(numerator(a))
      c = coeff(num, i - 1)
      z, zden = _rat_poly_add(z, zden, c*R.powers[i], R.powers_den[i])
   end
   a.num, a.den = _rat_poly_canonicalise(z, zden)
   return a
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (R::FunctionField{T})(p::Poly{S}, den::S) where
                                          {T <: FieldElement, S <: PolyElem{T}}
   z = FunctionFieldElem{T}(R, p, den)
   return z
end

function (R::FunctionField{T})() where T <: FieldElement
   p = zero(parent(R.powers[1]))
   den = one(parent(R.powers_den[1]))
   z = FunctionFieldElem{T}(R, p, den)
   return z
end

function (R::FunctionField{T})(a::Union{Rational, Integer}) where T <: FieldElement
   p = parent(R.powers[1])(a)
   den = one(parent(R.powers_den[1]))
   z = FunctionFieldElem{T}(R, p, den)
   return z
end

function (R::FunctionField{T})(a::T) where T <: FieldElem
   p = parent(R.powers[1])(a)
   den = one(parent(R.powers_den[1]))
   z = FunctionFieldElem{T}(R, p, den)
   return z
end

function (R::FunctionField{T})(a::Rat{T}) where T <: FieldElement
   p = parent(R.powers[1])([numerator(a)])
   den = parent(R.powers_den[1])(denominator(a))
   z = FunctionFieldElem{T}(R, p, den)
   return z
end

(R::FunctionField)(b::RingElem) = R(base_ring(R)(b))

###############################################################################
#
#   FunctionField constructor
#
###############################################################################

function powers_precompute(pol::Poly{W}, d::W) where {T <: FieldElement, W <: PolyElem{T}}
   len = length(pol)
   R = parent(d)
   S = parent(pol)
   U = typeof(pol)
   V = typeof(d)
   P = Vector{U}(undef, 2*len - 1)
   monic = isone(coeff(pol, len - 1))
   D = Vector{V}(undef, 2*len - 1)
   if monic
      pow = one(S)
      for i = 1:2*len - 3
         D[i] = one(R)
         if length(pow) == len
            pow -= pol*coeff(pow, len - 1)
            pow = truncate(pow, len - 1)
         end
         P[i] = pow
         pow = shift_left(pow, 1)
      end
   else
      pow = one(S)
      den = one(R)
      for i = 1:2*len - 3
         if length(pow) == len
            tden = den*coeff(pol, len - 1)
            t = pol*coeff(pow, len - 1)
            t = truncate(t, len - 1)
            t, tden = _rat_poly_canonicalise(t, tden)
            pow, den = _rat_poly_sub(pow, den, t, tden)
         end
         P[i] = pow
         D[i] = den
         pow = shift_left(pow, 1)
      end
   end
   return monic, P, D
end

function FunctionField(p::Poly{Rat{T}}, s::AbstractString; cached::Bool=true) where T <: FieldElement
   sym = Symbol(s)
   pol, den = _rat_poly(p, sym)
   
   par = FunctionField{T}(pol, den, sym, cached)
   par.monic, par.powers, par.powers_den = powers_precompute(pol, den)
   par.base_ring = base_ring(p)

   return par, gen(par)
end
