###############################################################################
#
#   FunctionField.jl : Generic univariate function fields (algebraic extension
#                      of rational function field)
#
###############################################################################

export FunctionField

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
function _rat_poly(p::Poly{Rat{T}}; cached::Bool=true) where T <: FieldElement
   K = base_ring(p)
   R = base_ring(fraction_field(K))
   S = elem_type(R)

   par = PolyRing{S}(R, parent(p).S, cached)

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
#   FunctionField constructor
#
###############################################################################

function powers_precompute(p::Poly{Rat{T}}) where T <: FieldElement
   pol, d = _rat_poly(p)
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
      for i = 1:2*len - 1
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
      for i = 1:2*len - 1
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