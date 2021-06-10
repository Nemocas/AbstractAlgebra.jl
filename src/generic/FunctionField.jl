###############################################################################
#
#   FunctionField.jl : Generic univariate function fields (algebraic extension
#                      of rational function field)
#
###############################################################################

export FunctionField, base_field

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

function _rat_canonicalise(n::PolyElem{T}, d::PolyElem{T}) where
                                                              T <: FieldElement
   g = gcd(n, d)
   if !isone(g)
      n = divexact(n, g)
      d = divexact(d, g)
   end
   u = canonical_unit(d)
   if !isone(u)
      n = divexact(n, u)
      d = divexact(d, u)
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

function _rat_poly_gcdx(a::Poly{U}, den_a::U,
                         b::Poly{U}, den_b::U) where
                         {T <: FieldElement, U <: PolyElem{T}}
   S = parent(a)
   R = parent(den_a)

   lena = length(a)
   lenb = length(b)

   if lena == 0 && lenb == 0
      return zero(S), one(R), zero(S), one(R), zero(S), one(R)
   end

   if lena == 0
      c = leading_coefficient(b)
      g, den_g = _rat_poly_canonicalise(b, c)
      t, den_t = _rat_poly_canonicalise(S(den_b), c)
      return g, den_g, zero(S), one(R), t, den_t
   end

   if lenb == 0
      c = leading_coefficient(a)
      g, den_g = _rat_poly_canonicalise(a, c)
      s, den_s = _rat_poly_canonicalise(S(den_a), c)
      return g, den_g, s, den_s, zero(S), one(R)
   end

   if lena == 1
      c = leading_coefficient(a)
      s, den_s = _rat_poly_canonicalise(S(den_a), c)
      return one(S), one(R), s, den_s, zero(S), one(R)
   end

   if lenb == 1
      c = leading_coefficient(b)
      t, den_t = _rat_poly_canonicalise(S(den_b), c)
      return one(S), one(R), zero(S), one(R), t, den_t
   end

   ca = content(a)
   cb = content(b)

   if !isone(ca)
      a = divexact(a, ca)
   end
   if !isone(cb)
      b = divexact(b, cb)
   end

   g = gcd(a, b)
   if length(g) > 1
      a = divexact(a, g)
      b = divexact(b, g)
   end

   den_g, s, t = resx(a, b)

   den_g *= leading_coefficient(g)

   s *= den_a
   den_s = ca*den_g

   t *= den_b
   den_t = cb*den_g

   s, den_s = _rat_poly_canonicalise(s, den_s)
   t, den_t = _rat_poly_canonicalise(t, den_t)

   return g, leading_coefficient(g), s, den_s, t, den_t   
end

function _rat_poly_resultant(a::Poly{U}, den_a::U,
                                    b::Poly{U}, den_b::U) where
                                          {T <: FieldElement, U <: PolyElem{T}}
   lena = length(a)
   lenb = length(b)
   R = parent(den_a)
   if lena == 0 || lenb == 0
      return zero(R), one(R)
   end
   ca = content(a)
   cb = content(b)
   a = divexact(a, ca)
   b = divexact(b, cb)
   g = gcd(a, b)
   if length(g) > 1 # not coprime, resultant is zero
      return zero(R), one(R)
   end
   rnum = resultant(a, b)
   if !isone(ca)
      rnum *= ca^(lenb - 1)
   end
   if !isone(cb)
      rnum *= cb^(lena - 1)
   end
   rden = den_a^(lenb - 1)*den_b^(lena - 1)
   return _rat_canonicalise(rnum, rden) 
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
      d = lcm(d, denominator(coeff(p, i - 1), false))
   end

   V = Vector{S}(undef, len)
   for i = 1:len
      c = coeff(p, i - 1)
      den_i = denominator(c, false)
      if den_i == d
         V[i] = deepcopy(numerator(c, false))
      else
         V[i] = divexact(d, den_i)*numerator(c, false)
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

function base_ring(R::FunctionField{T}) where T <: FieldElement
   return R.base_ring::RationalFunctionField{T}
end

base_ring(R::FunctionFieldElem) = base_ring(parent(R))

# For consistency with number fields in Hecke.jl
@doc Markdown.doc"""
    base_field(R::FunctionField)

Return the rational function field that the field `R` is an extension of.
Synonymous with `base_ring`.
"""
base_field(R::FunctionField) = base_ring(R)

parent(a::FunctionFieldElem) = a.parent

function isexact_type(a::Type{T}) where {S <: FieldElement, T <: FunctionFieldElem{S}}
   return isexact_type(S)
end

@doc Markdown.doc"""
    var(R::FunctionField)

Return the variable name of the generator of the function field `R` as a
symbol.
"""
var(R::FunctionField) = R.S

function check_parent(a::FunctionFieldElem{T}, b::FunctionFieldElem{T}, throw::Bool = true) where T <: FieldElement
   fl = parent(a) != parent(b)
   fl && throw && error("Incompatible function fields in function field operation")
   return !fl
end

@doc Markdown.doc"""
    characteristic(R::FunctionField)

Return the characteristic of the underlying rational function field.
"""
characteristic(R::FunctionField) = characteristic(base_ring(R))

###############################################################################
#
#   Basic manipulation
#
###############################################################################

@doc Markdown.doc"""
    defining_polynomial(R::FunctionField)
    modulus(R::FunctionField)

Return the original polynomial that was used to define the function field `R`.
"""
defining_polynomial(R::FunctionField) = R.pol

modulus(R::FunctionField) = defining_polynomial(R)

function power_precomp(R::FunctionField{T}, n::Int) where T <: FieldElement
   return R.powers[n + 1]::Poly{dense_poly_type(T)}
end

function power_precomp_den(R::FunctionField{T}, n::Int) where T <: FieldElement
   return R.powers_den[n + 1]::dense_poly_type(T)
end

function trace_precomp(R::FunctionField{T}, n::Int) where T <: FieldElement
   return R.traces[n + 1]::dense_poly_type(T)
end

function trace_precomp_den(R::FunctionField{T}) where T <: FieldElement
   return R.traces_den::dense_poly_type(T)
end

@doc Markdown.doc"""
    Base.numerator(R::FunctionField{T}, canonicalise::Bool=true) where T <: FieldElement
    Base.denominator(R::FunctionField{T}, canonicalise::Bool=true) where T <: FieldElement

Thinking of elements of the rational function field as fractions, put the
defining polynomial of the function field over a common denominator
and return the numerator/denominator respectively. Note that the resulting
polynomials belong to a different ring than the original defining polynomial.
The `canonicalise` is ignored, but exists for compatibility with the Generic
interface.
"""
function Base.numerator(R::FunctionField{T},
                               canonicalise::Bool=true) where T <: FieldElement
   # only used for type assert, so no need to canonicalise
   return R.num::Poly{dense_poly_type(T)}
end                 

function Base.denominator(R::FunctionField{T},
                               canonicalise::Bool=true) where T <: FieldElement
   # only used for type assert, so no need to canonicalise
   return R.den::dense_poly_type(T)
end                 

@doc Markdown.doc"""
    Base.numerator(a::FunctionFieldElem{T}, canonicalise::Bool=true) where T <: FieldElement
    Base.denominator(a::FunctionFieldElem{T}, canonicalise::Bool=true) where T <: FieldElement

Return the numerator and denominator of the function field element `a`.
Note that elements are stored in fraction free form so that the denominator
is a common denominator for the coefficients of the element `a`.
If `canonicalise` is set to `true` the fraction is first canonicalised.
"""
function Base.numerator(a::FunctionFieldElem{T},
                               canonicalise::Bool=true) where T <: FieldElement
   anum = a.num::Poly{dense_poly_type(T)}
   aden = a.den::dense_poly_type(T)
   if canonicalise
      u = canonical_unit(aden)
      return divexact(anum, u)
   else
      return anum
   end
end

function Base.denominator(a::FunctionFieldElem{T},
                               canonicalise::Bool=true) where T <: FieldElement
   aden = a.den::dense_poly_type(T)
   if canonicalise
      u = canonical_unit(aden)
      return divexact(aden, u)
   else
      return aden
   end
end

@doc Markdown.doc"""
    degree(S::FunctionField)

Return the degree of the defining polynomial of the function field, i.e. the
degree of the extension that the function field makes of the underlying
rational function field.
"""
degree(S::FunctionField) = degree(numerator(S))

zero(S::FunctionField) = S()

one(S::FunctionField) = S(1)

@doc Markdown.doc"""
    gen(S::FunctionField{T}) where T <: FieldElement

Return the generator of the function field returned by the function field
constructor.
"""
function gen(S::FunctionField{T}) where T <: FieldElement
   if degree(S) == 1
      return S(-coeff(modulus(S), 0)//coeff(modulus(S), 1))
   else
      return FunctionFieldElem{T}(S, deepcopy(power_precomp(S, 1)),
                                     deepcopy(power_precomp_den(S, 1)))
   end
end

iszero(a::FunctionFieldElem) = iszero(numerator(a, false))

isone(a::FunctionFieldElem) = numerator(a, false) == denominator(a, false)

isunit(a::FunctionFieldElem) = !iszero(a)

@doc Markdown.doc"""
    isgen(a::FunctionFieldElem)

Return `true` if `a` is the generator of the function field returned by the
function field constructor.
"""
function isgen(a::FunctionFieldElem)
   S = parent(a)
   if degree(S) == 1
      return a == S(-coeff(modulus(S), 0)//coeff(modulus(S), 1))
   else
      return isgen(numerator(a, false)) && isone(denominator(a, false)) ||
             isgen(numerator(a, true)) && isone(denominator(a, true))
   end
end

@doc Markdown.doc"""
    coeff(a::FunctionFieldElem, n::Int)

Return the degree `n` coefficient of the element `a` in its polynomial
representation in terms of the generator of the function field. The
coefficient is returned as an element of the underlying rational function
field.
"""
function coeff(a::FunctionFieldElem, n::Int)
   R = base_ring(a)
   n = coeff(numerator(a, false), n)
   d = denominator(a, false)
   return R(n//d)
end

@doc Markdown.doc"""
    num_coeff(a::FunctionFieldElem, n::Int)

Return the degree `n` coefficient of the numerator of the element `a` (in its
polynomial representation in terms of the generator of the function field,
rationalised as per `numerator/denominator` described above). The coefficient
will be an polynomial over the `base_ring` of the underlying rational function
field.
"""
function num_coeff(a::FunctionFieldElem, n::Int)
   return coeff(numerator(a, false), n)
end

function deepcopy_internal(a::FunctionFieldElem, dict::IdDict)
   S = parent(a)
   return S(deepcopy_internal(numerator(a, false), dict),
            deepcopy_internal(denominator(a, false), dict))
end

function Base.hash(a::FunctionFieldElem, h::UInt)
   b = 0x52fd76bf2694aa02%UInt
   b = xor(hash(denominator(a, false), h),
                                          xor(hash(numerator(a, false), h), b))
   return b
end

function _rat_poly(a::FunctionFieldElem)
   return numerator(a, false), denominator(a, false)
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

function show(io::IO, R::FunctionField)
   print(IOContext(io, :compact => true), "Function Field over ",
         base_ring(base_ring(R)), " with defining polynomial ",
         numerator(R))
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::FunctionFieldElem)
   R = parent(a)
   return R(-numerator(a, false), denominator(a, false))
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
   num = numerator(a, false)*b
   return R(_rat_poly_canonicalise(num, denominator(a, false))...)
end

*(a::Union{Integer, Rational}, b::FunctionFieldElem) = b*a

function *(a::FunctionFieldElem{T}, b::T) where T <: FieldElem
   R = parent(a)
   num = numerator(a, false)*b
   return R(_rat_poly_canonicalise(num, denominator(a, false))...)
end

*(a::T, b::FunctionFieldElem{T}) where T <: FieldElem = b*a

function *(a::FunctionFieldElem{T}, b::Rat{T}) where T <: FieldElement
   parent(b) != base_ring(a) && error("Could not coerce element")
   R = parent(a)
   num = numerator(a, false)*numerator(b, false)
   den = denominator(a, false)*denominator(b, false)
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
   if isgen(a) && b < 2*length(numerator(R)) - 3 # special case powers of generator
      return R(deepcopy(power_precomp(R, b)),
               deepcopy(power_precomp_den(R, b)))
   elseif b == 0
      return one(R)
   elseif iszero(a)
      return zero(R)
   elseif length(numerator(a, false)) == 1
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
#   Comparison
#
###############################################################################

function ==(a::FunctionFieldElem{T}, b::FunctionFieldElem{T}) where T <: FieldElement
   check_parent(a, b)
   aden = denominator(a, true)
   bden = denominator(b, true)
   if aden != bden
      return false
   end
   anum = numerator(a, true)
   bnum = numerator(b, true)
   return anum == bnum
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

function ==(a::FunctionFieldElem{T}, b::Rat{T}) where T <: FieldElement
   parent(b) != base_ring(a) && error("Unable to coerce element")
   if iszero(a) && iszero(b)
      return true
   elseif length(numerator(a, false)) != 1
      return false
   end
   return a == parent(a)(b)
end

==(a::Rat{T}, b::FunctionFieldElem{T}) where T <: FieldElement = b == a

==(a::FunctionFieldElem, b::Union{Integer, Rational}) = a == base_ring(a)(b)

==(a::Union{Integer, Rational}, b::FunctionFieldElem) = b == a

==(a::FunctionFieldElem, b::RingElem) = a == base_ring(a)(b)

==(a::RingElem, b::FunctionFieldElem) = b == a

###############################################################################
#
#   Inversion
#
###############################################################################

function Base.inv(a::FunctionFieldElem)
   R = parent(a)   
   anum = numerator(a, false)
   aden = denominator(a, false)

   G, G_den, S, S_den, T, T_den =
                   _rat_poly_gcdx(anum, aden, numerator(R), denominator(R))
   return R(S, S_den)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::FunctionFieldElem{T}, b::FunctionFieldElem{T}) where
                                                              T <: FieldElement
   return a*inv(b)
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(a::FunctionFieldElem, b::Union{Rational, Integer})
   S = parent(a)
   anum = numerator(a, false)
   aden = denominator(a, false)
   R = parent(aden)
   rnum, rden = _rat_poly_canonicalise(anum, R(b))
   return S(rnum, rden*aden)
end

function divexact(a::FunctionFieldElem{T}, b::T) where T <: FieldElem
   S = parent(a)
   anum = numerator(a, false)
   aden = denominator(a, false)
   R = parent(aden)
   rnum, rden = _rat_poly_canonicalise(anum, R(b))
   return S(rnum, rden*aden)
end

function divexact(a::FunctionFieldElem{T}, b::Rat{T}) where T <: FieldElement
   S = parent(a)
   base_ring(a) != parent(b) && error("Incompatible fields")
   bnum = numerator(b, false)
   bden = denominator(b, false)
   anum = numerator(a, false)
   aden = denominator(a, false)
   return S(_rat_poly_canonicalise(anum*bden, aden*bnum)...)
end

divexact(a::FunctionFieldElem, b::RingElem) = divexact(a, base_ring(a)(b))

divexact(a::Union{Rational, Integer}, b::FunctionFieldElem) = a*inv(b)

divexact(a::T, b::FunctionFieldElem{T}) where T <: FieldElement = a*inv(b)

divexact(a::Rat{T}, b::FunctionFieldElem{T}) where T <: FieldElement = a*inv(b)

divexact(a::RingElem, b::FunctionFieldElem) = a*inv(b)

###############################################################################
#
#   Norm
#
###############################################################################

@doc Markdown.doc"""
    norm(a::FunctionFieldElem)

Return the absolute norm of `a` as an element of the underlying rational
function field.
"""
function norm(a::FunctionFieldElem)
   R = base_ring(a)
   S = parent(a)
   if iszero(a)
      return zero(R)
   end
   anum = numerator(a, false)
   aden = denominator(a, false)
   alen = length(anum)
   pol = numerator(S, false)
   rnum, rden = _rat_poly_resultant(pol, one(base_ring(pol)), anum, aden)
   if !S.monic && length(anum) > 1
      lc = leading_coefficient(pol)
      pow = lc^(alen - 1)
      rden *= pow
      if alen == 2 || (alen > 2 && !isone(gcd(lc, rnum))) # must rationalise
         rnum, rden = _rat_canonicalise(rnum, rden)
      end
   end
   return R(rnum, rden)
end

###############################################################################
#
#   Trace
#
###############################################################################

@doc Markdown.doc"""
    tr(a::FunctionFieldElem)

Return the absolute trace of `a` as an element of the underlying rational
function field.
"""
function tr(a::FunctionFieldElem)
   R = base_ring(a)
   S = parent(a)
   if iszero(a)
      return zero(R)
   end
   anum = numerator(a, false)
   aden = denominator(a, false)
   alen = length(anum)
   rnum = coeff(anum, 0)*trace_precomp(S, 0)
   for i = 1:alen - 1
      rnum += coeff(anum, i)*trace_precomp(S, i)
   end
   rden = aden*trace_precomp_den(S)
   return R(rnum, rden)
end

###############################################################################
#
#   Unsafe operators
#
###############################################################################

function zero!(a::FunctionFieldElem)
   a.num = zero!(numerator(a, false))
   R = parent(denominator(a, false))
   a.den = one(R)
   return a
end

function setcoeff!(a::FunctionFieldElem{T}, n::Int, c::Rat{T}) where T <: FieldElement
   base_ring(a) != parent(c) && error("Unable to coerce element")
   n < 0 || n > degree(parent(a)) && error("Degree not in range")
   cnum = numerator(c.d, false)
   cden = denominator(c.d, false)
   anum = numerator(a, false)
   aden = denominator(a, false)
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
   parent(c) != parent(denominator(a, false)) && error("Unable to coerce element")
   n < 0 || n > degree(parent(a)) && error("Degree not in range")
   aden = denominator(a, false)
   if !isone(aden)
      c *= aden
   end
   anum = setcoeff!(numerator(a, false), n, c)
   a.num, a.den = _rat_poly_canonicalise(anum, aden)
   return a
end

function setcoeff!(a::FunctionFieldElem, n::Int, c::RingElement)
   return setcoeff!(a, n, base_ring(a)(c))
end

function reduce!(a::FunctionFieldElem{T}) where T <: FieldElement
   R = parent(a)
   len = length(numerator(R))
   num = numerator(a, false)
   den = denominator(a, false)
   z = truncate(num, len - 1)
   zden = den
   z, zden = _rat_poly_canonicalise(z, zden)
   for i = len:length(num)
      c = coeff(num, i - 1)
      t, tden = _rat_poly_canonicalise(power_precomp(R, i - 1)*c,
                                       power_precomp_den(R, i - 1)*den)
      z, zden = _rat_poly_add(z, zden, t, tden)
   end
   a.num, a.den = _rat_poly_canonicalise(z, zden)
   return a
end

function add!(c::FunctionFieldElem{T},
      a::FunctionFieldElem{T}, b::FunctionFieldElem{T}) where T <: FieldElement
   n1, d1 = _rat_poly(a)
   n2, d2 = _rat_poly(b)
   c.num, c.den = _rat_poly_add(n1, d1, n2, d2)
   return c
end

function addeq!(c::FunctionFieldElem{T}, a::FunctionFieldElem{T}) where
                                                              T <: FieldElement
   n1, d1 = _rat_poly(c)
   n2, d2 = _rat_poly(a)
   c.num, c.den = _rat_poly_add(n1, d1, n2, d2)
   return c
end

function mul!(c::FunctionFieldElem{T},
      a::FunctionFieldElem{T}, b::FunctionFieldElem{T}) where T <: FieldElement
   n1, d1 = _rat_poly(a)
   n2, d2 = _rat_poly(b)
   c.num, c.den = _rat_poly_mul(n1, d1, n2, d2)
   return reduce!(c)
end

###############################################################################
#
#   Random generation
#
###############################################################################

RandomExtensions.maketype(K::FunctionField, _) = elem_type(K)

function RandomExtensions.make(S::FunctionField, vs...)
   R = parent(numerator(S))
   n = degree(numerator(S))
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1]) # forward to default Make constructor
   else
      make(S, make(R, (n-1):(n-1), vs...))
   end
end

function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make2{<:FunctionFieldElem,
                                                                <:FunctionField}})
   K, v = sp[][1:end]
   r = v[3]
   S = parent(numerator(K))
   R = parent(denominator(K))
   return K(_rat_poly_canonicalise(rand(rng, v), rand(rng, r))...)
end

rand(rng::AbstractRNG, K::FunctionField, v...) = rand(rng, make(K, v...))

rand(K::FunctionField, v...) = rand(Random.GLOBAL_RNG, K, v...)

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{FunctionFieldElem{T}}, ::Type{FunctionFieldElem{T}}) where T <: FieldElement = FunctionFieldElem{T}

function promote_rule(::Type{FunctionFieldElem{T}}, ::Type{U}) where {T <: FieldElem, U <: RingElem}
   promote_rule(T, U) == T ? FunctionFieldElem{T} : Union{}
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
   p = zero(parent(power_precomp(R, 0)))
   den = one(parent(power_precomp_den(R, 0)))
   z = FunctionFieldElem{T}(R, p, den)
   return z
end

function (R::FunctionField{T})(a::Union{Rational, Integer}) where T <: FieldElement
   p = parent(power_precomp(R, 0))(a)
   den = one(parent(power_precomp_den(R, 0)))
   z = FunctionFieldElem{T}(R, p, den)
   return z
end

function (R::FunctionField{T})(a::T) where T <: FieldElem
   p = parent(power_precomp(R, 0))(a)
   den = one(parent(power_precomp_den(R, 0)))
   z = FunctionFieldElem{T}(R, p, den)
   return z
end

function (R::FunctionField{T})(a::Rat{T}) where T <: FieldElement
   p = parent(power_precomp(R, 0))([numerator(a, false)])
   den = parent(power_precomp_den(R, 0))(denominator(a, false))
   z = FunctionFieldElem{T}(R, p, den)
   return z
end

function (R::FunctionField{T})(a::FunctionFieldElem) where T <: FieldElement
   parent(a) !== R && error("Unable to coerce element")
   return a
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
   P = Vector{U}(undef, 2*len - 3)
   monic = isone(coeff(pol, len - 1))
   D = Vector{V}(undef, 2*len - 3)
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
            pow = truncate(pow, len - 1)
            pow, den = _rat_poly_sub(pow, den, t, tden)
         end
         P[i] = pow
         D[i] = den
         pow = shift_left(pow, 1)
      end
   end
   return monic, P, D
end

function traces_precompute(pol::Poly{W}, d::W) where {T <: FieldElement, W <: PolyElem{T}}
   len = length(pol)
   R = parent(d)
   S = parent(pol)
   U = typeof(pol)
   V = typeof(d)
   P = Vector{V}(undef, len - 1)
   Pden = one(R)
   lc = leading_coefficient(pol)
   for i = 1:len - 2
      P[i + 1] = coeff(pol, len - i - 1)*i
      for j = i - 1:-1:1
         P[i + 1] *= lc
         P[i + 1] += coeff(pol, len - j - 1)*P[i - j + 1]
      end
      P[i + 1] = -P[i + 1]
   end
   for i = 1:len - 2
      P[len - i] *= Pden
      Pden *= lc
   end
   P[1] = Pden*(len - 1)
   return P, Pden
end

function FunctionField(p::Poly{Rat{T}}, s::Symbol; cached::Bool=true) where T <: FieldElement
   length(p) < 2 && error("Polynomial must have degree at least 1")
   pol, den = _rat_poly(p, s)
   
   par = FunctionField{T}(pol, den, s, cached)
   par.monic, par.powers, par.powers_den = powers_precompute(pol, den)
   par.traces, par.traces_den = traces_precompute(pol, den)
   par.base_ring = base_ring(p)
   par.pol = p

   return par, gen(par)
end
