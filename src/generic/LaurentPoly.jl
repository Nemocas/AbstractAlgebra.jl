###############################################################################
#
#   LaurentPoly.jl : Generic Laurent polynomials over rings
#
###############################################################################

import AbstractAlgebra: monomials_degrees
using AbstractAlgebra: monomial_degree

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{LaurentPolyWrap{T, PE}}) where {T, PE} =
   LaurentPolyWrapRing{T, parent_type(PE)}

elem_type(::Type{LaurentPolyWrapRing{T, PR}}) where {T, PR} =
   LaurentPolyWrap{T, elem_type(PR)}

parent(p::LaurentPolyWrap) = LaurentPolyWrapRing(parent(p.poly))

base_ring(R::LaurentPolyWrapRing) = base_ring(R.polyring)

var(R::LaurentPolyWrapRing) = var(R.polyring)

symbols(R::LaurentPolyWrapRing) = symbols(R.polyring)

nvars(R::LaurentPolyWrapRing) = nvars(R.polyring)

characteristic(R::LaurentPolyWrapRing) = characteristic(R.polyring)


###############################################################################
#
#   Basic manipulation
#
###############################################################################

monomials_degrees(p::LaurentPolyWrap) = p.mindeg .+ (0:degree(p.poly))

coeff(p::LaurentPolyWrap, i::Int) =
   i < p.mindeg ? zero(base_ring(p)) : coeff(p.poly, i - p.mindeg)

###############################################################################
#
#   Comparisons
#
###############################################################################

==(p::LaurentPolyWrap, q::LaurentPolyWrap) =
   all(i -> coeff(p, i) == coeff(q, i), monomials_degrees(p, q))

###############################################################################
#
#   Unary and Binary operations
#
###############################################################################

-(p::LaurentPolyWrap) = LaurentPolyWrap(-p.poly, p.mindeg)

function +(p::LaurentPolyWrap, q::LaurentPolyWrap)
   if p.mindeg > q.mindeg
      p, q = q, p
   end
   p_, q_ = p.poly, q.poly
   if p.mindeg < q.mindeg
      q_ = shift_left(q_, q.mindeg - p.mindeg)
   end
   LaurentPolyWrap(p_ + q_, p.mindeg)
end

-(p::LaurentPolyWrap, q::LaurentPolyWrap) = p + (-q) # TODO: optimize

*(p::LaurentPolyWrap, q::LaurentPolyWrap) = LaurentPolyWrap(p.poly * q.poly, p.mindeg + q.mindeg)

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

*(p::LaurentPolyWrap, a::RingElement) = LaurentPolyWrap(p.poly * a, p.mindeg)
*(a::RingElement, p::LaurentPolyWrap) = p * a


###############################################################################
#
#   Powering
#
###############################################################################

function ^(p::LaurentPolyWrap, e::Integer)
   if e >= 0
      LaurentPolyWrap(p.poly^e, p.mindeg * e)
   else
      # p must be a monomial, whose coeff is invertible
      deg = monomial_degree(p)
      c = coeff(p, deg)
      # the following is to allow x^-3 even if 1^-3 is failing
      c = isone(c) ? c : c^e
      LaurentPolyWrap(c * one(p.poly), deg * e)
   end
end

###############################################################################
#
#   LaurentPolynomialRing constructor
#
###############################################################################

@doc doc"""
    LaurentPolynomialRing(R::AbstractAlgebra.Ring, s::AbstractString)
> Given a base ring `R` and string `s` specifying how the generator (variable)
> should be printed, return a tuple `S, x` representing the new Laurent polynomial
> ring $S = R[x, 1/x]$ and the generator $x$ of the ring. The parent
> object `S` will depend only on `R` and `x`.
"""
function LaurentPolynomialRing(R::AbstractAlgebra.Ring, s::AbstractString)
   P, x = AbstractAlgebra.PolynomialRing(R, s)
   LaurentPolyWrapRing(P), LaurentPolyWrap(x)
end
