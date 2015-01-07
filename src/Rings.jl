import Base: length, call, exp

export Ring, Field, RingElem, exp

export Poly, PolyRing, PolynomialRing, coeff, isgen, truncate, mullow, reverse,
       shift_left, shift_right, divexact, pseudorem, pseudodivrem, gcd, content, primpart,
       evaluate, compose, derivative, resultant, discriminant, bezout

abstract Ring

abstract Field <: Ring

abstract RingElem

abstract Poly <: RingElem

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

PolyID = ObjectIdDict()

type PolyRing{P <: Poly, S} <: Ring
   base_ring :: Ring

   function PolyRing(R::Ring)
      try
         PolyID[R, S]
      catch
         PolyID[R, S] = new(R)
      end
   end
end

include("ZZ.jl")

include("fmpz_poly.jl")
