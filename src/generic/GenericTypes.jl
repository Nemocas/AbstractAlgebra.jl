###############################################################################
#
#   GenericTypes.jl : Generic Types
#
###############################################################################

###############################################################################
#
#   PolynomialRing / Poly
#
###############################################################################

PolyID = ObjectIdDict()

type PolynomialRing{T <: RingElem} <: Ring
   base_ring :: Ring
   S::Symbol

   function PolynomialRing(R::Ring, s::Symbol)
      return try
         PolyID[R, s]
      catch
         PolyID[R, s] = new(R, s)
      end
   end
end

type Poly{T <: RingElem} <: PolyElem
   coeffs::Array{T, 1}
   length::Int
   parent::PolynomialRing{T}

   Poly() = new(Array(T, 0), 0)
   
   Poly(a::Array{T, 1}) = new(a, length(a))

   Poly(a::T) = a == 0 ? new(Array(T, 0), 0) : new([a], 1)
end

###############################################################################
#
#   ResidueRing / Residue
#
###############################################################################

ModulusDict = Dict{Tuple{Ring, RingElem}, Ring}()

type ResidueRing{T <: RingElem} <: Ring
   base_ring::Ring
   modulus::T

   function ResidueRing(modulus::T)
      return try
         ModulusDict[parent(modulus), modulus]
      catch
         ModulusDict[parent(modulus), modulus] = new(parent(modulus), modulus)
      end
   end
end

type Residue{T <: RingElem} <: RingElem
   data::T
   parent::ResidueRing

   Residue(a::T) = new(a)
end

###############################################################################
#
#   PowerSeriesRing / PowerSeries
#
###############################################################################

PowerSeriesID = ObjectIdDict()

type PowerSeriesRing{T <: RingElem} <: Ring
   base_ring::Ring
   prec_max::Int
   S::Symbol

   function PowerSeriesRing(R::Ring, prec::Int, s::Symbol)
      return try
         PowerSeriesID[R, prec, s]
      catch
         PowerSeriesID[R, prec, s] = new(R, prec, s)
      end
   end
end

type PowerSeries{T <: RingElem} <: PowerSeriesElem
   coeffs::Array{T, 1}
   length::Int
   prec::Int
   parent::PowerSeriesRing{T}

   PowerSeries(a::Array{T, 1}, length::Int, prec::Int) = new(a, length, prec)   
   PowerSeries(a::PowerSeries{T}) = a
end

###############################################################################
#
#   FractionField / Fraction
#
###############################################################################

FractionDict = ObjectIdDict()

type FractionField{T <:RingElem} <: Field
   base_ring::Ring

   function FractionField(R::Ring)
      return try
         FractionDict[R]
      catch
         FractionDict[R] = new(R)
      end
   end
end

type Fraction{T <: RingElem} <: FieldElem
   num::T
   den::T
   parent::FractionField{T}

   Fraction(num::T, den::T) = new(num, den) 
end
