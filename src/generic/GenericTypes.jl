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

const PolyID = ObjectIdDict()

type PolynomialRing{T <: RingElem} <: Ring{Generic}
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

type Poly{T <: RingElem} <: PolyElem{T}
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

const ModulusDict = Dict{Tuple{Ring, RingElem}, Ring}()

type ResidueRing{T <: RingElem} <: Ring{Generic}
   base_ring::Ring
   modulus::T

   function ResidueRing(modulus::T, cached=true)
      return try
         ModulusDict[parent(modulus), modulus]
      catch
         R = new(parent(modulus), modulus)
         cached ? ModulusDict[parent(modulus), modulus] = R : R
      end
   end
end

type Residue{T <: RingElem} <: ResidueElem{T}
   data::T
   parent::ResidueRing

   Residue(a::T) = new(a)
end

###############################################################################
#
#   PowerSeriesRing / PowerSeries
#
###############################################################################

const PowerSeriesID = ObjectIdDict()

type PowerSeriesRing{T <: RingElem} <: Ring{Generic}
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

type PowerSeries{T <: RingElem} <: SeriesElem{T}
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

const FractionDict = ObjectIdDict()

type FractionField{T <: RingElem} <: Field{Generic}
   base_ring::Ring

   function FractionField(R::Ring)
      return try
         FractionDict[R]
      catch
         FractionDict[R] = new(R)
      end
   end
end

type Fraction{T <: RingElem} <: FractionElem{T}
   num::T
   den::T
   parent::FractionField{T}

   Fraction(num::T, den::T) = new(num, den) 
end

###############################################################################
#
#   MatrixSpace / Matrix
#
###############################################################################

const MatrixDict = ObjectIdDict()

# not really a mathematical ring
type MatrixSpace{T <: RingElem} <: Ring{Generic}
   rows::Int
   cols::Int
   base_ring::Ring

   function MatrixSpace(R::Ring, r::Int, c::Int)
      return try
         MatrixDict[R, r, c]
      catch
         MatrixDict[R, r, c] = new(r, c, R)
      end
   end
end

type Mat{T <: RingElem} <: MatElem{T}
   entries::Array{T, 2}
   parent::MatrixSpace{T}

   Mat(a::Array{T, 2}) = new(a) 
end
