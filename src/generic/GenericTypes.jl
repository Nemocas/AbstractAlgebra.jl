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
      if haskey(PolyID, (R, s))
         return PolyID[R, s]::PolynomialRing{T}
      else 
         PolyID[R, s] = new{T}(R, s)
         return PolyID[R, s]::PolynomialRing{T}
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

   function ResidueRing(modulus::T)
      if haskey(ModulusDict, (parent(modulus), modulus))
         return ModulusDict[parent(modulus), modulus]::ResidueRing{T}
      else
         ModulusDict[parent(modulus), modulus] = new{T}(parent(modulus), modulus)
         return ModulusDict[parent(modulus), modulus]
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
      if haskey(PowerSeriesID, (R, prec, s))
         return PowerSeriesID[R, prec, s]::PowerSeriesRing{T}
      else
         PowerSeriesID[R, prec, s] = new{T}(R, prec, s)
         return PowerSeriesID[R, prec, s]
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
      if haskey(FractionDict, R)
         return FractionDict[R]::FractionField{T}
      else
         FractionDict[R] = new{T}(R)
         return FractionDict[R]
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
      if haskey(MatrixDict, (R, r, c))
         return MatrixDict[R, r, c]::MatrixSpace{T}
      else
         MatrixDict[R, r, c] = new{T}(r, c, R)::MatrixSpace{T}
      end
   end
end

type Mat{T <: RingElem} <: MatElem{T}
   entries::Array{T, 2}
   parent::MatrixSpace{T}

   Mat(a::Array{T, 2}) = new(a) 
end
