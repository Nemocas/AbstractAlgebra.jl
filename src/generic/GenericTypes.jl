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

   function PolynomialRing(R::Ring, s::Symbol, cached=true)
      if haskey(PolyID, (R, s))
         return PolyID[R, s]::PolynomialRing{T}
      else 
         z = new{T}(R, s)
         if cached
           PolyID[R, s] = z
         end
         return z
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
      if haskey(ModulusDict, (parent(modulus), modulus))
         return ModulusDict[parent(modulus), modulus]::ResidueRing{T}
      else
         z = new{T}(parent(modulus), modulus)
         if cached
            ModulusDict[parent(modulus), modulus] = z
         end
         return z
      end
   end
end

type Residue{T <: RingElem} <: ResidueElem{T}
   data::T
   parent::ResidueRing{T}

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

   function PowerSeriesRing(R::Ring, prec::Int, s::Symbol, cached=true)
      if haskey(PowerSeriesID, (R, prec, s))
         return PowerSeriesID[R, prec, s]::PowerSeriesRing{T}
      else
         z = new{T}(R, prec, s)
         if cached
            PowerSeriesID[R, prec, s] = z
         end
         return z
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

   function FractionField(R::Ring, cached=true)
      if haskey(FractionDict, R)
         return FractionDict[R]::FractionField{T}
      else
         z = new{T}(R)
         if cached
            FractionDict[R] = z
         end
         return z
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

   function MatrixSpace(R::Ring, r::Int, c::Int, cached=true)
      if haskey(MatrixDict, (R, r, c))
         return MatrixDict[R, r, c]::MatrixSpace{T}
      else
         z = new{T}(r, c, R)
         if cached
            MatrixDict[R, r, c] = z
         end
         return z
      end
   end
end

type Mat{T <: RingElem} <: MatElem{T}
   entries::Array{T, 2}
   parent::MatrixSpace{T}

   Mat(a::Array{T, 2}) = new(a) 
end

