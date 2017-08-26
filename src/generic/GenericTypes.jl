###############################################################################
#
#   GenericTypes.jl : Generic Types
#
###############################################################################

###############################################################################
#
#   Integers / BigInt
#
###############################################################################

mutable struct Integers
end

###############################################################################
#
#   MachineIntegers / Int
#
###############################################################################

mutable struct MachineIntegers
end

###############################################################################
#
#   Rationals / Rational
#
###############################################################################

mutable struct Rationals{T <: Union{Int, BigInt}}
end

###############################################################################
#
#   Unions of Nemo abstract types and Julia types
#
###############################################################################

const RingParent = Union{Ring, Integers, MachineIntegers, Rationals{Int}, Rationals{BigInt}}

const RingElement = Union{RingElem, Int, BigInt, Rational{Int}, Rational{BigInt}}

const FieldParent = Union{Field, Rationals{Int}, Rationals{BigInt}}

const FieldElement = Union{FieldElem, Rational{Int}, Rational{BigInt}}

###############################################################################
#
#   GenPolyRing / GenPoly
#
###############################################################################

mutable struct GenPolyRing{T <: RingElement} <: PolyRing{T}
   base_ring::RingParent
   S::Symbol

   function GenPolyRing{T}(R::RingParent, s::Symbol, cached::Bool = true) where T <: RingElement
      if haskey(GenPolyID, (R, s))
         return GenPolyID[R, s]::GenPolyRing{T}
      else 
         z = new{T}(R, s)
         if cached
           GenPolyID[R, s] = z
         end
         return z
      end
   end
end

const GenPolyID = Dict{Tuple{RingParent, Symbol}, Ring}()

mutable struct GenPoly{T <: RingElement} <: PolyElem{T}
   coeffs::Array{T, 1}
   length::Int
   parent::GenPolyRing{T}

   GenPoly{T}() where T <: RingElement = new{T}(Array{T}(0), 0)
   
   GenPoly{T}(a::Array{T, 1}) where T <: RingElement = new{T}(a, length(a))

   GenPoly{T}(a::T) where T <: RingElement = iszero(a) ? new{T}(Array{T}(0), 0) : new{T}([a], 1)
end

###############################################################################
#
#   GenMPolyRing / GenMPoly / Monomial
#
###############################################################################

# S is a Symbol which can take the values:
# :lex
# :revlex
# :deglex
# :degrevlex
# 
# T is an Int which is the number of variables
# (plus one if ordered by total degree)

mutable struct GenMPolyRing{T <: RingElement} <: PolyRing{T}
   base_ring::RingParent
   S::Array{Symbol, 1}
   ord::Symbol
   num_vars::Int
   N::Int

   function GenMPolyRing{T}(R::RingParent, s::Array{Symbol, 1}, ord::Symbol, N::Int,
                         cached::Bool = true) where T <: RingElement
      if haskey(GenMPolyID, (R, s, ord, N))
         return GenMPolyID[R, s, ord, N]::GenMPolyRing{T}
      else 
         z = new{T}(R, s, ord, length(s), N)
         if cached
           GenMPolyID[R, s, ord, N] = z
         end
         return z
      end
   end
end

const GenMPolyID = Dict{Tuple{RingParent, Array{Symbol, 1}, Symbol, Int}, RingParent}()

mutable struct GenMPoly{T <: RingElement} <: PolyElem{T}
   coeffs::Array{T, 1}
   exps::Array{UInt, 2}
   length::Int
   parent::GenMPolyRing{T}

   function GenMPoly{T}(R::GenMPolyRing) where T <: RingElement
      N = R.N
      return new{T}(Array{T}(0), Array{UInt}(N, 0), 0, R)
   end
   
   GenMPoly{T}(R::GenMPolyRing, a::Array{T, 1}, b::Array{UInt, 2}) where T <: RingElement = new{T}(a, b, length(a), R)

   function GenMPoly{T}(R::GenMPolyRing, a::T) where T <: RingElement
      N = R.N
      return iszero(a) ? new{T}(Array{T}(0), Array{UInt}(N, 0), 0, R) : 
                                          new{T}([a], zeros(UInt, N, 1), 1, R)
   end
end

###############################################################################
#
#   GenSparsePolyRing / GenSparsePoly
#
###############################################################################

mutable struct GenSparsePolyRing{T <: RingElement} <: Ring
   base_ring::RingParent
   S::Symbol
   num_vars::Int

   function GenSparsePolyRing{T}(R::RingParent, s::Symbol, cached::Bool = true) where T <: RingElement
      if haskey(GenSparsePolyID, (R, s))
         return GenSparsePolyID[R, s]::GenSparsePolyRing{T}
      else 
         z = new{T}(R, s)
         if cached
           GenSparsePolyID[R, s] = z
         end
         return z
      end
   end
end

const GenSparsePolyID = Dict{Tuple{RingParent, Symbol}, GenSparsePolyRing}()

mutable struct GenSparsePoly{T <: RingElement} <: RingElem
   coeffs::Array{T, 1}
   exps::Array{UInt}
   length::Int
   parent::GenSparsePolyRing{T}

   GenSparsePoly{T}() where T <: RingElement = new{T}(Array{T}(0), Array{UInt}(0), 0)
   
   GenSparsePoly{T}(a::Array{T, 1}, b::Array{UInt, 1}) where T <: RingElement = new{T}(a, b, length(a))

   GenSparsePoly{T}(a::T) where T <: RingElement = iszero(a) ? new{T}(Array{T}(0), Array{UInt}(0), 0) : 
                                               new{T}([a], [UInt(0)], 1)
end

###############################################################################
#
#   GenResRing / GenRes
#
###############################################################################

mutable struct GenResRing{T <: RingElement} <: ResRing{T}
   base_ring::RingParent
   modulus::T

   function GenResRing{T}(modulus::T, cached::Bool = true) where T <: RingElement
      if haskey(ModulusDict, (parent(modulus), modulus))
         return ModulusDict[parent(modulus), modulus]::GenResRing{T}
      else
         z = new{T}(parent(modulus), modulus)
         if cached
            ModulusDict[parent(modulus), modulus] = z
         end
         return z
      end
   end
end

const ModulusDict = Dict{Tuple{RingParent, RingElement}, RingParent}()

mutable struct GenRes{T <: RingElement} <: ResElem{T}
   data::T
   parent::GenResRing{T}

   GenRes{T}(a::T) where T <: RingElement = new{T}(a)
end

###############################################################################
#
#   GenRelSeriesRing / GenRelSeries
#
###############################################################################

mutable struct GenRelSeriesRing{T <: RingElement} <: SeriesRing{T}
   base_ring::RingParent
   prec_max::Int
   S::Symbol

   function GenRelSeriesRing{T}(R::RingParent, prec::Int, s::Symbol, cached::Bool = true) where T <: RingElement
      if haskey(GenRelSeriesID, (R, prec, s))
         return GenRelSeriesID[R, prec, s]::GenRelSeriesRing{T}
      else
         z = new{T}(R, prec, s)
         if cached
            GenRelSeriesID[R, prec, s] = z
         end
         return z
      end
   end
end

const GenRelSeriesID = Dict{Tuple{RingParent, Int, Symbol}, RingParent}()

mutable struct GenRelSeries{T <: RingElement} <: RelSeriesElem{T}
   coeffs::Array{T, 1}
   length::Int
   prec::Int
   val::Int
   parent::GenRelSeriesRing{T}

   function GenRelSeries{T}(a::Array{T, 1}, length::Int, prec::Int, val::Int) where T <: RingElement
      new{T}(a, length, prec, val)
   end

   GenRelSeries{T}(a::GenRelSeries{T}) where T <: RingElement = a
end

###############################################################################
#
#   GenAbsSeriesRing / GenAbsSeries
#
###############################################################################

mutable struct GenAbsSeriesRing{T <: RingElement} <: SeriesRing{T}
   base_ring::RingParent
   prec_max::Int
   S::Symbol

   function GenAbsSeriesRing{T}(R::RingParent, prec::Int, s::Symbol, cached::Bool = true) where T <: RingElement
      if haskey(GenAbsSeriesID, (R, prec, s))
         return GenAbsSeriesID[R, prec, s]::GenAbsSeriesRing{T}
      else
         z = new{T}(R, prec, s)
         if cached
            GenAbsSeriesID[R, prec, s] = z
         end
         return z
      end
   end
end

const GenAbsSeriesID = Dict{Tuple{RingParent, Int, Symbol}, RingParent}()

mutable struct GenAbsSeries{T <: RingElement} <: AbsSeriesElem{T}
   coeffs::Array{T, 1}
   length::Int
   prec::Int
   parent::GenAbsSeriesRing{T}

   GenAbsSeries{T}(a::Array{T, 1}, length::Int, prec::Int) where T <: RingElement = new{T}(a, length, prec)   
   GenAbsSeries{T}(a::GenAbsSeries{T}) where T <: RingElement = a
end

###############################################################################
#
#   GenFracField / GenFrac
#
###############################################################################

mutable struct GenFracField{T <: RingElem} <: FracField{T}
   base_ring::Ring

   function GenFracField{T}(R::Ring, cached::Bool = true) where {T}
      if haskey(GenFracDict, R)
         return GenFracDict[R]::GenFracField{T}
      else
         z = new{T}(R)
         if cached
            GenFracDict[R] = z
         end
         return z
      end
   end
end

const GenFracDict = Dict{Ring, Ring}()

mutable struct GenFrac{T <: RingElem} <: FracElem{T}
   num::T
   den::T
   parent::GenFracField{T}

   GenFrac{T}(num::T, den::T) where {T} = new{T}(num, den) 
end

###############################################################################
#
#   GenMatSpace / GenMat
#
###############################################################################

# not really a mathematical ring
mutable struct GenMatSpace{T <: RingElement} <: MatSpace{T}
   rows::Int
   cols::Int
   base_ring::RingParent

   function GenMatSpace{T}(R::RingParent, r::Int, c::Int, cached::Bool = true) where T <: RingElement
      if haskey(GenMatDict, (R, r, c))
         return GenMatDict[R, r, c]::GenMatSpace{T}
      else
         z = new{T}(r, c, R)
         if cached
            GenMatDict[R, r, c] = z
         end
         return z
      end
   end
end

const GenMatDict = Dict{Tuple{RingParent, Int, Int}, RingParent}()

mutable struct GenMat{T <: RingElement} <: MatElem{T}
   entries::Array{T, 2}
   base_ring::RingParent

   GenMat{T}(a::Array{T, 2}) where T <: RingElement = new{T}(a) 
end
