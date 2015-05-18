###############################################################################
#
#   Rings.jl : Generic rings
#
###############################################################################

import Base: length, call, exp, promote_rule, zero, one, show, divrem, mod, 
             hash, factor

export Ring, Field, RingElem

export PolyElem, PowerSeriesElem

###############################################################################
#
#   Abstract types
#
###############################################################################

abstract Ring

abstract Field <: Ring

abstract RingElem

abstract FieldElem <: RingElem

abstract PolyElem <: RingElem

abstract PowerSeriesElem <: RingElem

# not always mathematical ring elements
abstract MatElem <: RingElem

###############################################################################
#
#   Hashing (needed for hashing tuples)
#
###############################################################################

function hash(a::RingElem, b::UInt)
   h = hash(a) $ hash(b)
   h = (h << 1) | (h >> (sizeof(Int)*8 - 1))
   return h
end

###############################################################################
#
#   Generic catchall functions
#
###############################################################################

function +{S <: RingElem, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1
      +(x, parent(x)(y))
   elseif T == T1
      +(parent(y)(x), y)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function +{S <: RingElem, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1
      +(x, parent(x)(y))
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function +{S <: Integer, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if T == T1
      +(parent(y)(x), y)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function -{S <: RingElem, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1
      -(x, parent(x)(y))
   elseif T == T1
      -(parent(y)(x), y)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function -{S <: RingElem, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1
      -(x, parent(x)(y))
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function -{S <: Integer, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if T == T1
      -(parent(y)(x), y)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function *{S <: RingElem, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1
      *(x, parent(x)(y))
   elseif T == T1
      *(parent(y)(x), y)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function *{S <: RingElem, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1
      *(x, parent(x)(y))
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function *{S <: Integer, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if T == T1
      *(parent(y)(x), y)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function divexact{S <: RingElem, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1
      divexact(x, parent(x)(y))
   elseif T == T1
      divexact(parent(y)(x), y)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function divexact{S <: RingElem, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1
      divexact(x, parent(x)(y))
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function divexact{S <: Integer, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if T == T1
      divexact(parent(y)(x), y)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function =={S <: RingElem, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1
      ==(x, parent(x)(y))
   elseif T == T1
      ==(parent(y)(x), y)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function =={S <: RingElem, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1
      ==(x, parent(x)(y))
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function =={S <: Integer, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if T == T1
      ==(parent(y)(x), y)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

###############################################################################
#
#   Exponential function for generic rings
#
###############################################################################

function exp{T <: RingElem}(a::T)
   a != 0 && error("Exponential of nonzero element")
   return one(parent(a))
end

###############################################################################
#
#   Generic and specific rings and fields
#
###############################################################################

include("ZZ.jl")

include("Residue.jl")

include("Poly.jl")

include("fmpz_poly.jl")

include("PowerSeries.jl")

include("fmpz_series.jl")

include("fmpz_mod_series.jl")

include("fmpz_mat.jl")

include("nmod_mat.jl")

include("PariRings.jl")

include("Fields.jl")

include("PariFields.jl")

include("fmpq_poly.jl")

include("padic.jl")

include("fmpq_series.jl")

include("fq_series.jl")

include("fq_nmod_series.jl")

include("pari_poly2.jl")

include("NumberFields.jl")

include("MaximalOrders.jl")

include("PariIdeal.jl")

include("Factor.jl")

include("nmod_poly.jl")

###############################################################################
#
#   Test code
#
###############################################################################

include("../test/Rings-test.jl")

