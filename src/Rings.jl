###############################################################################
#
#   Rings.jl : Generic rings
#
###############################################################################

import Base: length, call, exp, promote_rule, zero, one, show, divrem, mod, 
             hash, factor

export Ring, Field, RingElem

export PolyElem

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

# not always mathematical ring elements
abstract MatElem <: RingElem

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
#   Generic and specific rings and fields
#
###############################################################################

include("ZZ.jl")

include("Residue.jl")

include("Poly.jl")

include("fmpz_poly.jl")

include("fmpz_mat.jl")

include("PariRings.jl")

# include("Fields.jl")

# include("PariFields.jl")

# include("Poly2.jl")

# include("pari_poly2.jl")

# include("NumberFields.jl")

# include("MaximalOrders.jl")

# include("PariIdeal.jl")

# include("Factor.jl")

###############################################################################
#
#   Test code
#
###############################################################################

include("../test/Rings-test.jl")

