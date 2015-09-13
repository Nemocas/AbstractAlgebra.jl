###############################################################################
#
#   Rings.jl : Generic rings
#
###############################################################################

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

include("flint/fmpz.jl")

include("generic/Residue.jl")

include("generic/Poly.jl")

include("flint/fmpz_poly.jl")

include("flint/nmod_poly.jl")

include("flint/fmpz_mod_poly.jl")

include("generic/PowerSeries.jl")

include("flint/fmpz_series.jl")

include("flint/fmpz_mod_series.jl")

include("generic/Matrix.jl")

include("flint/fmpz_mat.jl")

include("flint/nmod_mat.jl")

include("pari/pari_int.jl")

include("pari/pari_poly.jl")

include("pari/pari_polmod.jl")

include("pari/pari_vec.jl")

include("pari/PariFactor.jl")

include("Fields.jl")

include("pari/pari_frac.jl")

include("flint/fmpq_poly.jl")

include("flint/padic.jl")

include("flint/fmpq_series.jl")

include("flint/fq_series.jl")

include("flint/fq_nmod_series.jl")

include("flint/fq_poly.jl")

include("flint/fq_nmod_poly.jl")

include("pari/pari_poly2.jl")

include("pari/pari_maximal_order_elem.jl")

include("pari/PariIdeal.jl")

include("Factor.jl")

