###########################################################################################
#
#   Rings.jl : Generic rings
#
###########################################################################################

import Base: length, call, exp, promote_rule, zero, one, show, divrem, mod, hash, factor

export Ring, Field, RingElem

export PolyElem

abstract Ring

abstract Field <: Ring

abstract RingElem

abstract FieldElem <: RingElem

abstract PolyElem <: RingElem

# not always mathematical ring elements
abstract MatElem <: RingElem

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

include("ZZ.jl")

include("Residue.jl")

include("Poly.jl")

include("fmpz_poly.jl")

include("fmpz_mat.jl")

include("PariRings.jl")

include("Fields.jl")

include("PariFields.jl")

include("Poly2.jl")

include("pari_poly2.jl")

include("NumberFields.jl")

include("MaximalOrders.jl")

include("PariIdeal.jl")

###########################################################################################
#
#   Factorisation
#
###########################################################################################

type Factor{T <: Ring}
   d::PariFactor
   len::Int
   parent::T
end

function getindex(a::Factor{IntegerRing}, i::Int)
   p, n = a.d[i]
   return ZZ(p), n
end

function getindex{S}(a::Factor{FmpzPolyRing{S}}, i::Int)
   p, n = a.d[i]
   return a.parent(p), n
end

function getindex{S}(a::Factor{FmpqPolyRing{S}}, i::Int)
   p, n = a.d[i]
   return a.parent(p), n
end

function show(io::IO, a::Factor)
   print(io, "[")
   for i = 1:a.len
      print(io, a[i])
      if i != a.len
         print(io, ", ")
      end
   end
   print(io, "]")
end

function factor(n::BigInt)
   f = factor(pari(n))
   return Factor{IntegerRing}(f, f.len, ZZ)
end

function factor{S}(g::fmpz_poly{S})
   f = factor(pari(g))
   return Factor{FmpzPolyRing{S}}(f, f.len, g.parent)
end
function factor{S}(g::fmpq_poly{S})
   f = factor(pari(g))
   return Factor{FmpqPolyRing{S}}(f, f.len, g.parent)
end

include("../test/Rings-test.jl")

