import Base: length, call, exp, promote_rule, zero, one, show, divrem, mod

export Ring, Field, RingElem

export PolyElem

abstract Ring

abstract Field <: Ring

abstract RingElem

abstract FieldElem <: RingElem

abstract PolyElem <: RingElem

abstract Residue <: RingElem

function +{S <: RingElem, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      +(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function +{S <: RingElem, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      +(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function +{S <: Integer, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      +(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function -{S <: RingElem, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      -(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function -{S <: RingElem, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      -(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function -{S <: Integer, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      -(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function *{S <: RingElem, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      *(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function *{S <: RingElem, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      *(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function *{S <: Integer, T <: RingElem}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      *(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

include("ZZ.jl")

include("Poly.jl")

include("fmpz_poly.jl")
