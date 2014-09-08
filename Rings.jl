module Rings

import Base: exp

export Ring, Field, exp

abstract Ring

abstract Field <: Ring

function +{S <: Ring, T <: Ring}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      +(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function +{S <: Ring, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      +(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function +{S <: Integer, T <: Ring}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      +(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function -{S <: Ring, T <: Ring}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      -(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function -{S <: Ring, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      -(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function -{S <: Integer, T <: Ring}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      -(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function *{S <: Ring, T <: Ring}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      *(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function *{S <: Ring, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      *(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function *{S <: Integer, T <: Ring}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      *(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function =={S <: Ring, T <: Ring}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      ==(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function =={S <: Ring, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      ==(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function =={S <: Integer, T <: Ring}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      ==(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

isequal{T <: Ring}(a::T, b::T) = a == b

function divexact{S <: Ring, T <: Ring}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      divexact(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function divexact{S <: Ring, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      divexact(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function divexact{S <: Integer, T <: Ring}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      divexact(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function exp{T <: Ring}(a::T)
   a != 0 && error("Exponential of nonzero element")
   return one(T)
end

include("ZZ.jl")

include("Residue.jl")

include("Poly.jl")

include("PowerSeries.jl")

include("Rings-test.jl")

end # module
