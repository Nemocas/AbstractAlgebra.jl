module Rings

export Ring, Field

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

include("ZZ.jl")

include("Residue.jl")

include("Poly.jl")

end # module
