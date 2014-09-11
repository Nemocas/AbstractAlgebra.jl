module Fields

using Rings

function /{S <: Ring, T <: Ring}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      /(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function /{S <: Ring, T <: Integer}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      /(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

function /{S <: Integer, T <: Ring}(x::S, y::T) 
   T1 = promote_type(S, T)
   if S == T1 || T == T1
      /(promote(x, y)...)
   else
      error("Unable to promote ", S, " and ", T, " to common type")
   end
end

include("Fraction.jl")

include("FiniteFields.jl")

include("Padics2.jl")

include("Poly2.jl")

include("PowerSeries2.jl")

include("../test/Fields-test.jl")

end # module

