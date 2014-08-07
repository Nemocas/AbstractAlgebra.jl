export Residue, ResidueRing, modulus, copy, inv

import Base: convert, zero

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

type Residue{T <: Ring, S} <: Ring
   data::T
   Residue(a::Int) = new(mod(convert(T, a), eval(:($S))))
   Residue(a::T) = new(mod(a, eval(:($S))))
   Residue(a::Residue{T, S}) = a
   Residue() = new(T(0))
end

###########################################################################################
#
#   Basic manipulation
#
###########################################################################################

function modulus{T <: Ring, S}(::Type{Residue{T, S}})
   return eval(:($S))
end

zero{T <: Ring, S}(::Type{Residue{T, S}}) = Residue{T, S}(0)

one{T <: Ring, S}(::Type{Residue{T, S}}) = Residue{T, S}(1)

###########################################################################################
#
#   Unary operations
#
###########################################################################################

function -{T <: Ring, S}(a::Residue{T, S})
   Residue{T, S}(-a.data)
end

###########################################################################################
#
#   Comparisons
#
###########################################################################################

=={T, S}(x::Residue{T, S}, y::Residue{T, S}) = x.data == y.data

=={T, S}(x::Residue{T, S}, y::ZZ) = x.data == T(y)

=={T, S}(x::Residue{T, S}, y::Int) = x.data == T(y)

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{T <: Ring, S}(io::IO, x::Residue{T, S})
   print(io, x.data)
end

function show{T <: Ring, S}(io::IO, a::Type{Residue{T, S}})
   print(io, "Residue ring of ", T, " modulo ", modulus(a))
end

###########################################################################################
#
#   Binary operations and functions
#
###########################################################################################

+{T <: Ring, S}(a::Residue{T, S}, b::Residue{T, S}) = Residue{T, S}(a.data + b.data)

-{T <: Ring, S}(a::Residue{T, S}, b::Residue{T, S}) = Residue{T, S}(a.data - b.data)

*{T <: Ring, S}(a::Residue{T, S}, b::Residue{T, S}) = Residue{T, S}(a.data * b.data)

function divexact{T <: Ring, S}(a::Residue{T, S}, b::Residue{T, S})
   g, binv = gcdinv(b.data, eval(:($S)))
   if g != 1
      error("Impossible inverse in divexact")
   end
   Residue{T, S}(a.data * binv)
end

gcd{T <: Ring, S}(a::Residue{T, S}, b::Residue{T, S}) = Residue{T, S}(gcd(a.data, b.data))

###########################################################################################
#
#   Ad hoc binary operations
#
###########################################################################################

*{T <: Ring, S}(a::Residue{T, S}, b::Int) = Residue{T, S}(a.data * b)

*{T <: Ring, S}(a::Int, b::Residue{T, S}) = Residue{T, S}(a * b.data)

*{T <: Ring, S}(a::Residue{T, S}, b::ZZ) = Residue{T, S}(a.data * b)

*{T <: Ring, S}(a::ZZ, b::Residue{T, S}) = Residue{T, S}(a * b.data)

+{T <: Ring, S}(a::Residue{T, S}, b::Int) = Residue{T, S}(a.data + b)

+{T <: Ring, S}(a::Int, b::Residue{T, S}) = Residue{T, S}(a + b.data)

+{T <: Ring, S}(a::Residue{T, S}, b::ZZ) = Residue{T, S}(a.data + b)

+{T <: Ring, S}(a::ZZ, b::Residue{T, S}) = Residue{T, S}(a + b.data)

-{T <: Ring, S}(a::Residue{T, S}, b::Int) = Residue{T, S}(a.data - b)

-{T <: Ring, S}(a::Int, b::Residue{T, S}) = Residue{T, S}(a - b.data)

-{T <: Ring, S}(a::Residue{T, S}, b::ZZ) = Residue{T, S}(a.data - b)

-{T <: Ring, S}(a::ZZ, b::Residue{T, S}) = Residue{T, S}(a - b.data)

###########################################################################################
#
#   Powering
#
###########################################################################################

function ^{T <: Ring, S}(a::Residue{T, S}, b::Int)
   Residue{T, S}(powmod(a.data, b, eval(:($S))))
end

###########################################################################################
#
#   Inversion
#
###########################################################################################

function inv{T <: Ring, S}(a::Residue{T, S})
   g, ainv = gcdinv(a.data, eval(:($S)))
   if g != 1
      error("Impossible inverse in inv")
   end
   Residue{T, S}(ainv)
end

###########################################################################################
#
#   Conversions
#
###########################################################################################

Base.convert{T <: Ring, S}(::Type{Residue{T, S}}, a::T) = Residue{T, S}(a)

Base.convert{T <: Ring, S}(::Type{Residue{T, S}}, a::Int) = Residue{T, S}(a)

###########################################################################################
#
#   ResidueRing constructor
#
###########################################################################################

function ResidueRing{T <: Ring}(::Type{T}, el::T)
   el == 0 && throw(DivideError())
   S = gensym("residue")
   P = Residue{T, S}
   eval(:($S = $el))
   return P
end

function ResidueRing{T <: Ring}(::Type{T}, el::Int)
   el == 0 && throw(DivideError())
   S = gensym("residue")
   P = Residue{T, S}
   eval(:($S = $T($el)))
   return P
end