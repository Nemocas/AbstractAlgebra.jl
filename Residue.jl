export Residue, ResidueRing, modulus

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

=={T, S}(x::Residue{T, S}, y::ZZ) = x.data == convert(T, y)

=={T, S}(x::Residue{T, S}, y::Int) = x.data == convert(T, y)

###########################################################################################
#
#   String I/O
#
###########################################################################################

function show{T <: Ring, S}(io::IO, x::Residue{T, S})
   print(io, "Mod(", x.data, ", ", modulus(Residue{T, S}), ")")
end

function show{T <: Ring, S}(io::IO, a::Type{Residue{T, S}})
   print(io, "Residue ring of ", T, " modulo ", modulus(a))
end

###########################################################################################
#
#   Binary operations
#
###########################################################################################

+{T <: Ring, S}(a::Residue{T, S}, b::Residue{T, S}) = Residue{T, S}(a.data + b.data)

-{T <: Ring, S}(a::Residue{T, S}, b::Residue{T, S}) = Residue{T, S}(a.data - b.data)

*{T <: Ring, S}(a::Residue{T, S}, b::Residue{T, S}) = Residue{T, S}(a.data * b.data)

div{T <: Ring, S}(a::Residue{T, S}, b::Residue{T, S}) = Residue{T, S}(a.data * invmod(b.data, eval(:($S))))

###########################################################################################
#
#   Ad hoc binary operations
#
###########################################################################################

*{T <: Ring, S}(a::Residue{T, S}, b::Int) = Residue{T, S}(a.data * b)

*{T <: Ring, S}(a::Int, b::Residue{T, S}) = Residue{T, S}(a * b.data)

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
   S = gensym("residue")
   P = Residue{T, S}
   eval(:($S = $el))
   return P
end

function ResidueRing{T <: Ring}(::Type{T}, el::Int)
   S = gensym("residue")
   P = Residue{T, S}
   eval(:($S = $T($el)))
   return P
end