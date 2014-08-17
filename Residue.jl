export Residue, ResidueRing, modulus, copy, inv, canonical_unit

import Base: convert, zero

###########################################################################################
#
#   Data types and memory management
#
###########################################################################################

type Residue{T <: Ring, S} <: Ring
   data::T
   Residue(a::Int) = new(mod(convert(T, a), eval(:($S))))
   Residue(a::ZZ) = new(mod(convert(T, a), eval(:($S))))
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

=={T, S}(x::ZZ, y::Residue{T, S}) = T(x) == y.data

=={T, S}(x::Int, y::Residue{T, S}) = T(x) == y.data

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

needs_parentheses{T <: Ring, S}(x::Residue{T, S}) = needs_parentheses(x.data)

is_negative{T <: Ring, S}(x::Residue{T, S}) = is_negative(x.data)

show_minus_one{T <: Ring, S}(::Type{Residue{T, S}}) = true

###########################################################################################
#
#   Canonicalisation
#
###########################################################################################

canonical_unit{T, S}(a::Residue{T, S}) = a

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

gcd{T <: Ring, S}(a::Residue{T, S}, b::Residue{T, S}) = Residue{T, S}(gcd(gcd(a.data, eval(:($S))), b.data))

###########################################################################################
#
#   Unsafe operators and functions
#
###########################################################################################

function mul!{T <: Ring, S}(c::Residue{T, S}, a::Residue{T, S}, b::Residue{T, S})
   c.data = mod(a.data*b.data, eval(:($S)))
end

function addeq!{T <: Ring, S}(c::Residue{T, S}, a::Residue{T, S})
   c.data = mod(c.data + a.data, eval(:($S)))
end

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