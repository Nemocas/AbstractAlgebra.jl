###############################################################################
#
#   JuliaRational.jl : Additional Nemo functionality for Julia Rationals
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

JuliaQQ = Rationals{BigInt}()

qq = Rationals{Int}()

parent(a::Rational{T}) where T <: Union{Int, BigInt} = Rationals{T}()

elem_type(::Type{Rationals{T}}) where T <: Union{Int, BigInt} = Rational{T}
  
parent_type(::Type{Rational{T}}) where T <: Union{Int, BigInt} = Rationals{T}

base_ring(a::Rational{Int}) = zz

base_ring(a::Rational{BigInt}) = JuliaZZ

base_ring(a::Rationals{Int}) = zz

base_ring(a::Rationals{BigInt}) = JuliaZZ

###############################################################################
#
#   Basic manipulation
#
###############################################################################

isone(a::Rational{T}) where T <: Union{Int, BigInt} = a == 1

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::Rationals{T}) where T <: Union{Int, BigInt}
   print(io, "Rationals over ")
   show(io, base_ring(R))
end

needs_parentheses(::Rational{T}) where T <: Union{Int, BigInt} = false

isnegative(a::Rational{T}) where T <: Union{Int, BigInt} = a < 0

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::Rational{T}) where T <: Union{Int, BigInt}
   return Rational{T}(0)
end

function mul!(a::Rational{T}, b::Rational{T}, c::Rational{T}) where T <: Union{Int, BigInt}
   return b*c
end

function add!(a::Rational{T}, b::Rational{T}, c::Rational{T}) where T <: Union{Int, BigInt}
   return b + c
end

function addeq!(a::Rational{T}, b::Rational{T}) where T <: Union{Int, BigInt}
   return a + b
end

function addmul!(a::Rational{T}, b::Rational{T}, c::Rational{T}, d::Rational{T}) where T <: Union{Int, BigInt}
   return a + b*c
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (R::Rationals{Int})()
   return Rational{Int}(0)
end

function (R::Rationals{BigInt})()
   return Rational{BigInt}(0)
end

function (R::Rationals{Int})(b::Int)
   return Rational{Int}(b)
end

function (R::Rationals{BigInt})(b::Int)
   return Rational{BigInt}(b)
end

function (R::Rationals{BigInt})(b::BigInt)
   return Rational{BigInt}(b)
end

###############################################################################
#
#   FractionField constructor
#
###############################################################################

FractionField(R::MachineIntegers) = Rationals{Int}()

FractionField(R::Integers) = Rationals{BigInt}()