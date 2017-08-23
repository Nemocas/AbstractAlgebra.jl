###############################################################################
#
#   JuliaBigInt.jl : Additional Nemo functionality for Julia BigInts
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

JuliaZZ = Integers()

parent(a::BigInt) = JuliaZZ

elem_type(::Integers) = BigInt
 
parent_type(::Type{BigInt}) = Integers

base_ring(a::BigInt) = Union{}

base_ring(a::Integers) = Union{}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(::Integers) = BigInt(0)

one(::Integers) = BigInt(1)

isone(a::BigInt) = a == 1

canonical_unit(a::BigInt) = a < 0 ? BigInt(-1) : BigInt(1)

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::Integers)
   print(io, "Integers")
end

needs_parentheses(::BigInt) = false

isnegative(a::BigInt) = a < 0

show_minus_one(::Type{BigInt}) = false

###############################################################################
#
#   Divides
#
###############################################################################

function divides(a::BigInt, b::BigInt)
   q, r = divrem(a, b)
   return r == 0, q
end

###############################################################################
#
#   Exact division
#
###############################################################################

divexact(a::BigInt, b::Int) = div(a, b)

divexact(a::BigInt, b::BigInt) = div(a, b)

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::BigInt)
   ccall((:__gmpz_set_si, :libgmp), Void, (Ptr{BigInt}, Int), &a, 0)
   return a
end

function mul!(a::BigInt, b::BigInt, c::BigInt)
   ccall((:__gmpz_mul, :libgmp), Void, (Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}), &a, &b, &c)
   return a
end

function addeq!(a::BigInt, b::BigInt)
   ccall((:__gmpz_add, :libgmp), Void, (Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}), &a, &a, &b)
   return a
end

function add!(a::BigInt, b::BigInt, c::BigInt)
   ccall((:__gmpz_add, :libgmp), Void, (Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}), &a, &b, &c)
   return a
end

function addmul!(a::BigInt, b::BigInt, c::BigInt, d::BigInt)
   ccall((:__gmpz_addmul, :libgmp), Void, (Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}), &a, &b, &c)
   return a
end

function addmul!(a::BigInt, b::BigInt, c::BigInt) # special case, no temporary required
   ccall((:__gmpz_addmul, :libgmp), Void, (Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}), &a, &b, &c)
   return a
end

###############################################################################
#
#   Random generation
#
###############################################################################

function rand(R::Integers, n::UnitRange{Int})
   return R(rand(n))
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::Integers)()
   return BigInt(0)
end

function (a::Integers)(b::Int)
   return BigInt(b)
end

function (a::Integers)(b::BigInt)
   return b
end
