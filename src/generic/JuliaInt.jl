###############################################################################
#
#   JuliaInt.jl : Additional Nemo functionality for Julia Ints
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

zz = MachineIntegers()

parent(a::Int) = zz

elem_type(::MachineIntegers) = Int
 
parent_type(::Type{Int}) = MachineIntegers

base_ring(a::Int) = Union{}

base_ring(a::MachineIntegers) = Union{}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(::MachineIntegers) = 0

one(::MachineIntegers) = 1

isone(a::Int) = a == 1

canonical_unit(a::Int) = a < 0 ? -1 : 1

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::MachineIntegers)
   print(io, "Machine Integers")
end

needs_parentheses(::Int) = false

isnegative(a::Int) = a < 0

show_minus_one(::Type{Int}) = false

###############################################################################
#
#   Divides
#
###############################################################################

function divides(a::Int, b::Int)
   q, r = divrem(a, b)
   return r == 0, q
end

###############################################################################
#
#   Exact division
#
###############################################################################

divexact(a::Int, b::Int) = div(a, b)

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::Int)
   return 0
end

function mul!(a::Int, b::Int, c::Int)
   return b*c
end

function add!(a::Int, b::Int, c::Int)
   return b + c
end

function addeq!(a::Int, b::Int)
   return a + b
end

function addmul!(a::Int, b::Int, c::Int, d::Int)
   return a + b*c
end

function addmul!(a::Int, b::Int, c::Int) # special case, no temporary required
   return a + b*c
end

###############################################################################
#
#   Random generation
#
###############################################################################

function rand(R::MachineIntegers, n::UnitRange{Int})
   return R(rand(n))
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::MachineIntegers)()
   return 0
end

function (a::MachineIntegers)(b::Int)
   return b
end
