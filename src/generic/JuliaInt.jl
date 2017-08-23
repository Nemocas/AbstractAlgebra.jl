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

elem_type(::Type{MachineIntegers}) = Int
 
parent_type(::Type{Int}) = MachineIntegers

base_ring(a::Int) = Union{}

base_ring(a::MachineIntegers) = Union{}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

isone(a::Int) = a == 1

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

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function mul!(a::Int, b::Int, c::Int)
   return b*c
end

function addeq!(a::Int, b::Int)
   return a + b
end

function addmul!(a::Int, b::Int, c::Int, d::Int)
   return a + b*c
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
