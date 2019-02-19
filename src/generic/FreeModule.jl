###############################################################################
#
#   FreeModule.jl : Free modules over rings
#
###############################################################################

export FreeModule, free_module_elem

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{free_module_elem{T}}) where T <: RingElement = FreeModule{T}

base_ring(M::FreeModule{T}) where T <: RingElement = M.base_ring::parent_type(T)

base_ring(v::free_module_elem{T}) where T <: RingElement = base_ring(parent(v))

elem_type(::Type{FreeModule{T}}) where T <: RingElement = free_module_elem{T}

parent(m::free_module_elem{T}) where T <: RingElement = m.parent

function check_parent(m1::free_module_elem{T}, m2::free_module_elem{T}) where T <: RingElement
    parent(m1) != parent(m2) && ("Incompatible free modules")
end

function isdomain_type(::Type{free_module_elem{T}}) where T <: RingElement
   return isdomain_type(T)
end

function isexact_type(a::Type{free_module_elem{T}}) where T <: RingElement
   return isexact_type(T)
end

rank(M::FreeModule{T}) where T <: RingElement = M.rank

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, M::FreeModule)
   print(io, "Free module of rank ")
   print(io, M.rank)
   print(" over ")
   show(io, base_ring(M))
end

function show(io::IO, a::free_module_elem)
   print(io, "(")
   print(io, join(string.(a.v), ", "))
   print(io, ")")
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(m::free_module_elem)
   M = parent(m)
   return M(-m.v)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(m1::free_module_elem{T}, m2::free_module_elem{T}) where T <: RingElement
   check_parent(m1, m2)
   M = parent(m1)
   return M(m1.v + m2.v)
end

function -(m1::free_module_elem{T}, m2::free_module_elem{T}) where T <: RingElement
   check_parent(m1, m2)
   M = parent(m1)
   return M(m1.v - m2.v)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(m::free_module_elem{T}, c::T) where T <: RingElement
   parent(c) != base_ring(m) && error("Incompatible scalar")
   M = parent(m)
   return M([r*c for r in m.v])
end

function *(c::T, m::free_module_elem{T}) where T <: RingElement
   parent(c) != base_ring(m) && error("Incompatible scalar")
   M = parent(m)
   return M([c*r for r in m.v])
end

function *(m::free_module_elem{T}, c::U) where {T <: RingElement, U <: Integer}
   M = parent(m)
   return M([r*c for r in m.v])
end

function *(c::U, m::free_module_elem{T}) where {T <: RingElement, U <: Integer}
   M = parent(m)
   return M([c*r for r in m.v])
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(m1::free_module_elem{T}, m2::free_module_elem{T}) where T <: RingElement
   check_parent(m1, m2)
   return m1.v == m2.v
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (M::FreeModule{T})(a::Vector{T}) where T <: RingElement
   length(a) != rank(M) && error("Number of elements does not equal rank")
   z = free_module_elem{T}(a)
   z.parent = M
   return z
end

###############################################################################
#
#   FreeModule constructor
#
###############################################################################

function FreeModule(R::Ring, rank::Int; cached::Bool = true)
   T = elem_type(R)
   return FreeModule{T}(R, rank, cached)
end

