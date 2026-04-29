###############################################################################
#
#   FreeModule.jl : Free modules over rings
#
###############################################################################

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{FreeModuleElem{T}}) where T <: NCRingElement = FreeModule{T}

base_ring_type(::Type{FreeModule{T}}) where T <: NCRingElement = parent_type(T)

base_ring(M::FreeModule{T}) where T <: NCRingElement = M.base_ring::parent_type(T)

coefficient_ring_type(T::Type{<:FreeModule}) = base_ring_type(T)

coefficient_ring(M::FreeModule) = base_ring(M)

elem_type(::Type{FreeModule{T}}) where T <: NCRingElement = FreeModuleElem{T}

parent(m::FreeModuleElem{T}) where T <: NCRingElement = m.parent

function rels(M::FreeModule{T}) where T <: NCRingElement
   # there are no relations in a free module
   return Vector{dense_matrix_type(T)}(undef, 0)
end

is_free(M::FreeModule) = true

@doc raw"""
    rank(M::FreeModule{T}) where T <: NCRingElement

Return the rank of the given free module.
"""
rank(M::FreeModule{T}) where T <: NCRingElement = M.rank

@doc raw"""
    dim(M::FreeModule{T}) where T <: FieldElement

Return the dimension of the given vector space.
"""
dim(M::FreeModule{T}) where T <: FieldElement = M.rank
vector_space_dim(M::FreeModule{T}) where T <: FieldElement = M.rank

number_of_generators(M::FreeModule{T}) where T <: NCRingElement = M.rank

function gens(N::FreeModule{T}) where T <: NCRingElement
   return [gen(N, i) for i = 1:ngens(N)]
end

function gen(N::FreeModule{T}, i::Int) where T <: NCRingElement
   @boundscheck 1 <= i <= ngens(N) || throw(ArgumentError("generator index out of range"))
   R = base_ring(N)
   if N.is_row
     m = zero_matrix(R, 1, ngens(N))
     add_one!(m, 1, i)
   else
     m = zero_matrix(R, ngens(N), 1)
     add_one!(m, i, 1)
   end
   return N(m)
end

basis(N::FreeModule) = gens(N)

Base.hash(a::FreeModuleElem, h::UInt) = hash(a.v, h)

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, M::FreeModule{T}) where T <: NCRingElement
   @show_name(io, M)
   @show_special(io, M)
   if M.is_row
     print(io, "Free module of rank ")
   else
     print(io, "Free column-module of rank ")
   end
   print(io, rank(M))
   print(io, " over ")
   print(terse(pretty(io)), Lowercase(), base_ring(M))
end

function show(io::IO, M::FreeModule{T}) where T <: FieldElement
   @show_name(io, M)
   @show_special(io, M)
   print(io, "Vector space of dimension ")
   print(io, dim(M))
   print(io, " over ")
   print(terse(pretty(io)), Lowercase(), base_ring(M))
end

function show(io::IO, a::FreeModuleElem)
   print(io, "(")
   M = parent(a)
   if M.is_row
     for i = 1:rank(M) - 1
        print(IOContext(io, :compact => true), _matrix(a)[1, i])
        print(io, ", ")
     end
     if rank(M) > 0
        print(IOContext(io, :compact => true), _matrix(a)[1, rank(M)])
     end
   else
     for i = 1:rank(M) - 1
        print(IOContext(io, :compact => true), _matrix(a)[i, 1])
        print(io, ", ")
     end
     if rank(M) > 0
        print(IOContext(io, :compact => true), _matrix(a)[rank(M), 1])
     end
   end

   print(io, ")")
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (M::FreeModule{T})(a::Vector{T}) where T <: NCRingElement
   length(a) != rank(M) && error("Number of elements does not equal rank")
   R = base_ring(M)
   if M.is_row
     v = matrix(R, 1, length(a), a)
   else
     v = matrix(R, length(a), 1, a)
   end
   z = FreeModuleElem{T}(M, v)
   return z
end

#function (M::FreeModule{T})(a::Vector{<: Integer}) where T <: NCRingElement
#   length(a) != rank(M) && error("Number of elements does not equal rank")
#   R = base_ring(M)
#   v = matrix(R, 1, length(a), a)
#   z = FreeModuleElem{T}(M, v)
#   return z
#end

function (M::FreeModule{T})(a::Vector{S}) where {T <: NCRingElement, S <: RingElement}
   length(a) != rank(M) && error("Number of elements does not equal rank")
   R = base_ring(M)
   if M.is_row
     v = matrix(R, 1, length(a), a)
   else
     v = matrix(R, length(a), 1, a)
   end
   z = FreeModuleElem{T}(M, v)
   return z
end

function (M::FreeModule{T})(a::Vector{Any}) where T <: NCRingElement
   length(a) != 0 && error("Incompatible element")
   return M(T[])
end

function (M::FreeModule{T})(a::AbstractAlgebra.MatElem{T}) where T <: NCRingElement
  if M.is_row
     ncols(a) != rank(M) && error("Number of elements does not equal rank")
     nrows(a) != 1 && error("Matrix should have single row")
   else
     nrows(a) != rank(M) && error("Number of elements does not equal rank")
     ncols(a) != 1 && error("Matrix should have single column")
   end
   z = FreeModuleElem{T}(M, a)
   return z
end

function (M::FreeModule{T})(a::SubmoduleElem{T}) where T <: RingElement
   R = parent(a)
   base_ring(R) !== base_ring(M) && error("Incompatible modules")
   return M(R.map(a))
end

function (M::FreeModule{T})(a::FreeModuleElem{T}) where T <: RingElement
   R = parent(a)
   R !== M && error("Incompatible modules")
   return a
end

# Fallback for all other kinds of modules
function (M::FreeModule{T})(a::AbstractAlgebra.FPModuleElem{T}) where T <: RingElement
   error("Unable to coerce into given module")
end

###############################################################################
#
#   FreeModule constructor
#
###############################################################################

function FreeModule(R::NCRing, rank::Int; cached::Bool = true, is_row::Bool = true)
   T = elem_type(R)
   return FreeModule{T}(R, rank, cached; is_row)
end

