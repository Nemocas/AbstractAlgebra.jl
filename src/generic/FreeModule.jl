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

parent_type(::Type{FreeModuleElem{T}}) where T <: Union{RingElement, NCRingElem} = FreeModule{T}

base_ring_type(::Type{FreeModule{T}}) where T <: Union{RingElement, NCRingElem} = parent_type(T)

base_ring(M::FreeModule{T}) where T <: Union{RingElement, NCRingElem} = M.base_ring::parent_type(T)

elem_type(::Type{FreeModule{T}}) where T <: Union{RingElement, NCRingElem} = FreeModuleElem{T}

parent(m::FreeModuleElem{T}) where T <: Union{RingElement, NCRingElem} = m.parent

function rels(M::FreeModule{T}) where T <: RingElement
   # there are no relations in a free module
   return Vector{dense_matrix_type(T)}(undef, 0)
end

function check_parent(m1::FreeModuleElem{T}, m2::FreeModuleElem{T}) where T <: Union{RingElement, NCRingElem}
    parent(m1) !== parent(m2) && error("Incompatible free modules")
end

@doc raw"""
    rank(M::FreeModule{T}) where T <: Union{RingElement, NCRingElem}

Return the rank of the given free module.
"""
rank(M::FreeModule{T}) where T <: Union{RingElement, NCRingElem} = M.rank

@doc raw"""
    dim(M::FreeModule{T}) where T <: FieldElement

Return the dimension of the given vector space.
"""
dim(M::FreeModule{T}) where T <: FieldElement = M.rank

number_of_generators(M::FreeModule{T}) where T <: Union{RingElement, NCRingElem} = M.rank

function gens(N::FreeModule{T}) where T <: Union{RingElement, NCRingElem}
   return [gen(N, i) for i = 1:ngens(N)]
end

function gen(N::FreeModule{T}, i::Int) where T <: Union{RingElement, NCRingElem}
   @boundscheck 1 <= i <= ngens(N) || throw(ArgumentError("generator index out of range"))
   R = base_ring(N)
   m = zero_matrix(R, 1, ngens(N))
   m[1, i] = one(R)
   return N(m)
end

basis(N::FreeModule) = gens(N)

Base.hash(a::FreeModuleElem, h::UInt) = hash(a.v, h)

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, M::FreeModule{T}) where T <: Union{RingElement, NCRingElem}
   AbstractAlgebra.@show_name(io, M)
   AbstractAlgebra.@show_special(io, M)

   print(io, "Free module of rank ")
   print(io, rank(M))
   print(io, " over ")
   io = pretty(io)
   print(IOContext(io, :compact => true), Lowercase(), base_ring(M))
end

function show(io::IO, M::FreeModule{T}) where T <: FieldElement
   AbstractAlgebra.@show_name(io, M)
   AbstractAlgebra.@show_special(io, M)

   print(io, "Vector space of dimension ")
   print(io, dim(M))
   print(io, " over ")
   io = pretty(io)
   print(IOContext(io, :compact => true), Lowercase(), base_ring(M))
end

function show(io::IO, a::FreeModuleElem)
   print(io, "(")
   M = parent(a)
   for i = 1:rank(M) - 1
      print(IOContext(io, :compact => true), _matrix(a)[1, i])
      print(io, ", ")
   end
   if rank(M) > 0
      print(IOContext(io, :compact => true), _matrix(a)[1, rank(M)])
   end
   print(io, ")")
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (M::FreeModule{T})(a::Vector{T}) where T <: Union{RingElement, NCRingElem}
   length(a) != rank(M) && error("Number of elements does not equal rank")
   R = base_ring(M)
   v = matrix(R, 1, length(a), a)
   z = FreeModuleElem{T}(M, v)
   return z
end

#function (M::FreeModule{T})(a::Vector{<: Integer}) where T <: Union{RingElement, NCRingElem}
#   length(a) != rank(M) && error("Number of elements does not equal rank")
#   R = base_ring(M)
#   v = matrix(R, 1, length(a), a)
#   z = FreeModuleElem{T}(M, v)
#   return z
#end

function (M::FreeModule{T})(a::Vector{S}) where {T <: Union{RingElement, NCRingElem}, S <: RingElement}
   length(a) != rank(M) && error("Number of elements does not equal rank")
   R = base_ring(M)
   v = matrix(R, 1, length(a), a)
   z = FreeModuleElem{T}(M, v)
   return z
end

function (M::FreeModule{T})(a::Vector{Any}) where T <: Union{RingElement, NCRingElem}
   length(a) != 0 && error("Incompatible element")
   return M(T[])
end

function (M::FreeModule{T})(a::AbstractAlgebra.MatElem{T}) where T <: Union{RingElement, NCRingElem}
   ncols(a) != rank(M) && error("Number of elements does not equal rank")
   nrows(a) != 1 && error("Matrix should have single row")
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

function FreeModule(R::NCRing, rank::Int; cached::Bool = true)
   T = elem_type(R)
   return FreeModule{T}(R, rank, cached)
end

