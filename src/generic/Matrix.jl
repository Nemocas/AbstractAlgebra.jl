###############################################################################
#
#   Matrix.jl : generic matrices over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{S}) where {T <: NCRingElement, S <: Mat{T}} = MatSpace{T}

elem_type(::Type{MatSpace{T}}) where {T <: NCRingElement} = MatSpaceElem{T}

@doc raw"""
    parent(a::AbstractAlgebra.MatElem{T}) where T <: NCRingElement

Return the parent object of the given matrix.
"""
parent(a::Mat{T}) where T <: NCRingElement = MatSpace{T}(a.base_ring, size(a.entries)...)

@doc raw"""
    dense_matrix_type(::Type{T}) where T<:NCRingElement
    dense_matrix_type(::T) where T<:NCRingElement
    dense_matrix_type(::Type{S}) where S<:NCRing
    dense_matrix_type(::S) where S<:NCRing

Return the type of matrices with coefficients of type `T` respectively
`elem_type(S)`.
"""
dense_matrix_type(::T) where T <: NCRing = dense_matrix_type(elem_type(T))
dense_matrix_type(::T) where T <: NCRingElement = dense_matrix_type(T)
dense_matrix_type(::Type{T}) where T <: NCRing = dense_matrix_type(elem_type(T))

# default: MatSpaceElem
dense_matrix_type(::Type{T}) where T <: NCRingElement = MatSpaceElem{T}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

@doc raw"""
    nrows(a::MatSpace)

Return the number of rows of the given matrix space.
"""
nrows(a::MatSpace) = a.nrows

@doc raw"""
    ncols(a::MatSpace)

Return the number of columns of the given matrix space.
"""
ncols(a::MatSpace) = a.ncols

nrows(a::Union{Mat, MatAlgElem}) = size(a.entries, 1)

ncols(a::Union{Mat, MatAlgElem}) = size(a.entries, 2)

Base.@propagate_inbounds getindex(a::Union{Mat, MatAlgElem}, r::Int, c::Int) = a.entries[r, c]

Base.@propagate_inbounds function setindex!(a::Union{Mat, MatAlgElem}, d::NCRingElement,
                                            r::Int, c::Int)
    a.entries[r, c] = base_ring(a)(d)
end

Base.isassigned(a::Union{Mat,MatAlgElem}, i, j) = isassigned(a.entries, i, j)

################################################################################
#
#  Copy and deepcopy
#
################################################################################

function copy(d::MatSpaceElem{T}) where T <: NCRingElement
   z = similar(d)
   for i = 1:nrows(d)
      for j = 1:ncols(d)
         z[i, j] = d[i, j]
      end
   end
   return z
end

function deepcopy_internal(d::MatSpaceElem{T}, dict::IdDict) where T <: NCRingElement
   z = similar(d)
   for i = 1:nrows(d)
      for j = 1:ncols(d)
         z[i, j] = deepcopy_internal(d[i, j], dict)
      end
   end
   return z
end

function deepcopy_internal(d::MatSpaceView{T}, dict::IdDict) where T <: NCRingElement
   return MatSpaceView(deepcopy_internal(d.entries, dict), d.base_ring)
end

###############################################################################
#
#   getindex
#
###############################################################################

# linear indexing for row- or column- vectors
Base.@propagate_inbounds function getindex(M::MatElem, x::Integer)
   if nrows(M) == 1
      M[1, x]
   elseif ncols(M) == 1
      M[x, 1]
   else
      throw(ArgumentError("linear indexing not supported for non-vector matrices"))
   end
end

function Base.view(M::Mat{T}, rows::UnitRange{Int}, cols::UnitRange{Int}) where T <: NCRingElement
   return MatSpaceView(view(M.entries, rows, cols), M.base_ring)
end

################################################################################
#
#   Size, axes and is_square
#
################################################################################

is_square(a::MatElem) = (nrows(a) == ncols(a))

###############################################################################
#
#   Transpose
#
###############################################################################

function transpose(x::Mat{T}) where T <: NCRingElement
   MatSpaceElem{eltype(x)}(base_ring(x), permutedims(x.entries))
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{S}, ::Type{S}) where {T <: NCRingElement, S <: Mat{T}} = MatSpaceElem{T}

function promote_rule(::Type{S}, ::Type{U}) where {T <: NCRingElement, S <: Mat{T}, U <: NCRingElement}
   promote_rule(T, U) == T ? MatSpaceElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::MatSpace{T})() where {T <: NCRingElement}
   R = base_ring(a)
   entries = Matrix{T}(undef, a.nrows, a.ncols)
   for i = 1:a.nrows
      for j = 1:a.ncols
         entries[i, j] = zero(R)
      end
   end
   z = MatSpaceElem{T}(R, entries)
   return z
end

function (a::MatSpace{T})(b::S) where {S <: NCRingElement, T <: NCRingElement}
   R = base_ring(a)
   entries = Matrix{T}(undef, a.nrows, a.ncols)
   rb = R(b)
   for i = 1:a.nrows
      for j = 1:a.ncols
         if i != j
            entries[i, j] = zero(R)
         else
            entries[i, j] = rb
         end
      end
   end
   z = MatSpaceElem{T}(R, entries)
   return z
end

function (a::MatSpace{T})(b::Matrix{T}) where T <: NCRingElement
   R = base_ring(a)
   _check_dim(a.nrows, a.ncols, b)
   if !isempty(b)
      R != parent(b[1, 1]) && error("Unable to coerce matrix")
   end
   z = MatSpaceElem{T}(R, b)
   return z
end

function (a::MatSpace{T})(b::AbstractMatrix{S}) where {S <: NCRingElement, T <: NCRingElement}
   R = base_ring(a)
   _check_dim(a.nrows, a.ncols, b)
   entries = Matrix{T}(undef, a.nrows, a.ncols)
   for i = 1:a.nrows
      for j = 1:a.ncols
         entries[i, j] = R(b[i, j])
      end
   end
   z = MatSpaceElem{T}(R, entries)
   return z
end

function (a::MatSpace{T})(b::AbstractVector{S}) where {S <: NCRingElement, T <: NCRingElement}
   _check_dim(a.nrows, a.ncols, b)
   b = Matrix{S}(transpose(reshape(b, a.ncols, a.nrows)))
   z = a(b)
   return z
end

###############################################################################
#
#   matrix_space constructor
#
###############################################################################

function matrix_space(R::AbstractAlgebra.NCRing, r::Int, c::Int; cached::Bool = true)
   # TODO: the 'cached' argument is ignored and mainly here for backwards compatibility
   # (and perhaps future compatibility, in case we need it again)
   T = elem_type(R)
   return MatSpace{T}(R, r, c)
end
