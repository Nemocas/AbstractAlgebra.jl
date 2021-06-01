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

parent_type(::Type{S}) where {T <: RingElement, S <: Mat{T}} = MatSpace{T}

elem_type(::Type{MatSpace{T}}) where {T <: RingElement} = MatSpaceElem{T}

@doc Markdown.doc"""
    parent(a::AbstractAlgebra.MatElem{T}, cached::Bool = true) where T <: RingElement

Return the parent object of the given matrix.
"""
parent(a::Mat{T}, cached::Bool = true) where T <: RingElement =
    MatSpace{T}(a.base_ring, size(a.entries)..., cached)

dense_matrix_type(::Type{T}) where T <: RingElement = MatSpaceElem{T}

@doc Markdown.doc"""
    dense_matrix_type(R::Ring)
    
Return the type of matrices over the given ring.
"""
dense_matrix_type(R::Ring) = dense_matrix_type(elem_type(R))

###############################################################################
#
#   Basic manipulation
#
###############################################################################

@doc Markdown.doc"""
    nrows(a::MatSpace)

Return the number of rows of the given matrix space.
"""
nrows(a::MatSpace) = a.nrows

@doc Markdown.doc"""
    ncols(a::MatSpace)

Return the number of columns of the given matrix space.
"""
ncols(a::MatSpace) = a.ncols

nrows(a::Union{Mat, MatAlgElem}) = size(a.entries, 1)

ncols(a::Union{Mat, MatAlgElem}) = size(a.entries, 2)

Base.@propagate_inbounds getindex(a::Union{Mat, MatAlgElem}, r::Int, c::Int) = a.entries[r, c]

Base.@propagate_inbounds function setindex!(a::Union{Mat, MatAlgElem}, d::RingElement,
                                            r::Int, c::Int)
    a.entries[r, c] = base_ring(a)(d)
end

Base.isassigned(a::Union{Mat,MatAlgElem}, i, j) = isassigned(a.entries, i, j)

################################################################################
#
#  Copy and deepcopy
#
################################################################################

function copy(d::MatSpaceElem{T}) where T <: RingElement
   z = similar(d)
   for i = 1:nrows(d)
      for j = 1:ncols(d)
         z[i, j] = d[i, j]
      end
   end
   return z
end

function deepcopy_internal(d::MatSpaceElem{T}, dict::IdDict) where T <: RingElement
   z = similar(d)
   for i = 1:nrows(d)
      for j = 1:ncols(d)
         z[i, j] = deepcopy(d[i, j])
      end
   end
   return z
end

function deepcopy_internal(d::MatSpaceView{T}, dict::IdDict) where T <: RingElement
   return MatSpaceView(deepcopy(d.entries), d.base_ring)
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

function Base.view(M::Mat{T}, rows::UnitRange{Int}, cols::UnitRange{Int}) where T <: RingElement
   return MatSpaceView(view(M.entries, rows, cols), M.base_ring)
end

################################################################################
#
#   Size, axes and issquare
#
################################################################################

issquare(a::MatElem) = (nrows(a) == ncols(a))

###############################################################################
#
#   Transpose
#
###############################################################################

@doc Markdown.doc"""
    transpose(x::Mat)

Return the transpose of the given matrix.
"""
function transpose(x::Mat)
   y = MatSpaceElem{eltype(x)}(permutedims(x.entries))
   y.base_ring = x.base_ring
   y
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{S}, ::Type{S}) where {T <: RingElement, S <: Mat{T}} = MatSpaceElem{T}

function promote_rule(::Type{S}, ::Type{U}) where {T <: RingElement, S <: Mat{T}, U <: RingElement}
   promote_rule(T, U) == T ? MatSpaceElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::MatSpace{T})() where {T <: RingElement}
   R = base_ring(a)
   entries = Array{T}(undef, a.nrows, a.ncols)
   for i = 1:a.nrows
      for j = 1:a.ncols
         entries[i, j] = zero(R)
      end
   end
   z = MatSpaceElem{T}(entries)
   z.base_ring = R
   return z
end

function (a::MatSpace{T})(b::S) where {S <: RingElement, T <: RingElement}
   R = base_ring(a)
   entries = Array{T}(undef, a.nrows, a.ncols)
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
   z = MatSpaceElem{T}(entries)
   z.base_ring = R
   return z
end

function (a::MatSpace{T})(b::Mat{T}) where {T <: RingElement}
   parent(b) != a && error("Unable to coerce matrix")
   return b
end

function (a::MatSpace{T})(b::Array{T, 2}) where T <: RingElement
   R = base_ring(a)
   _check_dim(a.nrows, a.ncols, b)
   if !isempty(b)
      R != parent(b[1, 1]) && error("Unable to coerce matrix")
   end
   z = MatSpaceElem{T}(b)
   z.base_ring = R
   return z
end

function (a::MatSpace{T})(b::AbstractArray{S, 2}) where {S <: RingElement, T <: RingElement}
   R = base_ring(a)
   _check_dim(a.nrows, a.ncols, b)
   entries = Array{T}(undef, a.nrows, a.ncols)
   for i = 1:a.nrows
      for j = 1:a.ncols
         entries[i, j] = R(b[i, j])
      end
   end
   z = MatSpaceElem{T}(entries)
   z.base_ring = R
   return z
end

function (a::MatSpace{T})(b::AbstractArray{S, 1}) where {S <: RingElement, T <: RingElement}
   _check_dim(a.nrows, a.ncols, b)
   b = Array{S, 2}(transpose(reshape(b, a.ncols, a.nrows)))
   z = a(b)
   return z
end

###############################################################################
#
#   MatrixSpace constructor
#
###############################################################################

function MatrixSpace(R::AbstractAlgebra.Ring, r::Int, c::Int; cached::Bool = true)
   T = elem_type(R)
   return MatSpace{T}(R, r, c, cached)
end