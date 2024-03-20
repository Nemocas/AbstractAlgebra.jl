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

parent_type(::Type{<:MatElem{T}}) where {T <: NCRingElement} = MatSpace{T}

elem_type(::Type{MatSpace{T}}) where {T <: NCRingElement} = dense_matrix_type(T)

@doc raw"""
    parent(a::AbstractAlgebra.MatElem)

Return the parent object of the given matrix.
"""
parent(a::MatElem) = matrix_space(base_ring(a), nrows(a), ncols(a))

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
    number_of_rows(a::MatSpace)

Return the number of rows of the given matrix space.
"""
number_of_rows(a::MatSpace) = a.nrows

@doc raw"""
    number_of_columns(a::MatSpace)

Return the number of columns of the given matrix space.
"""
number_of_columns(a::MatSpace) = a.ncols

number_of_rows(a::Union{Mat, MatRingElem}) = size(a.entries, 1)

number_of_columns(a::Union{Mat,MatRingElem}) = size(a.entries, 2)

Base.@propagate_inbounds getindex(a::Union{Mat, MatRingElem}, r::Int, c::Int) = a.entries[r, c]

Base.@propagate_inbounds function setindex!(a::Union{Mat, MatRingElem}, d::NCRingElement,
                                            r::Int, c::Int)
    a.entries[r, c] = base_ring(a)(d)
end

Base.isassigned(a::Union{Mat,MatRingElem}, i, j) = isassigned(a.entries, i, j)

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

function Base.view(M::Mat{T}, rows::Union{Colon, AbstractVector{Int}}, cols::Union{Colon, AbstractVector{Int}}) where T <: NCRingElement
   return MatSpaceView(view(M.entries, rows, cols), M.base_ring)
end

function Base.view(M::Mat{T}, rows::Int, cols::Union{Colon, AbstractVector{Int}}) where T <: NCRingElement
   return MatSpaceVecView(view(M.entries, rows, cols), M.base_ring)
end

function Base.view(M::Mat{T}, rows::Union{Colon, AbstractVector{Int}}, cols::Int) where T <: NCRingElement
   return MatSpaceVecView(view(M.entries, rows, cols), M.base_ring)
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

# create a zero matrix
function (a::MatSpace{T})() where {T <: NCRingElement}
   return zero_matrix(base_ring(a), nrows(a), ncols(a))::dense_matrix_type(T)
end

# create a matrix with b on the diagonal
function (a::AbstractAlgebra.Generic.MatSpace)(b::NCRingElement)
   M = a()  # zero matrix
   R = base_ring(a)
   rb = R(b)
   for i in 1:min(nrows(a), ncols(a))
      M[i, i] = rb
   end
   return M
end

# convert a Julia matrix
function (a::MatSpace{T})(b::AbstractMatrix{S}) where {T <: NCRingElement, S}
   _check_dim(nrows(a), ncols(a), b)
   R = base_ring(a)

   # minor optimization for MatSpaceElem
   if S === T && dense_matrix_type(T) === MatSpaceElem{T} && all(x -> R === parent(x), b)
      return MatSpaceElem{T}(R, b)
   end

   # generic code
   M = a()  # zero matrix
   for i = 1:nrows(a), j = 1:ncols(a)
      M[i, j] = R(b[i, j])
   end
   return M
end

# convert a Julia vector
function (a::MatSpace{T})(b::AbstractVector) where T <: NCRingElement
   _check_dim(nrows(a), ncols(a), b)
   return a(transpose(reshape(b, a.ncols, a.nrows)))
end

###############################################################################
#
#   matrix_space constructor
#
###############################################################################

function matrix_space(R::AbstractAlgebra.NCRing, r::Int, c::Int)
   T = elem_type(R)
   return MatSpace{T}(R, r, c)
end

function AbstractAlgebra.sub!(A::Mat{T}, B::Mat{T}, C::Mat{T}) where T
  A.entries.= B.entries .- C.entries
  return A
end

#since type(view(MatElem{T})) != MatElem{T} which breaks
# sub!(A::T, B::T, C::T) where T  in AA 
function AbstractAlgebra.mul!(A::Mat{T}, B::Mat{T}, C::Mat{T}, f::Bool = false) where T
  if f
    A.entries .+= (B * C).entries
  else
    A.entries .= (B * C).entries
  end
  return A
end

Base.length(V::MatSpaceVecView) = length(V.entries)

Base.getindex(V::MatSpaceVecView, i::Int) = V.entries[i]

Base.setindex!(V::MatSpaceVecView{T}, z::T, i::Int) where {T} = (V.entries[i] = z)

Base.setindex!(V::MatSpaceVecView, z::RingElement, i::Int) = setindex!(V.entries, V.base_ring(z), i)

Base.size(V::MatSpaceVecView) = (length(V.entries), )
