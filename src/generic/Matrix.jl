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
number_of_rows(a::Generic.MatSpace) = a.nrows

@doc raw"""
    number_of_columns(a::MatSpace)

Return the number of columns of the given matrix space.
"""
number_of_columns(a::Generic.MatSpace) = a.ncols

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
function (a::Generic.MatSpace{T})() where {T <: NCRingElement}
   return zero_matrix(base_ring(a), nrows(a), ncols(a))::dense_matrix_type(T)
end

# create a matrix with b on the diagonal
function (a::Generic.MatSpace)(b::NCRingElement)
   M = a()  # zero matrix
   R = base_ring(a)
   rb = R(b)
   for i in 1:min(nrows(a), ncols(a))
      M[i, i] = rb
   end
   return M
end

# convert a Julia matrix
function (a::Generic.MatSpace{T})(b::AbstractMatrix{S}) where {T <: NCRingElement, S}
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
function (a::Generic.MatSpace{T})(b::AbstractVector) where T <: NCRingElement
   _check_dim(nrows(a), ncols(a), b)
   return a(transpose(reshape(b, a.ncols, a.nrows)))
end

###############################################################################
#
#   matrix_space constructor
#
###############################################################################

function matrix_space(R::AbstractAlgebra.NCRing, r::Int, c::Int; cached::Bool = true)
   # TODO: the 'cached' argument is ignored and mainly here for backwards compatibility
   # (and perhaps future compatibility, in case we need it again)
   (r < 0 || c < 0) && error("Dimensions must be non-negative")
   T = elem_type(R)
   return Generic.MatSpace{T}(R, r, c)
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

###############################################################################
#
#   InjProjMat
#
###############################################################################

function inj_proj_mat(R::NCRing, r::Int, c::Int, s::Int)
   @assert r >= 0 && c >= 0 && s > 0
   # Check whether there is space for a full identity matrix
   if r <= c
      @assert s + r - 1 <= c
   else
      @assert s + c - 1 <= r
   end
   return InjProjMat{elem_type(R)}(R, r, c, s)
end

AbstractAlgebra.nrows(K::InjProjMat) = K.n
AbstractAlgebra.ncols(K::InjProjMat) = K.m
AbstractAlgebra.base_ring_type(::Type{InjProjMat{T}}) where T = parent_type(T)
AbstractAlgebra.base_ring(K::InjProjMat{T}) where T = K.R::parent_type(T)

function AbstractAlgebra.matrix(K::InjProjMat)
  R = base_ring(K)
  if nrows(K) >= ncols(K)
    return [zero_matrix(R, K.s-1, ncols(K)) ; identity_matrix(R, ncols(K)) ; zero_matrix(R, nrows(K) - K.s - ncols(K) + 1, ncols(K))]
  else
    return [zero_matrix(R, nrows(K), K.s-1) identity_matrix(R, nrows(K)) zero_matrix(R, nrows(K), ncols(K)-K.s-nrows(K) + 1)]
  end
end

function Base.getindex(K::InjProjMat{T}, i::Int, j::Int) where T
  (1 <= i <= nrows(K) && 1 <= j <= ncols(K)) || error(BoundsError(K, (i, j)))
  nrows(K) >= ncols(K) && i - K.s + 1 == j && return one(base_ring(K))::T
  nrows(K) <= ncols(K) && i == j - K.s + 1 && return one(base_ring(K))::T
  return zero(base_ring(K))::T
end

function *(b::InjProjMat{T}, c::MatElem{T}) where {T <: NCRingElement}
  @assert ncols(b) == nrows(c)
  R = base_ring(b)
  @assert base_ring(c) === R
  if nrows(b) >= ncols(b)
    z = zero_matrix(R, nrows(b), ncols(c))
    z[b.s:b.s+nrows(c)-1, :] = c
    return z
  else
    return c[b.s:b.s+nrows(b)-1, :]
  end
end

function *(b::MatElem{T}, c::InjProjMat{T}) where {T <: NCRingElement}
  @assert ncols(b) == nrows(c)
  R = base_ring(b)
  @assert base_ring(c) === R
  if nrows(c) >= ncols(c)
    #c = [0 I 0]^t
    return b[:, c.s:c.s+ncols(c)-1]
  else
    z = zero_matrix(R, nrows(b), ncols(c))
    z[:, c.s:c.s+nrows(c)-1] = b
    return z
  end
end

function *(a::InjProjMat, b::InjProjMat)
   @assert base_ring(a) === base_ring(b)
   R = base_ring(a)
   @assert ncols(a) == nrows(b)
   m = nrows(a)
   n = ncols(a)
   k = ncols(b)
   s = a.s
   t = b.s

   # Easy cases
   if m >= n
      return reduce(vcat, [zero_matrix(R, s - 1, k), matrix(b), zero_matrix(R, m - n - s + 1, k)])
   end
   if m < n && n < k
      return reduce(hcat, [zero_matrix(R, m, t - 1), matrix(a), zero_matrix(R, m, k - n - t + 1)])
   end

   # Annoying case: m < n && n >= k
   c = zero_matrix(R, m, k)
   if s <= t
      offset = t - s
      for i in 1:min(m - offset, k)
         c[i + offset, i] = one(R)
      end
   else
      offset = s - t
      for i in 1:min(m, k - offset)
         c[i, i + offset] = one(R)
      end
   end
   return c
end

function +(b::MatElem{T}, c::InjProjMat{T}) where {T <: NCRingElement}
  @assert size(b) == size(c)
  R = base_ring(b)
  @assert base_ring(c) === R
  a = deepcopy(b)
  if nrows(c) >= ncols(c)
    for i in 1:ncols(c)
      add_one!(a, c.s+i-1, i)
    end
  else
    for i in 1:nrows(c)
      add_one!(a, i, c.s+i-1)
    end
  end
  return a
end
+(c::InjProjMat{T}, b::MatElem{T}) where {T <: NCRingElement} = b+c
+(c::InjProjMat{T}, b::InjProjMat{T}) where {T <: NCRingElement} = matrix(b) + c

function add_one!(a::MatElem, i::Int, j::Int)
  a[i, j] += 1
  return a
end
