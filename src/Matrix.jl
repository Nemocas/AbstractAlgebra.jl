###############################################################################
#
#   Matrix.jl : matrices over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

elem_type(::Type{MatSpace{T}}) where {T <: NCRingElement} = dense_matrix_type(T)

base_ring_type(::Type{MatSpace{T}}) where T <: NCRingElement = parent_type(T)

base_ring(a::MatSpace{T}) where {T <: NCRingElement} = a.base_ring::parent_type(T)

base_ring_type(::Type{<:MatElem{T}}) where T <: NCRingElement = parent_type(T)

parent_type(::Type{<:MatElem{T}}) where {T <: NCRingElement} = MatSpace{T}

@doc raw"""
    parent(a::MatElem)

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

Implementations of the ring interface only need to provide a method
for the argument a subtype of `NCRingElement`; the other variants are
implemented by calling that method.
"""
dense_matrix_type(::T) where T <: NCRing = dense_matrix_type(elem_type(T))
dense_matrix_type(::T) where T <: NCRingElement = dense_matrix_type(T)
dense_matrix_type(::Type{T}) where T <: NCRing = dense_matrix_type(elem_type(T))

# default: Generic.MatSpaceElem
dense_matrix_type(::Type{T}) where T <: NCRingElement = Generic.MatSpaceElem{T}


###############################################################################
#
#   Various helpers
#
###############################################################################

"""
    is_zero_initialized(T::Type{<:MatrixElem})
    is_zero_initialized(mat::T) where {T<:MatrixElem}

Specify whether the default-constructed matrices of type `T`, via the
`T(R::Ring, ::UndefInitializer, r::Int, c::Int)` constructor, are
zero-initialized. The default is `false`, and new matrix types should
specialize this method to return `true` if suitable, to enable optimizations.
"""
is_zero_initialized(::Type{<:MatrixElem}) = false
is_zero_initialized(::T) where {T<:MatrixElem} = is_zero_initialized(T)

function check_parent(a::MatrixElem, b::MatrixElem, throw::Bool = true)
  fl = (base_ring(a) != base_ring(b) || nrows(a) != nrows(b) || ncols(a) != ncols(b))
  fl && throw && error("Incompatible matrix spaces in matrix operation")
  return !fl
end

function _check_dim(r::Int, c::Int, arr::AbstractMatrix{T}, transpose::Bool = false) where {T}
  Base.require_one_based_indexing(arr)
  if !transpose
    size(arr) != (r, c) && throw(ErrorConstrDimMismatch(r, c, size(arr)...))
  else
    size(arr) != (c, r) && throw(ErrorConstrDimMismatch(r, c, (reverse(size(arr)))...))
  end
  return nothing
end

function _check_dim(r::Int, c::Int, arr::AbstractVector{T}) where {T}
  Base.require_one_based_indexing(arr)
  length(arr) != r*c && throw(ErrorConstrDimMismatch(r, c, length(arr)))
  return nothing
end

function _check_dim(r::Int, c::Int, a::MatrixElem)
  size(a) == (r, c) || throw(ErrorConstrDimMismatch(r, c, size(a)...))
  return nothing
end

_checkbounds(i::Int, j::Int) = 1 <= j <= i

_checkbounds(i::Int, j::AbstractVector{Int}) = all(jj -> 1 <= jj <= i, j)

function _checkbounds(A, i::Union{Int, AbstractVector{Int}}, j::Union{Int, AbstractVector{Int}})
  (_checkbounds(nrows(A), i) && _checkbounds(ncols(A), j)) ||
            Base.throw_boundserror(A, (i, j))
end

function _checkbounds(A, rows::AbstractArray{Int}, cols::AbstractArray{Int})
   Base.has_offset_axes(rows, cols) && throw(ArgumentError("offset arrays are not supported"))
   isempty(rows) || _checkbounds(nrows(A), first(rows)) && _checkbounds(nrows(A), last(rows)) ||
      throw(BoundsError(A, rows))
   isempty(cols) || _checkbounds(ncols(A), first(cols)) && _checkbounds(ncols(A), last(cols)) ||
      throw(BoundsError(A, cols))
end

function check_square(A::MatrixElem{T}) where T <: NCRingElement
   is_square(A) || throw(DomainError(A, "matrix must be square"))
   A
end

function check_square(S::MatSpace)
   nrows(S) == ncols(S) || throw(DomainError(S, "matrices must be square"))
   S
end

check_square(S::MatRing) = S

###############################################################################
#
#   Parent object call overload
#
###############################################################################

# create a zero matrix
function (s::MatSpace{T})() where {T <: NCRingElement}
  return zero_matrix(base_ring(s), nrows(s), ncols(s))::eltype(s)
end

function (s::MatSpace{T})(a::MatrixElem{T}) where {T <: NCRingElement}
  _check_dim(nrows(s), ncols(s), a)
  base_ring(s) == base_ring(a) || throw(DomainError((s, a), "Base rings do not match."))
  a isa eltype(s) && return a
  return matrix(base_ring(s), a)
end

# create a matrix with b on the diagonal
function (s::MatSpace)(b::NCRingElement)
  R = base_ring(s)
  return diagonal_matrix(R(b), nrows(s), ncols(s))
end

# convert a Julia matrix or vector
function (a::MatSpace{T})(b::AbstractVecOrMat) where T <: NCRingElement
  return matrix(base_ring(a), nrows(a), ncols(a), b)
end


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

vector_space_dim(a::MatSpace{T}) where {T <: Union{FieldElem, Rational{BigInt}}} = a.nrows * a.ncols

function Base.hash(a::MatElem, h::UInt)
   b = 0x3e4ea81eb31d94f4%UInt
   for i in 1:nrows(a)
      for j in 1:ncols(a)
         b = xor(b, xor(hash(a[i, j], h), h))
         b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
      end
   end
   return b
end

@doc raw"""
    number_of_rows(a::MatrixElem{T}) where T <: NCRingElement

Return the number of rows of the given matrix.
"""
number_of_rows(a::MatrixElem{T}) where T <: NCRingElement

@doc raw"""
    number_of_columns(a::MatrixElem{T}) where T <: NCRingElement

Return the number of columns of the given matrix.
"""
number_of_columns(a::MatrixElem{T}) where {T<:NCRingElement}

@doc raw"""
    length(a::MatrixElem{T}) where T <: NCRingElement

Return the number of entries in the given matrix.
"""
length(a::MatrixElem{T}) where T <: NCRingElement = nrows(a) * ncols(a)

@doc raw"""
    isempty(a::MatrixElem{T}) where T <: NCRingElement

Return `true` if `a` does not contain any entry (i.e. `length(a) == 0`), and `false` otherwise.
"""
isempty(a::MatrixElem{T}) where T <: NCRingElement = (nrows(a) == 0) | (ncols(a) == 0)

Base.eltype(::Type{<:MatrixElem{T}}) where {T <: NCRingElement} = T

function Base.isassigned(a::MatrixElem{T}, i, j) where T <: NCRingElement
    try
        a[i, j]
        true
    catch e
        if isa(e, BoundsError) || isa(e, UndefRefError)
            return false
        else
            rethrow()
        end
    end
end

@doc raw"""
    zero(a::MatSpace)

Return the zero matrix in the given matrix space.
"""
zero(a::MatSpace) = a()

@doc raw"""
    one(a::MatSpace)

Return the identity matrix of given matrix space. The matrix space must contain
square matrices or else an error is thrown.
"""
one(a::MatSpace) = check_square(a)(1)

@doc raw"""
    one(a::MatrixElem{T}) where T <: NCRingElement

Return the identity matrix in the same matrix space as $a$. If the space does
not contain square matrices, an error is thrown.
"""
one(a::MatrixElem{T}) where T <: NCRingElement = identity_matrix(a)

function iszero(a::MatrixElem{T}) where T <: NCRingElement
   for i = 1:nrows(a)
      for j = 1:ncols(a)
         if !is_zero_entry(a, i, j)
            return false
         end
      end
  end
  return true
end

function isone(a::MatrixElem{T}) where T <: NCRingElement
   is_square(a) || return false
   for i = 1:nrows(a)
      for j = 1:ncols(a)
         if i == j
            if !isone(a[i, j])
               return false
            end
         else
            if !is_zero_entry(a, i, j)
               return false
            end
         end
      end
   end
   return true
end

@doc raw"""
    is_zero_entry(M::Union{Matrix,MatrixElem}, i::Int, j::Int)

Return `true` if $M_{i,j}$ is zero.
"""
@inline is_zero_entry(M::Union{Matrix,MatrixElem}, i::Int, j::Int) = iszero(M[i,j])

@doc raw"""
    is_positive_entry(M::Union{Matrix,MatrixElem}, i::Int, j::Int)

Return `is_positive(M[i,j])`, but with a possibly more efficient implementation.
"""
@inline is_positive_entry(M::Union{Matrix,MatrixElem}, i::Int, j::Int) = is_positive(M[i,j])

@doc raw"""
    is_negative_entry(M::Union{Matrix,MatrixElem}, i::Int, j::Int)

Return `is_negative(M[i,j])`, but with a possibly more efficient implementation.
"""
@inline is_negative_entry(M::Union{Matrix,MatrixElem}, i::Int, j::Int) = is_negative(M[i,j])

@doc raw"""
    is_zero_row(M::Union{Matrix,MatrixElem}, i::Int)

Return `true` if the $i$-th row of the matrix $M$ is zero.
"""
function is_zero_row(M::Union{Matrix,MatrixElem}, i::Int)
  @boundscheck 1 <= i <= nrows(M) || Base.throw_boundserror(M, (i, 1:ncols(M)))
  for j in 1:ncols(M)
    @inbounds if !is_zero_entry(M, i, j)
      return false
    end
  end
  return true
end

@doc raw"""
    is_zero_column(M::Union{Matrix,MatrixElem}, j::Int)

Return `true` if the $j$-th column of the matrix $M$ is zero.
"""
function is_zero_column(M::Union{Matrix,MatrixElem}, j::Int)
  @boundscheck 1 <= j <= ncols(M) || Base.throw_boundserror(M, (1:nrows(M), j))
  for i in 1:nrows(M)
    @inbounds if !is_zero_entry(M, i, j)
      return false
    end
  end
  return true
end

###############################################################################
#
#   Block diagonal matrices    
#
###############################################################################

@doc raw"""
    block_diagonal_matrix(V::Vector{<:MatElem{T}}) where T <: NCRingElement

Create the block diagonal matrix whose blocks are given by the matrices in `V`.
There must be at least one matrix in V.
"""
function block_diagonal_matrix(V::Vector{<:MatElem{T}}) where T <: NCRingElement
   @req !isempty(V) "Cannot infer base ring from empty vector; consider passing the desired base ring as first argument to `block_diagonal_matrix`"
   rows = sum(nrows(N) for N in V)
   cols = sum(ncols(N) for N in V)
   R = base_ring(V[1])
   M = similar(V[1], rows, cols)
   start_row = 1
   start_col = 1
   for i = 1:length(V)
      end_row = start_row + nrows(V[i]) - 1
      end_col = start_col + ncols(V[i]) - 1
      for j = start_row:end_row
         for k = 1:start_col - 1
            M[j, k] = zero(R)
         end
         for k = start_col:end_col
            M[j, k] = V[i][j - start_row + 1, k - start_col + 1]
         end
         for k = end_col + 1:cols
            M[j, k] = zero(R)
         end
      end
      start_row = end_row + 1
      start_col = end_col + 1
   end
   return M
end

@doc raw"""
    block_diagonal_matrix(R::NCRing, V::Vector{<:Matrix{T}}) where T <: NCRingElement

Create the block diagonal matrix over the ring `R` whose blocks are given
by the matrices in `V`. Entries are coerced into `R` upon creation.
"""
function block_diagonal_matrix(R::NCRing, V::Vector{<:Matrix{T}}) where T <: NCRingElement
   if length(V) == 0
      return zero_matrix(R, 0, 0)
   end
   rows = sum(size(N)[1] for N in V)
   cols = sum(size(N)[2] for N in V)
   M = zero_matrix(R, rows, cols)
   start_row = 1
   start_col = 1
   for i = 1:length(V)
      end_row = start_row + size(V[i])[1] - 1
      end_col = start_col + size(V[i])[2] - 1
      for j = start_row:end_row
         for k = start_col:end_col
            M[j, k] = R(V[i][j - start_row + 1, k - start_col + 1])
         end
      end
      start_row = end_row + 1
      start_col = end_col + 1
   end
   return M
end

###############################################################################
#
#   Similar
#
###############################################################################

@doc raw"""
    similar(x::MatElem{T}, R::NCRing, r::Int, c::Int) where T <: NCRingElement
    similar(x::MatElem{T}, R::NCRing) where T <: NCRingElement
    similar(x::MatElem{T}, r::Int, c::Int) where T <: NCRingElement
    similar(x::MatElem{T}) where T <: NCRingElement

Create an uninitialized matrix over the given ring and dimensions,
with defaults based upon the given source matrix `x`.
"""
similar(x::MatElem, R::NCRing, r::Int, c::Int) = dense_matrix_type(R)(R, undef, r, c)
  
similar(x::MatElem, R::NCRing) = similar(x, R, nrows(x), ncols(x))

similar(x::MatElem, r::Int, c::Int) = similar(x, base_ring(x), r, c)

similar(x::MatElem) = similar(x, nrows(x), ncols(x))

@doc raw"""
    zero(x::MatElem{T}, R::NCRing, r::Int, c::Int) where T <: NCRingElement
    zero(x::MatElem{T}, r::Int, c::Int) where T <: NCRingElement
    zero(x::MatElem{T}, R::NCRing) where T <: NCRingElement
    zero(x::MatElem{T}) where T <: NCRingElement

Create an zero matrix over the given ring and dimensions,
with defaults based upon the given source matrix `x`.
"""
zero(x::MatElem{T}, R::NCRing) where T <: NCRingElement = zero(x, R, nrows(x), ncols(x))
zero(x::MatElem{T}) where T <: NCRingElement = zero(x, nrows(x), ncols(x))

function zero(x::MatElem{T}, R::NCRing, r::Int, c::Int) where T <: NCRingElement
  y = similar(x, R, r, c)
  return is_zero_initialized(y) ? y : zero!(y)
end

function zero(x::MatElem{T}, r::Int, c::Int) where T <: NCRingElement
  y = similar(x, r, c)
  return is_zero_initialized(y) ? y : zero!(y)
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::MatrixElem{T}) where T <: NCRingElement = canonical_unit(a[1, 1])

###############################################################################
#
#   getindex
#
###############################################################################

@doc raw"""
    Base.getindex(M::MatElem, rows, cols)

When `rows` and `cols` are specified as an `AbstractVector{Int}`, return a copy of
the submatrix $A$ of $M$ defined by `A[i,j] = M[rows[i], cols[j]]`
for `i=1,...,length(rows)` and `j=1,...,length(cols)`.
Instead of a vector, `rows` and `cols` can also be:
* an integer `i`, which is  interpreted as `i:i`, or
* `:`, which is interpreted as `1:nrows(M)` or `1:ncols(M)` respectively.
"""
getindex(M::MatElem, r::AbstractVector{<:Integer}, c::AbstractVector{<:Integer}) = sub(M, r, c)

function getindex(M::MatElem, i::Int, cols::AbstractVector{Int})
   _checkbounds(M, i, cols)
   A = Vector{elem_type(base_ring(M))}(undef, length(cols))
   for j in eachindex(cols)
     A[j] = deepcopy(M[i, cols[j]])
   end
   return A
end

function getindex(M::MatElem, rows::AbstractVector{Int}, j::Int)
   _checkbounds(M, rows, j)
   A = Vector{elem_type(base_ring(M))}(undef, length(rows))
   for i in eachindex(rows)
     A[i] = deepcopy(M[rows[i], j])
   end
   return A
 end

getindex(M::MatElem, ::Colon, cols) = getindex(M, 1:nrows(M), cols)

getindex(M::MatElem, rows, ::Colon) = getindex(M, rows, 1:ncols(M))

getindex(M::MatElem, ::Colon, ::Colon) = getindex(M, 1:nrows(M), 1:ncols(M))

function sub(M::MatElem, rows::AbstractVector{Int}, cols::AbstractVector{Int})
   _checkbounds(M, rows, cols)
   A = similar(M, length(rows), length(cols))
   for i in 1:length(rows)
      for j in 1:length(cols)
         A[i, j] = deepcopy(M[rows[i], cols[j]])
      end
   end
   return A
end

# fallback method that converts Colons to UnitRanges
function Base.view(M::MatElem, rows, cols)
   # indirection to avoid ambiguities
   return _view(M, rows, cols)
end

function _view(M::MatElem, ::Colon, cols)
   return view(M, 1:nrows(M), cols)
end

function _view(M::MatElem, rows, ::Colon)
   return view(M, rows, 1:ncols(M))
end

function _view(M::MatElem, ::Colon, ::Colon)
   return view(M, 1:nrows(M), 1:ncols(M))
end

Base.firstindex(M::MatrixElem{T}, i::Int) where T <: NCRingElement = 1

function Base.lastindex(M::MatrixElem{T}, i::Int) where T <: NCRingElement
   if i == 1
      return nrows(M)
   elseif i == 2
      return ncols(M)
   else
      error("Dimension in lastindex must be 1 or 2 (got $i)")
   end
end

###############################################################################
#
#   Array interface
#
###############################################################################

Base.ndims(::MatrixElem{T}) where T <: NCRingElement = 2

# Cartesian indexing

Base.eachindex(a::MatrixElem{T}) where T <: NCRingElement = CartesianIndices((nrows(a), ncols(a)))

Base.@propagate_inbounds Base.getindex(a::MatrixElem{T}, I::CartesianIndex) where T <: NCRingElement =
   a[I[1], I[2]]

Base.@propagate_inbounds function Base.setindex!(a::MatrixElem{T}, x, I::CartesianIndex) where T <: NCRingElement
   a[I[1], I[2]] = x
   a
end

# linear indexing for row- or column- vectors
Base.@propagate_inbounds function getindex(M::MatrixElem, i::Integer)
   if nrows(M) == 1
      M[1, i]
   elseif ncols(M) == 1
      M[i, 1]
   else
      throw(ArgumentError("linear indexing not supported for non-vector matrices"))
   end
end

Base.@propagate_inbounds function setindex!(M::MatrixElem, x, i::Integer)
   if nrows(M) == 1
      M[1, i] = x
      return M
   elseif ncols(M) == 1
      M[i, 1] = x
      return M
   else
      throw(ArgumentError("linear indexing not supported for non-vector matrices"))
   end
end

# iteration

function Base.iterate(a::MatrixElem{T}, ij=(0, 1)) where T <: NCRingElement
   i, j = ij
   i += 1
   if i > nrows(a)
      iszero(nrows(a)) && return nothing
      i = 1
      j += 1
   end
   j > ncols(a) && return nothing
   a[i, j], (i, j)
end

Base.IteratorSize(::Type{<:MatrixElem}) = Base.HasShape{2}()

Base.keys(M::MatElem) = CartesianIndices(axes(M))
Base.pairs(M::MatElem) = Base.pairs(IndexCartesian(), M)
Base.pairs(::IndexCartesian, M::MatElem) = Base.Iterators.Pairs(M, CartesianIndices(axes(M)))

###############################################################################
#
#   Block replacement
#
###############################################################################

function setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix}, r::AbstractUnitRange{Int}, c::AbstractUnitRange{Int}) where T <: NCRingElement
    _checkbounds(a, r, c)
    size(b) == (length(r), length(c)) || throw(DimensionMismatch("tried to assign a $(size(b, 1))x$(size(b, 2)) matrix to a $(length(r))x$(length(c)) destination"))
    startr = first(r)
    startc = first(c)
    for i in r
        for j in c
            a[i, j] = b[i - startr + 1, j - startc + 1]
        end
    end
end

function setindex!(a::MatrixElem{T}, b::Vector, r::AbstractUnitRange{Int}, c::AbstractUnitRange{Int}) where T <: NCRingElement
    _checkbounds(a, r, c)
    if !((length(r) == 1 && length(c) == length(b)) || length(c) == 1 && length(r) == length(b))
      throw(DimensionMismatch("tried to assign vector of length $(length(b)) to a $(length(r))x$(length(c)) destination"))
    end
    startr = first(r)
    startc = first(c)
    for i in r
        for j in c
            a[i, j] = b[i - startr + 1 + j - startc]
        end
    end
end

# AbstractUnitRange{Int}, Colon
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::AbstractUnitRange{Int}, ::Colon) where T <: NCRingElement = setindex!(a, b, r, 1:ncols(a))

# Colon, AbstractUnitRange{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, ::Colon, c::AbstractUnitRange{Int}) where T <: NCRingElement = setindex!(a, b, 1:nrows(a), c)

# Colon, Colon
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, ::Colon, ::Colon) where T <: NCRingElement = setindex!(a, b, 1:nrows(a), 1:ncols(a))

# Int, AbstractUnitRange{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Int, c::AbstractUnitRange{Int}) where T <: NCRingElement = setindex!(a, b, r:r, c)

# AbstractUnitRange{Int}, Int
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::AbstractUnitRange{Int}, c::Int) where T <: NCRingElement = setindex!(a, b, r, c:c)

# Int, Colon
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Int, ::Colon) where T <: NCRingElement = setindex!(a, b, r:r, 1:ncols(a))

# Colon, Int
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, ::Colon, c::Int) where T <: NCRingElement = setindex!(a, b, 1:nrows(a), c:c)

function _setindex!(a::MatrixElem{T}, b, r, c) where T <: NCRingElement
   for (i, i2) in enumerate(r)
      for (j, j2) in enumerate(c)
         a[i2, j2] = b[i, j]
      end
   end
end

function _setindex!(a::MatrixElem{T}, b::Vector, r, c) where T <: NCRingElement
   for (i, i2) in enumerate(r)
      for (j, j2) in enumerate(c)
         a[i2, j2] = b[i + j - 1]
      end
   end
end

# Vector{Int}, Vector{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Vector{Int}, c::Vector{Int}) where T <: NCRingElement = _setindex!(a, b, r, c)

# Vector{Int}, AbstractUnitRange{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Vector{Int}, c::AbstractUnitRange{Int}) where T <: NCRingElement = _setindex!(a, b, r, c)

# AbstractUnitRange{Int}, Vector{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::AbstractUnitRange{Int}, c::Vector{Int}) where T <: NCRingElement = _setindex!(a, b, r, c)

# Vector{Int}, Colon
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Vector{Int}, ::Colon) where T <: NCRingElement = _setindex!(a, b, r, 1:ncols(a))

# Colon, Vector{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, ::Colon, c::Vector{Int}) where T <: NCRingElement = _setindex!(a, b, 1:nrows(a), c)

# Int, Vector{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Int, c::Vector{Int}) where T <: NCRingElement = setindex!(a, b, r:r, c)

# Vector{Int}, Int
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Vector{Int}, c::Int) where T <: NCRingElement = setindex!(a, b, r, c:c)

################################################################################
#
#   Size, axes and is_square
#
################################################################################

size(x::MatrixElem{T}) where T <: NCRingElement = (nrows(x), ncols(x))

size(t::MatrixElem{T}, d::Integer) where T <: NCRingElement = d <= 2 ? size(t)[d] : 1

axes(t::MatrixElem{T}) where T <: NCRingElement = Base.OneTo.(size(t))

axes(t::MatrixElem{T}, d::Integer) where T <: NCRingElement = Base.OneTo(size(t, d))

###############################################################################
#
#   eachrow / eachcol
#
###############################################################################

Base.eachrow(a::MatrixElem) = Slices(a, (1, :), (axes(a, 1),))
Base.eachcol(a::MatrixElem) = Slices(a, (:, 1), (axes(a, 2),))

###############################################################################
#
#   Matrix spaces iteration
#
###############################################################################

function Base.iterate(M::MatSpace)
   R = base_ring(M)
   d = nrows(M) * ncols(M)
   p = ProductIterator(fill(R, d); inplace=true)
   if d == 0
      # handle this carefully to preserve type stability
      a = elem_type(R)[]
      state_type = typeof(iterate(R)[2])
      st = (elem_type(R)[], state_type[])
   else
      a, st = iterate(p)::Tuple{Any, Any} # R is presumably not empty
   end
   M(a), (p, st)
end

function Base.iterate(M::MatSpace, (p, st))
   nrows(M) * ncols(M) == 0 && return nothing
   a_st = iterate(p, st)
   a_st === nothing && return nothing
   M(first(a_st)), (p, last(a_st))
end

Base.eltype(::Type{M}) where {M<:MatSpace} = elem_type(M)
Base.length(M::MatSpace) = BigInt(length(base_ring(M)))^(nrows(M)*ncols(M))

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::MatrixElem{T}; context = nothing) where T <: NCRingElement
   r = nrows(a)
   c = ncols(a)
   isempty(a) && return "$r by $c empty matrix"
   mat = Expr(:vcat)
   for i in 1:r
      row = Expr(:row)
      for j in 1:c
         if isassigned(a, i, j)
            push!(row.args, expressify(a[i, j], context = context))
         else
            push!(row.args,Base.undef_ref_str)
         end
      end
      push!(mat.args, row)
   end
   return mat
end

@enable_all_show_via_expressify MatrixElem

function Base.show(io::IO, ::MIME"text/plain", a::MatrixElem{T}) where T <: NCRingElement
   r = nrows(a)
   c = ncols(a)

   if isempty(a)
      print(io, "$r by $c empty matrix")
      return
   end

   # preprint each element to know the widths so as to align the columns
   strings = String[sprint(print, isassigned(a, i, j) ? a[i, j] : Base.undef_ref_str,
                           context = :compact => true) for i=1:r, j=1:c]
   maxs = maximum(length, strings, dims=1)

   for i = 1:r
      print(io, "[")
      for j = 1:c
         s = strings[i, j]
         s = ' '^(maxs[j] - length(s)) * s
         print(io, s)
         if j != c
            print(io, "   ")
         end
      end
      print(io, "]")
      if i != r
         println(io, "")
      end
   end
end

function show(io::IO, ::MIME"text/plain", a::MatSpace)
  print(io, "Matrix space of ")
  print(io, ItemQuantity(nrows(a), "row"), " and ", ItemQuantity(ncols(a), "column"))
  println(io)
  io = pretty(io)
  print(io, Indent(), "over ")
  print(io, Lowercase(), base_ring(a))
  print(io, Dedent())
end

function show(io::IO, a::MatSpace)
   if is_terse(io)
      print(io, "Matrix space")
   else
      io = pretty(io)
      print(io, "Matrix space of ")
      print(io, ItemQuantity(nrows(a), "row"), " and ", ItemQuantity(ncols(a), "column"))
      print(io, " over ")
      print(terse(io), Lowercase(), base_ring(a))
   end
end

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(x::MatrixElem{T}) where T <: NCRingElement
   z = similar(x)
   for i in 1:nrows(x)
      for j in 1:ncols(x)
         z[i, j] = -x[i, j]
      end
   end
   return z
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(x::T, y::T) where {T <: MatElem}
   check_parent(x, y)
   r = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         r[i, j] = x[i, j] + y[i, j]
      end
   end
   return r
end

function -(x::T, y::T) where {T <: MatElem}
   check_parent(x, y)
   r = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         r[i, j] = x[i, j] - y[i, j]
      end
   end
   return r
end

function *(x::MatElem{T}, y::MatElem{T}) where {T}
   ncols(x) != nrows(y) && error("Incompatible matrix dimensions")
   base_ring(x) !== base_ring(y) && error("Base rings do not match")
   A = similar(x, nrows(x), ncols(y))
   C = base_ring(x)()
   for i = 1:nrows(x)
      for j = 1:ncols(y)
         A[i, j] = base_ring(x)()
         for k = 1:ncols(x)
            A[i, j] = addmul_delayed_reduction!(A[i, j], x[i, k], y[k, j], C)
         end
         A[i, j] = reduce!(A[i, j])
      end
   end
   return A
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(x::MatrixElem{T}) where T <: NCRingElement
   R = base_ring(x)
   for i = 1:nrows(x), j = 1:ncols(x)
      x[i, j] = zero(R)
   end
   return x
end

function add!(c::MatrixElem{T}, a::MatrixElem{T}, b::MatrixElem{T}) where T <: NCRingElement
   check_parent(a, b)
   check_parent(a, c)
   for i = 1:nrows(c)
      for j = 1:ncols(c)
         c[i, j] = add!(c[i, j], a[i, j], b[i, j])
      end
   end
   return c
end

function mul!(c::MatElem{T}, a::MatElem{T}, b::MatElem{T}) where T <: NCRingElement
   @assert base_ring(a) === base_ring(b) && base_ring(a) === base_ring(c)
   ncols(a) != nrows(b) && error("Incompatible matrix dimensions")
   nrows(c) != nrows(a) && error("Incompatible matrix dimensions")
   ncols(c) != ncols(b) && error("Incompatible matrix dimensions")

   if c === a || c === b
      d = parent(a)()
      return mul!(d, a, b)
   end

   t = base_ring(a)()
   for i = 1:nrows(a)
      for j = 1:ncols(b)
         c[i, j] = zero!(c[i, j])
         for k = 1:ncols(a)
            c[i, j] = addmul_delayed_reduction!(c[i, j], a[i, k], b[k, j], t)
         end
         c[i, j] = reduce!(c[i, j])
      end
   end
   return c
end

function sub!(c::MatrixElem{T}, a::MatrixElem{T}, b::MatrixElem{T}) where T <: NCRingElement
   check_parent(a, b)
   check_parent(a, c)
   for i = 1:nrows(c)
      for j = 1:ncols(c)
         c[i, j] = sub!(c[i, j], a[i, j], b[i, j])
      end
   end
   return c
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(x::JuliaRingElement, y::MatrixElem{T}) where T <: NCRingElement
   z = similar(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         z[i, j] = x*y[i, j]
      end
   end
   return z
end

function *(x::T, y::MatrixElem{T}) where {T <: NCRingElem}
   z = similar(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         z[i, j] = x*y[i, j]
      end
   end
   return z
end

function *(x::MatrixElem{T}, y::JuliaRingElement) where T <: NCRingElement
   z = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         z[i, j] = x[i, j]*y
      end
   end
   return z
end

function *(x::MatrixElem{T}, y::T) where {T <: NCRingElem}
   z = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         z[i, j] = x[i, j]*y
      end
   end
   return z
end

function +(x::JuliaRingElement, y::MatrixElem{T}) where T <: NCRingElement
   z = similar(y)
   R = base_ring(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         if i != j
            z[i, j] = deepcopy(y[i, j])
         else
            z[i, j] = y[i, j] + R(x)
         end
      end
   end
   return z
end

+(x::MatrixElem{T}, y::JuliaRingElement) where T <: NCRingElement = y + x

@doc raw"""
    +(x::NCRingElement, y::MatrixElem{<:NCRingElement})

Return $S(x) + y$ where $S$ is the parent of $y$.
"""
function +(x::T, y::MatrixElem{T}) where {T <: NCRingElem}
   z = similar(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         if i != j
            z[i, j] = deepcopy(y[i, j])
         else
            z[i, j] = y[i, j] + x
         end
      end
   end
   return z
end

@doc raw"""
    +(x::MatrixElem{<:NCRingElement}, y::NCRingElement)

Return $x + S(y)$, where $S$ is the parent of $a$.
"""
+(x::MatrixElem{T}, y::T) where {T <: NCRingElem} = y + x

function -(x::JuliaRingElement, y::MatrixElem{T}) where T <: NCRingElement
   z = similar(y)
   R = base_ring(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         if i != j
            z[i, j] = -y[i, j]
         else
            z[i, j] = R(x) - y[i, j]
         end
      end
   end
   return z
end

function -(x::MatrixElem{T}, y::JuliaRingElement) where T <: NCRingElement
   z = similar(x)
   R = base_ring(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if i != j
            z[i, j] = deepcopy(x[i, j])
         else
            z[i, j] = x[i, j] - R(y)
         end
      end
   end
   return z
end

@doc raw"""
    -(x::NCRingElement, y::MatrixElem{<:NCRingElement})

Return $S(x) - y$ where $S$ is the parent of $y$.
"""
function -(x::T, y::MatrixElem{T}) where {T <: NCRingElem}
   z = similar(y)
   R = base_ring(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         if i != j
            z[i, j] = -y[i, j]
         else
            z[i, j] = x - y[i, j]
         end
      end
   end
   return z
end

@doc raw"""
    -(x::MatrixElem{<:NCRingElem}, y::NCRingElement)

Return $x - S(y)$, where $S$ is the parent of $a$.
"""
function -(x::MatrixElem{T}, y::T) where {T <: NCRingElem}
   z = similar(x)
   R = base_ring(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if i != j
            z[i, j] = deepcopy(x[i, j])
         else
            z[i, j] = x[i, j] - y
         end
      end
   end
   return z
end

function mul!(z::Vector{T}, x::MatrixElem{T}, y::Vector{T}) where T <: NCRingElement
   n = min(ncols(x), length(y))
   tmp = base_ring(x)()
   for i in 1:nrows(x)
      if n > 0
         z[i] = mul!(z[i], x[i, 1], y[1])
         for j in 2:n
            z[i] = addmul_delayed_reduction!(z[i], x[i, j], y[j], tmp)
         end
         z[i] = reduce!(z[i])
      else
         z[i] = zero!(z[i])
      end
   end
   return z
end

function *(x::MatrixElem{T}, y::Vector{T}) where T <: NCRingElement
   ncols(x) == length(y) || error("Incompatible dimensions")
   return mul!(T[base_ring(x)() for i in 1:nrows(x)], x, y)
end

function mul!(z::Vector{T}, x::Vector{T}, y::MatrixElem{T}) where T <: NCRingElement
   m = min(length(x), nrows(y))
   tmp = base_ring(y)()
   for j in 1:ncols(y)
      if m > 0
         z[j] = mul!(z[j], x[1], y[1, j])
         for i in 2:m
            z[j] = addmul_delayed_reduction!(z[j], x[i], y[i, j], tmp)
         end
         z[j] = reduce!(z[j])
      else
         z[j] = zero!(z[j])
      end
   end
   return z
end

function *(x::Vector{T}, y::MatrixElem{T}) where T <: NCRingElement
   length(x) == nrows(y) || error("Incompatible dimensions")
   return mul!(T[base_ring(y)() for j in 1:ncols(y)], x, y)
end

function mul!(c::MatrixElem{T}, a::MatrixElem{T}, b::T) where T <: NCRingElement
   @assert base_ring(a) === parent(b) && base_ring(a) === base_ring(c)
   nrows(c) != nrows(a) && error("Incompatible matrix dimensions")
   ncols(c) != ncols(a) && error("Incompatible matrix dimensions")

   if c === a
      d = parent(a)()
      return mul!(d, a, b)
   end

   t = base_ring(a)()
   for i = 1:nrows(a)
      for j = 1:ncols(a)
         c[i, j] = mul!(c[i, j], a[i, j], b)
      end
   end
   return c
end

################################################################################
#
#  Promotion
#
################################################################################

function Base.promote(x::MatrixElem{S},
                      y::MatrixElem{T}) where {S <: NCRingElement,
                                               T <: NCRingElement}
   U = promote_rule_sym(S, T)
   if U === S
      return x, change_base_ring(base_ring(x), y)
   elseif U === T
      return change_base_ring(base_ring(y), x), y
   else
      error("Cannot promote to common type")
   end
end

*(x::MatElem, y::MatElem) = *(promote(x, y)...)

+(x::MatElem, y::MatElem) = +(promote(x, y)...)

-(x::MatElem, y::MatElem) = -(promote(x, y)...)

==(x::MatElem, y::MatElem) = ==(promote(x, y)...)

# matrix * vec and vec * matrx
function Base.promote(x::MatrixElem{S},
                      y::Vector{T}) where {S <: NCRingElement,
                                               T <: NCRingElement}
   U = promote_rule_sym(S, T)
   if U === S
      return x, map(base_ring(x), y)::Vector{S}  # Julia needs help here
   elseif U === T && length(y) != 0
      return change_base_ring(parent(y[1]), x), y
   else
      error("Cannot promote to common type")
   end
end

function Base.promote(x::Vector{S},
                      y::MatrixElem{T}) where {S <: NCRingElement,
                                               T <: NCRingElement}
   yy, xx = promote(y, x)
   return xx, yy
end

*(x::MatrixElem, y::Vector) = *(promote(x, y)...)

*(x::Vector, y::MatrixElem) = *(promote(x, y)...)

function Base.promote(x::MatElem{S}, y::T) where {S <: NCRingElement, T <: NCRingElement}
   U = promote_rule_sym(S, T)
   if U === S
      return x, base_ring(x)(y)
   elseif U === T
      return change_base_ring(parent(y), x), y
   else
      error("Cannot promote to common type")
   end
end

function Base.promote(x::S, y::MatElem{T}) where {S <: NCRingElement, T <: NCRingElement}
   u, v = Base.promote(y, x)
   return v, u
end

*(x::MatElem, y::NCRingElem) = *(promote(x, y)...)

*(x::NCRingElem, y::MatElem) = *(promote(x, y)...)

+(x::MatElem, y::NCRingElem) = +(promote(x, y)...)

+(x::NCRingElem, y::MatElem) = +(promote(x, y)...)

==(x::MatElem, y::NCRingElem) = ==(promote(x, y)...)

==(x::NCRingElem, y::MatElem) = ==(promote(x, y)...)

divexact(x::MatElem, y::NCRingElem; check::Bool = true) =
    divexact(promote(x, y)...; check = check)

divexact_left(x::MatElem, y::NCRingElem; check::Bool = true) =
    divexact_left(promote(x, y)...; check = check)

divexact_right(x::MatElem, y::NCRingElem; check::Bool = true) =
    divexact_right(promote(x, y)...; check = check)

###############################################################################
#
#   Powering
#
###############################################################################

Base.literal_pow(::typeof(^), x::T, ::Val{p}) where {p, U <: NCRingElement, T <: MatrixElem{U}} = x^p

@doc raw"""
    ^(a::MatrixElem{T}, b::Int) where T <: NCRingElement

Return $a^b$. We require that the matrix $a$ is square.
"""
function ^(a::MatrixElem{T}, b::Int) where T <: NCRingElement
   !is_square(a) && error("Incompatible matrix dimensions in power")
   if b < 0
      return inv(a)^(-b)
   end
   # special case powers of x for constructing polynomials efficiently
   if b == 0
      return identity_matrix(a)
   elseif b == 1
      return deepcopy(a)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit != 0
         z = z*z
         if (UInt(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end

###############################################################################
#
#   Comparisons
#
###############################################################################

@doc raw"""
    ==(x::MatrixElem{T}, y::MatrixElem{T}) where {T <: NCRingElement}

Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(x::MatrixElem{T}, y::MatrixElem{T}) where {T <: NCRingElement}
   b = check_parent(x, y, false)
   !b && return false
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if x[i, j] != y[i, j]
            return false
         end
      end
   end
   return true
end

@doc raw"""
    isequal(x::MatrixElem{T}, y::MatrixElem{T}) where {T <: NCRingElement}

Return `true` if $x == y$ exactly, otherwise return `false`. This function is
useful in cases where the entries of the matrices are inexact, e.g. power
series. Only if the power series are precisely the same, to the same precision,
are they declared equal by this function.
"""
function isequal(x::MatrixElem{T}, y::MatrixElem{T}) where {T <: NCRingElement}
   b = check_parent(x, y, false)
   !b && return false
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if !isequal(x[i, j], y[i, j])
            return false
         end
      end
   end
   return true
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

function ==(x::MatrixElem{T}, y::JuliaRingElement) where T <: NCRingElement
   for i = 1:min(nrows(x), ncols(x))
      if x[i, i] != y
         return false
      end
   end
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if i != j && !is_zero_entry(x, i, j)
            return false
         end
      end
   end
   return true
end

==(x::JuliaRingElement, y::MatrixElem{T}) where T <: NCRingElement = y == x

@doc raw"""
    ==(x::MatrixElem{<:NCRingElement}, y::NCRingElement)

Return `true` if $x == S(y)$ arithmetically, where $S$ is the parent of $x$,
otherwise return `false`.
"""
function ==(x::MatrixElem{T}, y::T) where {T <: NCRingElem}
   for i = 1:min(nrows(x), ncols(x))
      if x[i, i] != y
         return false
      end
   end
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if i != j && !is_zero_entry(x, i, j)
            return false
         end
      end
   end
   return true
end

@doc raw"""
    ==(x::NCRingElement, y::MatrixElem{<:NCRingElement})

Return `true` if $S(x) == y$ arithmetically, where $S$ is the parent of $y$,
otherwise return `false`.
"""
==(x::T, y::MatrixElem{T}) where {T <: NCRingElem} = y == x

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::MatrixElem{T}, y::JuliaRingElement; check::Bool=true) where T <: NCRingElement
   z = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         z[i, j] = divexact(x[i, j], y; check=check)
      end
   end
   return z
end

function divexact(x::MatrixElem{T}, y::T; check::Bool=true) where {T <: RingElem}
   z = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         z[i, j] = divexact(x[i, j], y; check=check)
      end
   end
   return z
end

function divexact_left(x::MatrixElem{T}, y::T; check::Bool=true) where {T <: NCRingElem}
   z = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         z[i, j] = divexact_left(x[i, j], y; check=check)
      end
   end
   return z
end

function divexact_right(x::MatrixElem{T}, y::T; check::Bool=true) where {T <: NCRingElem}
   z = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         z[i, j] = divexact_right(x[i, j], y; check=check)
      end
   end
   return z
end

###############################################################################
#
#   Symmetry
#
###############################################################################

"""
    is_symmetric(M::MatrixElem)

Return `true` if the given matrix is symmetric with respect to its main
diagonal, i.e., `transpose(M) == M`, otherwise return `false`.

Alias for `LinearAlgebra.issymmetric`.

# Examples

```jldoctest
julia> M = matrix(ZZ, [1 2 3; 2 4 5; 3 5 6])
[1   2   3]
[2   4   5]
[3   5   6]

julia> is_symmetric(M)
true

julia> N = matrix(ZZ, [1 2 3; 4 5 6; 7 8 9])
[1   2   3]
[4   5   6]
[7   8   9]

julia> is_symmetric(N)
false
```
"""
function is_symmetric(M::MatrixElem)
   n = nrows(M)
   n == ncols(M) || return false
   for i in 2:n, j in 1:i-1
      M[i, j] == M[j, i] || return false
   end
   return true
end


@doc raw"""
    transpose(x::MatElem)
    transpose(x::MatRingElem)

Return the transpose of `x`.

# Examples

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> A = matrix(R, [t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t + 1       t             1]
[  t^2       t             t]
[   -2   t + 2   t^2 + t + 1]

julia> B = transpose(A)
[t + 1   t^2            -2]
[    t     t         t + 2]
[    1     t   t^2 + t + 1]

```
"""
function transpose(x::MatElem)
  z = similar(base_ring(x), ncols(x), nrows(x))
  return transpose!(z, x)
end

@doc raw"""
    transpose!(x::MatElem)
    transpose!(x::MatRingElem)
    transpose!(z::T, x::T) where T <: MatElem
    transpose!(z::T, x::T) where T <: MatRingElem

Return the transpose of `x`, possibly modifying the object `z` in the process.
Aliasing is permitted.

The unary version only is supported if `x` is a square matrix.
"""
function transpose!(x::MatElem)
  @req is_square(x) "Matrix must be a square matrix"
  return transpose!(x, x)
end

function transpose!(z::T, x::T) where T <: MatElem
  if z === x
    n = nrows(x)
    for i in 1:n, j in i+1:n
      x[i, j], x[j, i] = x[j, i], x[i, j]
    end
  else
    for i in 1:nrows(x), j in 1:ncols(x)
      z[j, i] = x[i, j]
    end
  end
  return z
end

###############################################################################
#
#   Kronecker product
#
###############################################################################

function kronecker_product(x::MatrixElem{T}, y::MatrixElem{T}) where {T <: RingElement}
    base_ring(parent(x)) == base_ring(parent(y)) || error("Incompatible matrix spaces in matrix operation")
    z = similar(x, nrows(x)*nrows(y), ncols(x)*ncols(y))
    for ix in 1:nrows(x)
       ixr = (ix - 1)*nrows(y)
       for jx in 1:ncols(x)
          jxc = (jx - 1)*ncols(y)
          for iy in 1:nrows(y)
             for jy in 1:ncols(y)
               z[ixr + iy, jxc + jy] = x[ix, jx]*y[iy, jy]
             end
          end
       end
    end
    return z
end

###############################################################################
#
#   Gram matrix
#
###############################################################################

@doc raw"""
    gram(x::MatElem)

Return the Gram matrix of $x$, i.e. if $x$ is an $r\times c$ matrix return
the $r\times r$ matrix whose entries $i, j$ are the dot products of the
$i$-th and $j$-th rows, respectively.

# Examples

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S = matrix_space(R, 3, 3)
Matrix space of 3 rows and 3 columns
  over univariate polynomial ring in t over rationals

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t + 1       t             1]
[  t^2       t             t]
[   -2   t + 2   t^2 + t + 1]

julia> B = gram(A)
[2*t^2 + 2*t + 2   t^3 + 2*t^2 + t                   2*t^2 + t - 1]
[t^3 + 2*t^2 + t       t^4 + 2*t^2                       t^3 + 3*t]
[  2*t^2 + t - 1         t^3 + 3*t   t^4 + 2*t^3 + 4*t^2 + 6*t + 9]

```
"""
function gram(x::MatElem)
   z = similar(x, nrows(x), nrows(x))
   for i = 1:nrows(x)
      for j = 1:nrows(x)
         z[i, j] = zero(base_ring(x))
         for k = 1:ncols(x)
            z[i, j] += x[i, k] * x[j, k]
         end
      end
   end
   return z
end

###############################################################################
#
#   Trace
#
###############################################################################

@doc raw"""
    tr(x::MatrixElem{T}) where T <: NCRingElement

Return the trace of the matrix $a$, i.e. the sum of the diagonal elements. We
require the matrix to be square.

# Examples

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S = matrix_space(R, 3, 3)
Matrix space of 3 rows and 3 columns
  over univariate polynomial ring in t over rationals

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t + 1       t             1]
[  t^2       t             t]
[   -2   t + 2   t^2 + t + 1]

julia> b = tr(A)
t^2 + 3*t + 2

```
"""
function tr(x::MatrixElem{T}) where T <: NCRingElement
   !is_square(x) && error("Not a square matrix in trace")
   d = zero(base_ring(x))
   for i = 1:nrows(x)
      d = add!(d, x[i, i])
   end
   return d
end

###############################################################################
#
#   Content
#
###############################################################################

@doc raw"""
    content(x::MatrixElem{T}) where T <: RingElement

Return the content of the matrix $a$, i.e. the greatest common divisor of all
its entries, assuming it exists.

# Examples

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S = matrix_space(R, 3, 3)
Matrix space of 3 rows and 3 columns
  over univariate polynomial ring in t over rationals

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t + 1       t             1]
[  t^2       t             t]
[   -2   t + 2   t^2 + t + 1]

julia> b = content(A)
1

```
"""
function content(x::MatrixElem{T}) where T <: RingElement
  d = zero(base_ring(x))
  for i = 1:nrows(x)
     for j = 1:ncols(x)
        d = gcd!(d, x[i, j])
     end
  end
  return d
end

###############################################################################
#
#   Permutation
#
###############################################################################

@doc raw"""
    *(P::Perm, x::MatrixElem{T}) where T <: NCRingElement

Apply the pemutation $P$ to the rows of the matrix $x$ and return the result.

# Examples

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S = matrix_space(R, 3, 3)
Matrix space of 3 rows and 3 columns
  over univariate polynomial ring in t over rationals

julia> G = SymmetricGroup(3)
Full symmetric group over 3 elements

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t + 1       t             1]
[  t^2       t             t]
[   -2   t + 2   t^2 + t + 1]

julia> P = G([1, 3, 2])
(2,3)

julia> B = P*A
[t + 1       t             1]
[   -2   t + 2   t^2 + t + 1]
[  t^2       t             t]

```
"""
function *(P::Perm, x::MatrixElem{T}) where T <: NCRingElement
   z = similar(x)
   m = nrows(x)
   n = ncols(x)
   for i = 1:m
      for j = 1:n
         z[P[i], j] = x[i, j]
      end
   end
   return z
end

@doc raw"""
    *(x::MatrixElem{T}, P::Perm) where T <: NCRingElement

Apply the pemutation $P$ to the columns of the matrix $x$ and return the result.

# Examples

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S = matrix_space(R, 3, 3)
Matrix space of 3 rows and 3 columns
  over univariate polynomial ring in t over rationals

julia> G = SymmetricGroup(3)
Full symmetric group over 3 elements

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t + 1       t             1]
[  t^2       t             t]
[   -2   t + 2   t^2 + t + 1]

julia> P = G([1, 3, 2])
(2,3)

julia> B = A*P
[t + 1             1       t]
[  t^2             t       t]
[   -2   t^2 + t + 1   t + 2]

```
"""
function *(x::MatrixElem{T}, P::Perm) where T <: NCRingElement
   z = similar(x)
   m = nrows(x)
   n = ncols(x)
   for i = 1:m
      for j = 1:n
         z[i, P[j]] = x[i, j]
      end
   end
   return z
end

###############################################################################
#
#   LU factorisation
#
###############################################################################

function lu!(P::Perm, A::MatrixElem{T}) where {T <: FieldElement}
   m = nrows(A)
   n = ncols(A)
   rank = 0
   r = 1
   c = 1
   R = base_ring(A)
   t = R()
   while r <= m && c <= n
      if c != 1
         # reduction of lower right square was delayed, reduce left col now
         for i = r:m
            A[i, c] = reduce!(A[i, c])
         end
      end
      if is_zero_entry(A, r, c)
         i = r + 1
         while i <= m
            if !is_zero_entry(A, i, c)
               for j = 1:n
                  A[i, j], A[r, j] = A[r, j], A[i, j]
               end
               P[r], P[i] = P[i], P[r]
               break
            end
            i += 1
         end
         if i > m
            c += 1
            continue
         end
      end
      rank += 1
      d = -inv(A[r, c])
      # reduction of lower right square was delayed, reduce top row now
      if c != 1
         for j = c + 1:n
            A[r, j] = reduce!(A[r, j])
         end
      end
      for i = r + 1:m
         q = A[i, c]*d
         for j = c + 1:n
            t = mul_red!(t, A[r, j], q, false)
            A[i, j] = add!(A[i, j], t)
         end
         A[i, c] = R()
         A[i, rank] = -q
      end
      r += 1
      c += 1
   end
   P = inv!(P)
   return rank
end

@doc raw"""
    lu(A::MatrixElem{T}, P = SymmetricGroup(nrows(A))) where {T <: FieldElement}

Return a tuple $r, p, L, U$ consisting of the rank of $A$, a permutation
$p$ of $A$ belonging to $P$, a lower triangular matrix $L$ and an upper
triangular matrix $U$ such that $p(A) = LU$, where $p(A)$ stands for the
matrix whose rows are the given permutation $p$ of the rows of $A$.
"""
function lu(A::MatrixElem{T}, P = SymmetricGroup(nrows(A))) where {T <: FieldElement}
   m = nrows(A)
   n = ncols(A)
   P.n != m && error("Permutation does not match matrix")
   p = one(P)
   R = base_ring(A)
   U = deepcopy(A)
   L = similar(A, m, m)
   rank = lu!(p, U)
   for i = 1:m
      for j = 1:n
         if i > j
            L[i, j] = U[i, j]
            U[i, j] = R()
         elseif i == j
            L[i, j] = one(R)
         elseif j <= m
            L[i, j] = R()
         end
      end
   end
   for i = 1:m
      for j = n + 1:m
         L[i, j] = R()
      end
   end
   return rank, p, L, U
end

function fflu!(P::Perm, A::MatrixElem{T}) where {T <: RingElement}
   if !is_domain_type(T)
      error("Not implemented")
   end
   m = nrows(A)
   n = ncols(A)
   rank = 0
   r = 1
   c = 1
   R = base_ring(A)
   d = one(R)
   d2 = one(R)
   if m == 0 || n == 0
      return 0, d
   end
   t = R()
   while r <= m && c <= n
      if is_zero_entry(A, r, c)
         i = r + 1
         while i <= m
            if !is_zero_entry(A, i, c)
               for j = 1:n
                  A[i, j], A[r, j] = A[r, j], A[i, j]
               end
               P[r], P[i] = P[i], P[r]
               break
            end
            i += 1
         end
         if i > m
            c += 1
            continue
         end
      end
      rank += 1
      q = -A[r, c]
      for i = r + 1:m
         for j = c + 1:n
            A[i, j] = mul_red!(A[i, j], A[i, j], q, false)
            t = mul_red!(t, A[i, c], A[r, j], false)
            A[i, j] = add!(A[i, j], t)
            A[i, j] = reduce!(A[i, j])
            if r > 1
               A[i, j] = divexact(A[i, j], d)
            else
               A[i, j] = -A[i, j]
            end
         end
      end
      d = -A[r, c]
      d2 = A[r, c]
      r += 1
      c += 1
   end
   P = inv!(P)
   return rank, d2
end

function fflu!(P::Perm, A::MatrixElem{T}) where {T <: Union{FieldElement, ResElem}}
   m = nrows(A)
   n = ncols(A)
   rank = 0
   r = 1
   c = 1
   R = base_ring(A)
   d = one(R)
   d2 = one(R)
   if m == 0 || n == 0
      return 0, d
   end
   t = R()
   while r <= m && c <= n
      if is_zero_entry(A, r, c)
         i = r + 1
         while i <= m
            if !is_zero_entry(A, i, c)
               for j = 1:n
                  A[i, j], A[r, j] = A[r, j], A[i, j]
               end
               P[r], P[i] = P[i], P[r]
               break
            end
            i += 1
         end
         if i > m
            c += 1
            continue
         end
      end
      rank += 1
      q = -A[r, c]
      for i = r + 1:m
         for j = c + 1:n
            A[i, j] = mul_red!(A[i, j], A[i, j], q, false)
            t = mul_red!(t, A[i, c], A[r, j], false)
            A[i, j] = add!(A[i, j], t)
            A[i, j] = reduce!(A[i, j])
            if r > 1
               A[i, j] = mul!(A[i, j], A[i, j], d)
            else
               A[i, j] = -A[i, j]
            end
         end
      end
      d = -inv(A[r, c])
      d2 = A[r, c]
      r += 1
      c += 1
   end
   P = inv!(P)
   return rank, d2
end

@doc raw"""
    fflu(A::MatrixElem{T}, P = SymmetricGroup(nrows(A))) where {T <: RingElement}

Return a tuple $r, d, p, L, U$ consisting of the rank of $A$, a
denominator $d$, a permutation $p$ of $A$ belonging to $P$, a lower
triangular matrix $L$ and an upper triangular matrix $U$ such that
$p(A) = LDU$, where $p(A)$ stands for the matrix whose rows are the
given permutation $p$ of the rows of $A$ and such that $D$ is the diagonal
matrix diag$(p_1, p_1p_2, \ldots, p_{n-2}p_{n-1}, p_{n-1}p_n)$ where the $p_i$
are the inverses of the diagonal entries of $L$. The denominator $d$ is set to
$\pm \mathrm{det}(S)$ where $S$ is an appropriate submatrix of $A$ ($S = A$ if
$A$ is square and nonsingular) and the sign is decided by the parity of the
permutation.
"""
function fflu(A::MatrixElem{T}, P = SymmetricGroup(nrows(A))) where {T <: RingElement}
   m = nrows(A)
   n = ncols(A)
   P.n != m && error("Permutation does not match matrix")
   p = one(P)
   R = base_ring(A)
   U = deepcopy(A)
   L = similar(A, m, m)
   rank, d = fflu!(p, U)
   i = 1
   j = 1
   k = 1
   while i <= m && j <= n
      if !is_zero_entry(U, i, j)
         L[i, k] = U[i, j]
         for l = 1:i - 1
            L[l, k] = 0
         end
         for l = i + 1:m
            L[l, k] = U[l, j]
            U[l, j] = 0
         end
         i += 1
         k += 1
      end
      j += 1
   end

   while k <= m
      for l = 1:k - 1
         L[l, k] = 0
      end
      L[k, k] = 1
      for l = k + 1: m
         L[l, k] = 0
      end
      k += 1
   end

   return rank, d, p, L, U
end

###############################################################################
#
#   Reduced row-echelon form
#
###############################################################################

function rref_rational!(A::MatrixElem{T}) where {T <: RingElement}
   m = nrows(A)
   n = ncols(A)
   R = base_ring(A)
   P = one(SymmetricGroup(m))
   rank, d = fflu!(P, A)
   for i = rank + 1:m
      for j = 1:n
         A[i, j] = R()
      end
   end
   if rank > 1
      t = R()
      q = R()
      d = -d
      pivots = zeros(Int, n)
      np = rank
      j = k = 1
      for i = 1:rank
         while is_zero_entry(A, i, j)
            pivots[np + k] = j
            j += 1
            k += 1
         end
         pivots[i] = j
         j += 1
      end
      while k <= n - rank
         pivots[np + k] = j
         j += 1
         k += 1
      end
      for k = 1:n - rank
         for i = rank - 1:-1:1
            t = mul_red!(t, A[i, pivots[np + k]], d, false)
            for j = i + 1:rank
               t = addmul_delayed_reduction!(t, A[i, pivots[j]], A[j, pivots[np + k]], q)
            end
            t = reduce!(t)
            A[i, pivots[np + k]] = divexact(-t, A[i, pivots[i]])
         end
      end
      d = -d
      for i = 1:rank
         for j = 1:rank
            if i == j
               A[j, pivots[i]] = d
            else
               A[j, pivots[i]] = R()
            end
         end
      end
   end
   return rank, d
end

@doc raw"""
    rref_rational(M::MatrixElem{T}) where {T <: RingElement}

Return a tuple $(r, A, d)$ consisting of the rank $r$ of $M$ and a
denominator $d$ in the base ring of $M$ and a matrix $A$ such that $A/d$ is
the reduced row echelon form of $M$. Note that the denominator is not usually
minimal.
"""
function rref_rational(M::MatrixElem{T}) where {T <: RingElement}
   A = deepcopy(M)
   r, d = rref_rational!(A)
   return r, A, d
end

function rref!(A::MatrixElem{T}) where {T <: FieldElement}
   m = nrows(A)::Int
   n = ncols(A)::Int
   R = base_ring(A)
   P = one(SymmetricGroup(m))
   rnk = lu!(P, A)
   if rnk == 0
      return 0
   end
   for i = 1:m
      for j = 1:min(rnk, i - 1)
         A[i, j] = R()
      end
   end
   U = zero_matrix(R, rnk, rnk)
   V = zero_matrix(R, rnk, n - rnk)
   pivots = zeros(Int, n)
   np = rnk
   j = k = 1
   for i = 1:rnk
      while is_zero_entry(A, i, j)
         pivots[np + k] = j
         j += 1
         k += 1
      end
      pivots[i] = j
      j += 1
   end
   while k <= n - rnk
      pivots[np + k] = j
      j += 1
      k += 1
   end
   for i = 1:rnk
      for j = 1:i
         U[j, i] = A[j, pivots[i]]
      end
      for j = i + 1:rnk
         U[j, i] = R()
      end
   end
   for i = 1:n - rnk
      for j = 1:rnk
         V[j, i] = A[j, pivots[np + i]]
      end
   end
   V = _solve_triu(U, V; unipotent = false, side = :right)
   for i = 1:rnk
      for j = 1:i
         A[j, pivots[i]] = i == j ? one(R) : R()
      end
   end
   for i = 1:n - rnk
      for j = 1:rnk
         A[j, pivots[np + i]] = V[j, i]
      end
   end
   return rnk
end

@doc raw"""
    rref(M::MatrixElem{T}) where {T <: FieldElement}

Return a tuple $(r, A)$ consisting of the rank $r$ of $M$ and a reduced row
echelon form $A$ of $M$.
"""
function rref(M::MatrixElem{T}) where {T <: FieldElement}
   A = deepcopy(M)
   r = rref!(A)
   return r, A
end

@doc raw"""
    is_rref(M::MatrixElem{T}) where {T <: RingElement}

Return `true` if $M$ is in reduced row echelon form, otherwise return
`false`.
"""
function is_rref(M::MatrixElem{T}) where {T <: RingElement}
   m = nrows(M)
   n = ncols(M)
   c = 1
   for r = 1:m
      for i = 1:c - 1
         if !is_zero_entry(M, r, i)
            return false
         end
      end
      while c <= n && is_zero_entry(M, r, c)
         c += 1
      end
      if c <= n
         for i = 1:r - 1
            if !is_zero_entry(M, i, c)
               return false
            end
         end
      end
   end
   return true
end

@doc raw"""
    is_rref(M::MatrixElem{T}) where {T <: FieldElement}

Return `true` if $M$ is in reduced row echelon form, otherwise return
`false`.
"""
function is_rref(M::MatrixElem{T}) where {T <: FieldElement}
   m = nrows(M)
   n = ncols(M)
   c = 1
   for r = 1:m
      for i = 1:c - 1
         if !is_zero_entry(M, r, i)
            return false
         end
      end
      while c <= n && is_zero_entry(M, r, c)
         c += 1
      end
      if c <= n
         if !isone(M[r, c])
            return false
         end
         for i = 1:r - 1
            if !is_zero_entry(M, i, c)
               return false
            end
         end
      end
   end
   return true
end

# Reduce row m of the matrix A, assuming the first m - 1 rows are in Gauss
# form. However those rows may not be in order. The i-th entry of the array
# P is the row of A which has a pivot in the i-th column. If no such row
# exists, the entry of P will be 0. The function returns the column in which
# the m-th row has a pivot after reduction. This will always be chosen to be
# the first available column for a pivot from the left. This information is
# also updated in P. The i-th entry of the array L contains the last column
# of A for which the row i is nonzero. This speeds up reduction in the case
# that A is chambered on the right. Otherwise the entries can all be set to
# the number of columns of A. The entries of L must be monotonic increasing.

function reduce_row!(A::MatrixElem{T}, P::Vector{Int}, L::Vector{Int}, m::Int) where {T <: FieldElement}
   R = base_ring(A)
   n = ncols(A)
   t = R()
   for i = 1:n
      # reduction of row was delayed, reduce next element now
      if i != 1
         A[m, i] = reduce!(A[m, i])
      end
      if !is_zero_entry(A, m, i)
         h = -A[m, i]
         r = P[i]
         if r != 0
            A[m, i] = R()
            for j = i + 1:L[r]
               t = mul_red!(t, A[r, j], h, false)
               A[m, j] = add!(A[m, j], t)
            end
         else
            # reduce remainder of row for return
            for j = i + 1:L[m]
               A[m, j] = reduce!(A[m, j])
            end
            h = inv(A[m, i])
            A[m, i] = one(R)
            for j = i + 1:L[m]
               A[m, j] = mul!(A[m, j], A[m, j], h)
            end
            P[i] = m
            return i
         end
      end
   end
   return 0
end

function reduce_row!(A::MatrixElem{T}, P::Vector{Int}, L::Vector{Int}, m::Int) where {T <: RingElement}
   R = base_ring(A)
   n = ncols(A)
   t = R()
   c = one(R)
   c1 = 0
   for i = 1:n
      if !is_zero_entry(A, m, i)
         h = -A[m, i]
         r = P[i]
         if r != 0
            d = A[r, i]
            A[m, i] = R()
            for j = i + 1:L[r]
               t = mul_red!(t, A[r, j], h, false)
               A[m, j] = mul_red!(A[m, j], A[m, j], d, false)
               A[m, j] = add!(A[m, j], t)
               A[m, j] = reduce!(A[m, j])
            end
            for j = L[r] + 1:L[m]
               A[m, j] = mul!(A[m, j], A[m, j], d)
            end
            if c1 > 0 && P[c1] < P[i]
               for j = i + 1:L[m]
                  A[m, j] = divexact(A[m, j], c)
               end
            end
            c1 = i
            c = d
         else
            P[i] = m
            return i
         end
      else
         r = P[i]
         if r != 0
            for j = i + 1:L[m]
               A[m, j] = mul!(A[m, j], A[m, j], A[r, i])
            end
         end
      end
   end
   return 0
end

###############################################################################
#
#   Determinant
#
###############################################################################

function det_clow(M::MatElem{T}) where {T <: RingElement}
   R = base_ring(M)
   n = nrows(M)
   if n == 0
      return one(R)
   end
   A = Matrix{T}(undef, n, n)
   B = Matrix{T}(undef, n, n)
   C = R()
   for i = 1:n
      for j = 1:n
         A[i, j] = i == j ? one(R) : zero(R)
         B[i, j] = R()
      end
   end
   for k = 1:n - 1
      for i = 1:n
         for j = 1:i
            if !iszero(A[i, j])
               for m = j + 1:n
                  C = mul!(C, A[i, j], M[i, m])
                  B[m, j] = add!(B[m, j], C)
               end
               for m = j + 1:n
                  C = mul!(C, A[i, j], M[i, j])
                  B[m, m] = add!(B[m, m], -C)
               end
            end
         end
      end
      Temp = A
      A = B
      B = Temp
      if k != n - 1
         for i = 1:n
            for j = 1:i
               B[i, j] = R()
            end
         end
      end
   end
   D = R()
   for i = 1:n
      for j = 1:i
         if !iszero(A[i, j])
            D -= A[i, j]*M[i, j]
         end
      end
   end
   return isodd(n) ? -D : D
end

function det_df(M::MatElem{T}) where {T <: RingElement}
   R = base_ring(M)
   S = poly_ring(R)
   n = nrows(M)
   p = charpoly(S, M)
   d = coeff(p, 0)
   return isodd(n) ? -d : d
end

function det_fflu(M::MatElem{T}) where {T <: RingElement}
   n = nrows(M)
   if n == 0
      return base_ring(M)()
   end
   A = deepcopy(M)
   P = one(SymmetricGroup(n))
   r, d = fflu!(P, A)
   return r < n ? base_ring(M)() : (parity(P) == 0 ? d : -d)
end

function det(M::MatElem{T}) where {T <: FieldElement}
   !is_square(M) && error("Not a square matrix in det")
   if nrows(M) == 0
      return one(base_ring(M))
   end
   return det_fflu(M)
end

@doc raw"""
    det(M::MatrixElem{T}) where {T <: RingElement}

Return the determinant of the matrix $M$. We assume $M$ is square.

# Examples

```jldoctest
julia> R, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over rationals, x)

julia> A = R[x 1; 1 x^2];

julia> d = det(A)
x^3 - 1
```
"""
function det(M::MatElem{T}) where {T <: RingElement}
   !is_square(M) && error("Not a square matrix in det")
   if nrows(M) == 0
      return one(base_ring(M))
   end
   try
      return det_fflu(M)
   catch
      return det_df(M)
   end
end

function det_interpolation(M::MatElem{T}) where {T <: PolyRingElem}
   n = nrows(M)
   !is_domain_type(elem_type(base_ring(base_ring(M)))) &&
          error("Generic interpolation requires a domain type")
   R = base_ring(M)
   if n == 0
      return R()
   end
   maxlen = 0
   for i = 1:n
      for j = 1:n
         maxlen = max(maxlen, length(M[i, j]))
      end
   end
   if maxlen == 0
      return R()
   end
   bound = n*(maxlen - 1) + 1
   x = Vector{elem_type(base_ring(R))}(undef, bound)
   d = Vector{elem_type(base_ring(R))}(undef, bound)
   X = zero_matrix(base_ring(R), n, n)
   b2 = div(bound, 2)
   pt1 = base_ring(R)(1 - b2)
   for i = 1:bound
      x[i] = base_ring(R)(i - b2)
      (x[i] == pt1 && i != 1) && error("Not enough interpolation points in ring")
      for j = 1:n
         for k = 1:n
            X[j, k] = evaluate(M[j, k], x[i])
         end
      end
      d[i] = det(X)
   end
   return interpolate(R, x, d)
end

function det(M::MatElem{T}) where {S <: FinFieldElem, T <: PolyRingElem{S}}
   !is_square(M) && error("Not a square matrix in det")
   if nrows(M) == 0
      return one(base_ring(M))
   end
   return det_popov(M)
end

function det(M::MatElem{T}) where {T <: PolyRingElem}
   !is_square(M) && error("Not a square matrix in det")
   if nrows(M) == 0
      return one(base_ring(M))
   end
   try
      return det_interpolation(M)
   catch
      # no point trying fflu, since it probably fails
      # for same reason as det_interpolation
      return det_df(M)
   end
end

###############################################################################
#
#   Minors
#
###############################################################################

@doc raw"""
    combinations(n::Int, k::Int)

Return an array consisting of k-combinations of {1,...,n} as arrays.
"""
combinations(n::Int, k::Int) = combinations(1:n, k)

@doc raw"""
    combinations(v::AbstractVector, k::Int)

Return an array consisting of k-combinations of a given vector v as arrays.
"""
function combinations(v::AbstractVector{T}, k::Int) where T
   n = length(v)
   ans = Vector{T}[]
   k > n && return ans
   _combinations_dfs!(ans, Vector{T}(undef, k), v, n, k)
   return ans
end
function _combinations_dfs!(ans::Vector{Vector{T}}, comb::Vector{T}, v::AbstractVector{T}, n::Int, k::Int) where T
   k < 1 && (pushfirst!(ans, comb[:]); return)
   for m in n:-1:k
      comb[k] = v[m]
      _combinations_dfs!(ans, comb, v, m - 1, k - 1)
   end
end

@doc raw"""
    minors(A::MatElem, k::Int)

Return an array consisting of the `k`-minors of `A`.

# Examples

```jldoctest
julia> A = ZZ[1 2 3; 4 5 6]
[1   2   3]
[4   5   6]

julia> minors(A, 2)
3-element Vector{BigInt}:
 -3
 -6
 -3

```
"""
minors(A::MatElem, k::Int) = collect(minors_iterator(A, k))

@doc raw"""
    minors_with_position(A::MatElem, k::Int)

Return an array consisting of the `k`-minors of `A` and the respective data on the rows and columns involved.

# Examples

```jldoctest
julia> A = ZZ[1 2 3; 4 5 6]
[1   2   3]
[4   5   6]

julia> minors_with_position(A, 2)
3-element Vector{Tuple{BigInt, Vector{Int64}, Vector{Int64}}}:
 (-3, [1, 2], [1, 2])
 (-6, [1, 2], [1, 3])
 (-3, [1, 2], [2, 3])

```
"""
minors_with_position(A::MatElem, k::Int) = collect(minors_iterator_with_position(A,k))

@doc raw"""
    minors_iterator(A::MatElem, k::Int)

Return an iterator that computes the `k`-minors of `A`.

# Examples

```jldoctest
julia> A = ZZ[1 2 3; 4 5 6]
[1   2   3]
[4   5   6]

julia> first(minors_iterator(A, 2))
-3

julia> collect(minors_iterator(A, 2))
3-element Vector{BigInt}:
 -3
 -6
 -3

```
"""
function minors_iterator(M::MatElem, k::Int)
  row_indices = combinations(nrows(M), k)
  col_indices = combinations(ncols(M), k)
  return (det(M[rows, cols]) for rows in row_indices for cols in col_indices)
end

@doc raw"""
    minors_iterator_with_position(A::MatElem, k::Int)

Return an iterator that computes the `k`-minors of `A` also specifying the row and column indices of the minor.

# Examples

```jldoctest
julia> A = ZZ[1 2 3; 4 5 6]
[1   2   3]
[4   5   6]

julia> first(minors_iterator_with_position(A, 2))
(-3, [1, 2], [1, 2])

```
"""
function minors_iterator_with_position(M::MatElem, k::Int)
  row_indices = combinations(nrows(M), k)
  col_indices = combinations(ncols(M), k)
  return ((det(M[rows, cols]), rows, cols) for rows in row_indices for cols in col_indices)
end

@doc raw"""
    exterior_power(A::MatElem, k::Int) -> MatElem

Return the `k`-th exterior power of `A`.

# Examples

```jldoctest
julia> A = matrix(ZZ, 3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9]);

julia> exterior_power(A, 2)
[-3    -6   -3]
[-6   -12   -6]
[-3    -6   -3]
```
"""
function exterior_power(A::MatElem, k::Int)
  ri = combinations(nrows(A), k)
  n = length(ri)
  res = similar(A, n, n)
   for i in 1:n
     for j in 1:n
       res[i, j] = det(A[ri[i], ri[j]])
     end
   end
   return res 
end

###############################################################################
#
#   Pfaffian
#
###############################################################################

"""
    is_alternating(M::MatrixElem)

Return whether the form corresponding to the matrix `M` is alternating,
i.e. `M = -transpose(M)` and `M` has zeros on the diagonal.
Return `false` if `M` is not a square matrix.
"""
function is_alternating(M::MatElem)
  is_skew_symmetric(M) || return false
  for i in 1:nrows(M)
    is_zero_entry(M, i, i) || return false
  end
  return true
end

"""
    is_skew_symmetric(M::MatrixElem)

Return `true` if the given matrix is skew symmetric with respect to its main
diagonal, i.e., `transpose(M) == -M`, otherwise return `false`.

# Examples
```jldoctest
julia> M = matrix(ZZ, [0 -1 -2; 1 0 -3; 2 3 0])
[0   -1   -2]
[1    0   -3]
[2    3    0]

julia> is_skew_symmetric(M)
true

```
"""
function is_skew_symmetric(M::MatElem)
   n = nrows(M)
   n == ncols(M) || return false
   for i in 1:n, j in 1:i
      M[i, j] == -M[j, i] || return false
   end
   return true
end

function check_skew_symmetric(M::MatElem)
   is_skew_symmetric(M) || throw(DomainError(M, "matrix must be skew-symmetric"))
   return M
end

@doc raw"""
    pfaffian(M::MatElem)

Return the Pfaffian of a skew-symmetric matrix `M`.
"""
function pfaffian(M::MatElem)
   check_skew_symmetric(M)
   # when the matrix is big, try use the BFL algorithm
   if ncols(M) > 10
      try
         return pfaffian_bfl_bsgs(M)
      catch
      end
   end
   # fallback to using recursion
   return pfaffian_r(M)
end

@doc raw"""
    pfaffians(M::MatElem, k::Int)

Return a vector consisting of the `k`-Pfaffians of a skew-symmetric matrix `M`.
"""
function pfaffians(M::MatElem, k::Int)
   check_skew_symmetric(M)
   indices = combinations(ncols(M), k)
   pfs = elem_type(base_ring(M))[]
   for i in indices
      push!(pfs, pfaffian(M[i, i]))
   end
   return pfs
end

# using recursion
pfaffian_r(M::MatElem) = _pfaffian(M, collect(1:ncols(M)), ncols(M))

function _pfaffian(M::MatElem, idx::Vector{Int}, k::Int)
   R = base_ring(M)
   k == 0 && return one(R)
   isodd(k) && return R()
   k == 2 && return M[idx[1], idx[2]]
   ans = R()
   sig = false
   for i in k - 1:-1:1
      idx[i], idx[k - 1] = idx[k - 1], idx[i]
      sig = !sig
      g = M[idx[k - 1], idx[k]]
      if !iszero(g)
         ans = (sig ? (+) : (-))(ans, g * _pfaffian(M, idx, k - 2))
      end
   end
   x = idx[k - 1]
   deleteat!(idx, k - 1)
   pushfirst!(idx, x) # restore idx
   return ans
end

# using the algorithm of Baer-Faddeev-LeVerrier
# the base ring of M should allow divisions of small integers
# (specifically, 2,4,6,...,n).
function pfaffian_bfl(M::MatElem)
   R = base_ring(M)
   n = ncols(M)
   characteristic(R) == 0 || characteristic(R) > n || throw(DomainError(M, "base ring must allow divisions of small integers"))
   n == 0 && return one(R)
   isodd(n) && return R()
   n == 2 && return M[1, 2]
   N = deepcopy(M)
   for i in 1:2:n
      for j in 1:n
         N[j, i], N[j, i + 1] = N[j, i + 1], -N[j, i]
      end
   end
   P = deepcopy(N)
   half_n = div(n, 2)
   for i in 1:half_n - 1
      P -= inv(R(2i)) * tr(P)
      P *= N
   end
   return (-1)^(half_n + 1) * inv(R(n)) * tr(P)
end

function trace_of_prod(M::MatElem, N::MatElem)
   is_square(M) && is_square(N) || error("Not a square matrix in trace")
   d = zero(base_ring(M))
   for i = 1:nrows(M)
      d += (M[i:i, :] * N[:, i:i])[1, 1]
   end
   return d
end

# use baby-step giant-step
# see https://arxiv.org/abs/2011.12573
function pfaffian_bfl_bsgs(M::MatElem)
   R = base_ring(M)
   n = ncols(M)
   characteristic(R) == 0 || characteristic(R) > n || throw(DomainError(M, "base ring must allow divisions of small integers"))
   n == 0 && return one(R)
   isodd(n) && return zero(R)
   n == 2 && return M[1, 2]
   N = deepcopy(M)
   for i in 1:2:n
      for j in 1:n
         N[j, i], N[j, i + 1] = N[j, i + 1], -N[j, i]
      end
   end

   # precompute the powers of N and their traces
   m = isqrt(n)
   N_power = [N]
   for i in 1:m - 1
      push!(N_power, N_power[end] * N)
   end
   t = tr.(N_power)

   P = identity_matrix(R, n)
   c = Vector{elem_type(R)}(undef, m)
   i = 1
   half_n = div(n, 2)
   while i <= half_n - 1
      m = min(m, half_n - i)
      # compute the coefficient c[m - j] before each N^j
      for j in 1:m
         # when i = 1, P = Id, so tr(N^j) is already known
         c[j] = (i == 1) ? t[j] : trace_of_prod(N_power[j], P)
         for k in 1:j - 1
            c[j] += t[k] * c[j - k]
         end
         c[j] *= -inv(R(2(i + j - 1)))
      end
      P *= N_power[m]
      for j in 1:m - 1
         P += c[m - j] * N_power[j]
      end
      P += c[m]
      i += m
   end

   return (-1)^(half_n + 1) * inv(R(n)) * trace_of_prod(N, P)
end

###############################################################################
#
#   Rank
#
###############################################################################

@doc raw"""
    rank(M::MatrixElem{T}) where {T <: RingElement}

Return the rank of the matrix $M$.

# Examples

```jldoctest
julia> A = QQ[1 2; 3 4];

julia> d = rank(A)
2
```
"""
function rank(M::MatElem{T}) where {T <: RingElement}
   n = nrows(M)
   if n == 0
      return 0
   end
   A = deepcopy(M)
   P = one(SymmetricGroup(n))
   r, d = fflu!(P, A)
   return r
end

function rank(M::MatElem{T}) where {T <: FieldElement}
   n = nrows(M)
   if n == 0
      return 0
   end
   A = deepcopy(M)
   P = one(SymmetricGroup(n))
   return lu!(P, A)
end

@doc raw"""
    rank_interpolation(M::MatrixElem{T}) where {T <: PolyRingElem} -> Int
    rank_interpolation(M::MatrixElem{T}) where {T <: MPolyRingElem} -> Int
    rank_interpolation(M::MatrixElem{T}) where {T <: AbstractAlgebra.Generic.RationalFunctionFieldElem} -> Int

Returns the rank of $A$ using an interpolation-like method. 

# Examples

```jldoctest
julia> Qx, x = polynomial_ring(QQ, :x);

julia> AbstractAlgebra.rank_interpolation(matrix(Qx, 2, 2, [x 2x; x^2 1]))
2

julia> Qy, y = rational_function_field(QQ, :y);

julia> AbstractAlgebra.rank_interpolation(matrix(Qy, 2, 2, [1//y -y^2; 3 2-y]))
2
```
"""
function rank_interpolation(M::MatElem{T}) where {T <: PolyRingElem}
   n = nrows(M)
   m = ncols(M)
   if is_zero(n) || is_zero(m)
      return 0
   end
   Kx = base_ring(M)
   K = base_ring(Kx)
   #The maximum degree of det(M') is calculated where M' is an arbitrary quadratic submatrix of M.
   min_ = min(n, m)
   if (min_ == n)
      maxdetdeg = sum(maximum(degree(M[i, j]) for j in 1:m) for i in 1:n)
   else 
      maxdetdeg = sum(maximum(degree(M[i, j]) for i in 1:n) for j in 1:m)
   end
   r = 0
   eval_set = evaluation_points(K, maxdetdeg+1)
   if !is_empty(eval_set)
      #rank(M) is calculated by computing the rank of det_deg+1 matrices, evaluated in an element of eval_set, respectively.
      for elem in eval_set
         M_eval = map_entries(p -> evaluate(p, elem), M)
         r = max(rank_interpolation(M_eval), r)
         if r == min_
            break
         end
      end
      return r
   else
      #function get_eval_set returns an empty set if and only if K is a finite field and order(K) < det_deg+1 holds. 
      #In this case a field extension L of K such that order(L) >= det_deg+1 is constructed.
      #d is the smallest natural number such that order(K)^d >= det_deg+1.
      d = Int(ceil(Base.log(BigInt(order(K)), maxdetdeg + 1, )))
      @assert order(K)^d >= maxdetdeg + 1
      #d = clog(order(K), ZZ(maxdetdeg+1))
      L, l = ext_of_degree(K, d)
      Lx, _ = polynomial_ring(L, var(Kx)) 
      #The given matrix M is embedded into the space of matrices over Lx
      A = matrix(Lx, n, m, [map_coefficients(l, M[i, j]; parent = Lx) for i in 1:n, j in 1:m])
      return rank_interpolation(A)
   end
end

function rank_interpolation(M::MatElem{T}) where {T <: MPolyRingElem}
   n = nrows(M)
   m = ncols(M)
   if is_zero(n) || is_zero(m)
      return 0
   end
   Kx = base_ring(M)
   K = base_ring(Kx)
   num_vars = number_of_variables(Kx)
   #The maximum degree of det(M') is calculated where M' is an arbitrary quadratic submatrix of M.
   min_ = min(n, m)
   if min_ == n
      maxdetdeg = [sum(maximum(degree(M[i, j], k) for j in 1:m) for i in 1:n) for k in 1:num_vars]
   else
      maxdetdeg = [sum(maximum(degree(M[i, j], k) for i in 1:n) for j in 1:m) for k in 1:num_vars] 
   end
   r = 0
   eval_set = Vector{Vector{elem_type(K)}}(undef, num_vars)
   for i = 1:num_vars
      eval_set[i] = evaluation_points(K, maxdetdeg[i]+1)
      if is_empty(eval_set[i])
         #function get_eval_set returns an empty set if and only if K is a finite field and order(K) < det_deg+1 holds. 
         #In this case a field extension L of K such that order(L) >= det_deg+1 is constructed.
         #d is the smallest natural number such that order(K)^d >= det_deg+1.
         d = Int(ceil(Base.log(BigInt(order(K)), maximum(maxdetdeg)+1)))
         @assert order(K)^d >= maximum(maxdetdeg) + 1
         #d = clog(order(K), ZZ(maximum(maxdetdeg)+1))
         L, l = ext_of_degree(K, d)
         Lx, _ = polynomial_ring(L, symbols(Kx)) 
         #The given matrix M is embedded into the space of matrices over Lx
         A = matrix(Lx, n, m, [map_coefficients(l, M[i, j]; parent = Lx) for i in 1:n, j in 1:m])
         return rank_interpolation(A)
      end
   end
   #rank(M) is calculated by computing the rank of (det_deg+1)^number_of_variables(Kx) matrices, evaluated in elements of eval_set.
   for tup in ProductIterator(eval_set)
      M_eval = map_entries(p -> evaluate(p, tup), M)
      r = max(rank_interpolation(M_eval), r)
      if r == min_
         break
      end
   end
   return r
end

function rank_interpolation(M::MatElem{T}) where {T <: RingElement}
   return rank(M)
end

@doc raw"""
    rank_interpolation_mc(M::MatrixElem{T}, ϵ::Float64) where {T <: PolyRingElem} -> Int
    rank_interpolation_mc(M::MatrixElem{T}, ϵ::Float64) where {T <: MPolyRingElem} -> Int
    rank_interpolation_mc(M::MatrixElem{T}, ϵ::Float64) where {T <: AbstractAlgebra.Generic.RationalFunctionFieldElem} -> Int

Returns the rank of $A$ with error probability < $ϵ$ using an interpolation-like method.

# Examples

```jldoctest
julia> Qx, x = polynomial_ring(QQ, :x);

julia> AbstractAlgebra.rank_interpolation_mc(matrix(Qx, 2, 2, [x 2x; x^2 1]), 0.01)
2

julia> Qy, y = rational_function_field(QQ, :y);

julia> AbstractAlgebra.rank_interpolation_mc(matrix(Qy, 2, 2, [1//y -y^2; 3 2-y]), 0.00001)
2
```
"""
function rank_interpolation_mc(M::MatElem{T}, err::Float64) where {T <: PolyRingElem}
   n = nrows(M)
   m = ncols(M)
   if is_zero(n) || is_zero(m)
      return 0
   end
   Kx = base_ring(M)
   K = base_ring(Kx)
   min_ = min(n, m)
   #The maximum degree of det(M') is calculated where M' is an arbitrary quadratic submatrix of M.
   maxdetdeg = min_*maximum(degree(M[i, j]) for i in 1:n, j in 1:n)
   S = evaluation_points(K, 10*maxdetdeg)
   #k is the minimum amount of evaluations of M needed to compute the correct rank of M with error probability < err
   k = ceil(Base.log(10, 1/err))
   r = 0
   if !is_empty(S)
      #rank(M) is calculated by computing the rank of k matrices, evaluated in elements of eval_set, and taking the maximum
      for _ = 1:k
         a = rand(S)
         M_eval = map_entries(p -> evaluate(p, a), M)
         r = max(rank(M_eval), r)
         if r == min_
            return r
         end
      end
   else
      #function get_eval_set returns an empty set if and only if K is a finite field and order(K) < det_deg+1 holds. 
      #In this case a field extension L of K such that order(L) >= det_deg+1 is constructed.
      #d is the smallest natural number such that order(K)^d >= det_deg+1.
      d = Int(ceil(Base.log(BigInt(order(K)), maxdetdeg*10)))
      #d = clog(order(K), ZZ(maxdetdeg*10))
      @assert order(K)^d >= maxdetdeg*10
      L, l = ext_of_degree(K, d)
      Lx, _ = polynomial_ring(L, var(Kx)) 
      #The given matrix M is embedded into the space of matrices over Lx
      A = matrix(Lx, n, m, [map_coefficients(l, M[i, j]; parent = Lx) for i in 1:n, j in 1:m])
      return rank_interpolation_mc(A, err)
   end
   return r
end

function rank_interpolation_mc(M::MatElem{T}, err::Float64) where {T <: MPolyRingElem}
   n = nrows(M)
   m = ncols(M)
   if is_zero(n) || is_zero(m)
      return 0
   end
   Kx = base_ring(M)
   K = base_ring(Kx)
   num_vars = number_of_variables(Kx)
   min_ = min(n, m)
   #The maximum degree of det(M') is calculated where M' is an arbitrary quadratic submatrix of M.
   maxdetdeg = min_*maximum(degree(M[i, j], k) for i in 1:n, j in 1:m, k in 1:num_vars)
   S = evaluation_points(K, 10*maxdetdeg)
   #k is the minimum amount of evaluations of M needed to compute the correct rank of M with error probability < err
   k = ceil(Base.log(10, 1/err))
   r = 0
   if !is_empty(S)
      #rank(M) is calculated by computing the rank of k matrices, evaluated in elements of eval_set, and taking the maximum
      for _ = 1:k
         a = Vector{elem_type(K)}(undef, num_vars)
         for i = 1:num_vars
            a[i] = rand(S)
         end
         M_eval = map_entries(p -> evaluate(p, a), M)
         r = max(rank(M_eval), r)
         if r == min_
            return r
         end
      end
      return r
   else
      #function get_eval_set returns an empty set if and only if K is a finite field and order(K) < det_deg*s holds. 
      #In this case a field extension L of K is constructed such that order(L) >= det_deg*s.
      #d is the smallest natural number such that order(K)^d >= det_deg*s.
      d = Int(ceil(Base.log(BigInt(order(K)), maxdetdeg*10)))
      #d = clog(order(K), ZZ(maxdetdeg*10))
      @assert order(K)^d >= maxdetdeg*10
      L, l = ext_of_degree(K, d)
      Lx, _ = polynomial_ring(L, symbols(Kx)) 
      #The given matrix M is embedded into the space of matrices over Lx
      A = matrix(Lx, n, m, [map_coefficients(l, M[i, j]; parent = Lx) for i in 1:n, j in 1:m])
      return rank_interpolation_mc(A, err)
   end
end



###############################################################################
#
#   Linear solving
#
###############################################################################

# Return `flag, y, d` where flag is set to `true` if `Ay = bd` has a solution
# in the base ring. If not, it returns false. The matrix A can be non-square
# and singular. If a solution exists, `y` is set to one such solution. The
# value `d` is set to an appropriate denominator so that `Ay = bd` is a
# solution over the ring. This implements the fraction free decomposition of
# David J. Jeffrey, see "LU Factoring of non-invertible matrices" in ACM
# Communications in Computer Algebra, July 2010. Note that we handle column
# permutations implicitly and add units along the diagonal of the lower
# triangular matrix L instead of removing columns (and corresponding rows of
# the upper triangular matrix U). We also set free variables to zero.
function _can_solve_with_solution_fflu(A::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   base_ring(A) != base_ring(b) && error("Base rings don't match in _can_solve_with_solution_fflu")
   nrows(A) != nrows(b) && error("Dimensions don't match in _can_solve_with_solution_fflu")
   FFLU = deepcopy(A)
   p = one(SymmetricGroup(nrows(A)))
   rank, d = fflu!(p, FFLU)
   flag, y = _solve_fflu_precomp(p, FFLU, b)
   n = nrows(A)
   if flag && rank < n
      b2 = p*b
      A2 = p*A
      A3 = A2[rank + 1:n, :]
      flag = A3*y == b2[rank + 1:n, :]*d
   end
   return flag, y, d
end

# Given a fraction free LU decomposition `LdU` of `p(A)` over an integral
# domain with `L` invertible over fraction field, this function will return a
# pair `flag, y` such that `Ay = bd` and `flag` is set to `true` if a solution
# exists. If not, `flag` may be set to `false`, however it is not required to
# be. If `r` is the rank of `A` then the first `r` rows of `p(A)y = p(b)d` will
# hold iff `flag` is `true`. The remaining rows must be checked by the caller.
function _solve_fflu_precomp(p::Perm, FFLU::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   x = p * b
   n = nrows(x)
   m = ncols(x)
   R = base_ring(FFLU)
   c = ncols(FFLU)
   t = base_ring(b)()
   s = base_ring(b)()
   minus_one = R(-1)
   y = similar(x, c, m)
   diag = Vector{elem_type(R)}(undef, n)
   piv = Vector{Int}(undef, n)

   rnk = 0
   l = 0
   for i = 1:n
      l += 1
      while l <= c && is_zero_entry(FFLU, i, l)
         l += 1
      end
      piv[i] = l
      if l <= c
         diag[i] = FFLU[i, l]
         rnk += 1
      end
   end

   for k in 1:m
      for i in 1:(n - 1)
         t = mul!(t, x[i, k], minus_one)
         for j in (i + 1):n
            if i <= c && piv[i] <= c
               if i == 1
                  x[j, k] = mul_red!(R(), x[j, k], diag[i], false)
               else
                  x[j, k] = mul_red!(x[j, k], x[j, k], diag[i], false)
               end
            else
               x[j, k] = deepcopy(x[j, k])
            end
            if i <= c && piv[i] <= c
               s = mul_red!(s, FFLU[j, piv[i]], t, false)
               x[j, k] = add!(x[j, k], s)
            end
            x[j, k] = reduce!(x[j, k])
            if i > 1
               if i <= c && piv[i - 1] <= c
                  flag, x[j, k] = divides(x[j, k], diag[i - 1])
                  if !flag
                     return false, x
                  end
               end
            end
         end
      end

      l = rnk
      for i in c:-1:1
         if l > 0 && i == piv[l]
            if rnk != 0
               y[i, k] = x[l, k]*diag[rnk]
            else
               y[i, k] = deepcopy(x[l, k])
            end
            for j in (piv[l] + 1):c
               t = mul!(t, y[j, k], FFLU[l, j])
               t = mul!(t, t, minus_one)
               y[i, k] = add!(y[i, k], t)
            end

            flag, y[i, k] = divides(y[i, k], diag[l])
            if !flag
               return false, x
            end

            l -= 1
         else
            y[i, k] = R()
         end
      end
   end

   return true, y
end

# Return `flag, y` where flag is set to `true` if `Ay = b` has a solution
# in the base field. If not, it returns false. The matrix A can be non-square
# and singular. If a solution exists, `y` is set to one such solution.
# This implements the LU decomposition, for non-invertible matrices, of
# David J. Jeffrey, see "LU Factoring of non-invertible matrices" in ACM
# Communications in Computer Algebra, July 2010. Note that we handle column
# permutations implicitly and add units along the diagonal of the lower
# triangular matrix L instead of removing columns (and corresponding rows of
# the upper triangular matrix U). We also set free variables to zero.
function _can_solve_with_solution_lu(A::MatElem{T}, b::MatElem{T}) where {T <: FieldElement}
   base_ring(A) != base_ring(b) && error("Base rings don't match in can_solve_with_solution_lu")
   nrows(A) != nrows(b) && error("Dimensions don't match in can_solve_with_solution_lu")

   if nrows(A) == 0
      return true, zero_matrix(base_ring(A), ncols(A), ncols(b))
   end

   if ncols(A) == 0
      return iszero(b), zero_matrix(base_ring(A), ncols(A), ncols(b))
   end

   LU = deepcopy(A)
   p = one(SymmetricGroup(nrows(A)))
   rank = lu!(p, LU)

   y = _solve_lu_precomp(p, LU, b)

   n = nrows(A)
   flag = true
   if rank < n
      b2 = p*b
      A2 = p*A
      A3 = A2[rank + 1:n, :]
      flag = A3*y == b2[rank + 1:n, :]
   end

   return flag, y
end

# Given an LU decomposition `LU` of `p(A)` over a field with `L` invertible,
# this function will return `y` such that the first `r` rows of `p(A)y = p(b)`
# hold, where `r` is the rank of `A`. The remaining rows must be checked by
# the caller.
function _solve_lu_precomp(p::Perm, LU::MatElem{T}, b::MatElem{T}) where {T <: FieldElement}
   x = p * b
   n = nrows(x)
   m = ncols(x)
   R = base_ring(LU)
   c = ncols(LU)
   t = base_ring(b)()
   s = base_ring(b)()
   y = similar(x, c, m)

   diag = Vector{elem_type(R)}(undef, n)
   piv = Vector{Int}(undef, n)

   l = 0
   rnk = 0
   for i = 1:n
      l += 1
      while l <= c && is_zero_entry(LU, i, l)
         l += 1
      end
      piv[i] = l
      if l <= c
         diag[i] = LU[i, l]
         rnk += 1
      end
   end

   for k in 1:m
      x[1, k] = deepcopy(x[1, k])
      for i in 2:n
         for j in 1:(i - 1)
            # x[i, k] = x[i, k] - LU[i, j] * x[j, k]
            if j <= c
               t = mul_red!(t, -LU[i, j], x[j, k], false)
               if j == 1
                  x[i, k] = x[i, k] + t # LU[i, j] * x[j, k]
               else
                  x[i, k] = add!(x[i, k], t)
               end
            else
               x[i, k] = deepcopy(x[i, k])
            end
         end
         x[i, k] = reduce!(x[i, k])
      end

      # Now every entry of x is a proper copy, so we can change the entries
      # as much as we want.

      l = rnk
      for i in c:-1:1
         if l > 0 && i == piv[l]
            y[i, k] = x[l, k]

            for j in (piv[l] + 1):c
               # x[i, k] = x[i, k] - x[j, k] * LU[l, j]
               t = mul_red!(t, y[j, k], -LU[l, j], false)
               y[i, k] = add!(y[i, k], t)
            end

            y[i, k] = reduce!(y[i, k])
            y[i, k] = divexact(y[i, k], diag[l])

            l -= 1
         else
            y[i, k] = R()
         end
      end
   end

   return y
end

function _solve_ff(M::MatrixElem{T}, b::MatrixElem{T}) where {T <: FieldElement}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   m = nrows(M)
   flag, x, d = _can_solve_with_solution_fflu(M, b)
   !flag && error("System not solvable in _solve_ff")
   for i in 1:nrows(x)
      for j in 1:ncols(x)
         x[i, j] = divexact(x[i, j], d)
      end
   end
   return x
end

function _can_solve_with_solution_with_det(M::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   # We cannot use _solve_fflu directly, since it forgot about the (parity of
   # the) permutation.
   R = base_ring(M)
   FFLU = deepcopy(M)
   p = one(SymmetricGroup(nrows(M)))
   rank, d = fflu!(p, FFLU)
   pivots = zeros(Int, nrows(M))
   c = 1
   for r = 1:nrows(M)
      while c <= ncols(M)
         if !is_zero_entry(FFLU, r, c)
            pivots[r] = c
            c += 1
            break
         end
         c += 1
      end
   end
   flag, x = _solve_fflu_precomp(p, FFLU, b)
   n = nrows(M)
   if flag && rank < n
      b2 = p*b
      A2 = p*M
      A3 = A2[rank + 1:n, :]
      flag = A3*x == b2[rank + 1:n, :]*d
   end
   if !flag
      return false, rank, p, pivots, x, d
   end
   # Now M*x = d*b, but d is only sign(P) * det(M)
   if parity(p) != 0
      minus_one = R(-1)
      for k in 1:ncols(x)
         for i in 1:nrows(x)
            # We are allowed to modify x in-place.
            x[i, k] = mul!(x[i, k], x[i, k], minus_one)
         end
      end
      d = mul!(d, d, minus_one)
   end
   return true, rank, p, pivots, x, d
end

function _can_solve_with_solution_with_det(M::MatElem{T}, b::MatElem{T}) where {T <: PolyRingElem}
   flag, r, p, pivots, x, d = can_solve_with_solution_interpolation_inner(M, b)
   return flag, r, p, pivots, x, d
end

# This can be removed once Nemo implements _can_solve_with_solution_with_det
# It's here now only because Nemo overloads it
function _solve_with_det(M::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   flag, r, p, piv, x, d = _can_solve_with_solution_with_det(M, b)
   !flag && error("System not solvable in _solve_with_det")
   return x, d
end

function _solve_ff(M::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   m = nrows(M)
   n = ncols(M)
   if m == 0
      return zero_matrix(base_ring(M), ncols(M), ncols(b)), base_ring(M)()
   end
   if n == 0
      b != 0 && error("System not soluble in _solve_ff")
      return zero_matrix(base_ring(M), ncols(M), ncols(b)), base_ring(M)()
   end
   flag, S, d = _can_solve_with_solution_fflu(M, b)
   !flag && error("System not soluble in _solve_ff")
   return S, d
end

function can_solve_with_solution_interpolation_inner(M::MatElem{T}, b::MatElem{T}) where {T <: PolyRingElem}
   m = nrows(M)
   h = ncols(b)
   c = ncols(M)
   R = base_ring(M)
   prm = one(SymmetricGroup(nrows(M)))
   pivots = zeros(Int, nrows(M))
   if m == 0
      return true, 0, prm, pivots, zero_matrix(R, c, h), one(R)
   end
   R = base_ring(M)
   maxlen = 0
   for i = 1:m
      for j = 1:c
         maxlen = max(maxlen, length(M[i, j]))
      end
   end
   if maxlen == 0
      return iszero(b), 0, prm, pivots, zero_matrix(R, c, h), one(R)
   end
   maxlenb = 0
   for i = 1:m
      for j = 1:h
         maxlenb = max(maxlenb, length(b[i, j]))
      end
   end
   # bound from xd = (M*)b where d is the det
   bound = (maxlen - 1)*(max(m, c) - 1) + max(maxlenb, maxlen)
   tmat = matrix(base_ring(R), 0, 0, elem_type(base_ring(R))[])
   V = Vector{typeof(tmat)}(undef, bound)
   d = Vector{elem_type(base_ring(R))}(undef, bound)
   y = Vector{elem_type(base_ring(R))}(undef, bound)
   bj = Vector{elem_type(base_ring(R))}(undef, bound)
   X = similar(tmat, m, c)
   Y = similar(tmat, m, h)
   x = similar(b, c, h)
   b2 = div(bound, 2)
   pt1 = base_ring(R)(1 - b2)
   l = 1
   i = 1
   pt = 1
   rnk = -1
   firstprm = true
   while l <= bound
      y[l] = base_ring(R)(pt - b2)
      # Running out of interpolation points doesn't imply there is no solution
      (y[l] == pt1 && pt != 1) && error("Not enough interpolation points in ring")
      bad_evaluation = false
      for j = 1:m
         for k = 1:c
            X[j, k] = evaluate(M[j, k], y[l])
            if is_zero_entry(X, j, k) && !is_zero_entry(M, j, k)
               bad_evaluation = true
               break
            end
         end
         if bad_evaluation
            break
         end
         for k = 1:h
            Y[j, k] = evaluate(b[j, k], y[l])
            if is_zero_entry(Y, j, k) && !is_zero_entry(b, j, k)
               bad_evaluation = true
               break
            end
         end
         if bad_evaluation
            break
         end
      end
      try
         if bad_evaluation
            error("Bad evaluation point")
         end
         flag, r, p, pv, Vl, dl = _can_solve_with_solution_with_det(X, Y)
         if !flag
            return flag, r, p, pv, zero(x), zero(R)
         end
         p = inv!(p)
         # Check that new solution has the same pivots as previous ones
         if r != rnk || p != prm || pv != pivots
            if r < rnk # rank is too low: reject
               pt += 1
               continue
            elseif r > rnk # rank has increased: restart
               l = 1
               rnk = r
               prm = p
               pivots = pv
               firstprm = false
            elseif p != prm || pv != pivots # pivot structure different
               reset_prm = false
               for j = 1:length(p.d)
                  # If earlier pivots or row swaps are encountered, restart
                  if firstprm || (p[j] < prm[j] && pv[j] <= pivots[j]) ||
                                 (p[j] == prm[j] && pv[j] < pivots[j] && pv[j] != 0)
                     prm = p
                     pivots = pv
                     l = 1
                     firstprm = false
                     break
                  elseif p[j] > prm[j] || pv[j] > pivots[j] # worse pivots/row swaps: reject
                     reset_prm = true
                     break
                  end
               end
               if reset_prm
                  pt += 1
                  continue
               end
            end
         end
         V[l] = Vl
         d[l] = dl
         y[l] = base_ring(R)(pt - b2)
         l += 1
      catch e
         if !(e isa ErrorException)
            rethrow(e)
         end
         i = i + 1
      end

      # We tested bound evaluation points and an impossible inverse was
      # encountered for all the values.

      if i > bound && l == 1
         # impossible inverse doesn't imply no solution
         error("Impossible inverse or too many failures in can_solve_with_solution_interpolation")
      end

      pt = pt + 1
   end
   # Interpolate
   for k = 1:h
      for i = 1:c
         for j = 1:bound
            bj[j] = V[j][i, k]
         end
         try
            x[i, k] = interpolate(R, y, bj)
         catch e
            if !(e isa ErrorException)
               rethrow(e)
            end
            return false, rnk, prm, pivots, zero(x), zero(R)
         end
      end
   end
   return true, rnk, prm, pivots, x, interpolate(R, y, d)
end

function _can_solve_with_solution_interpolation(M::MatElem{T}, b::MatElem{T}) where {T <: PolyRingElem}
   flag, r, p, pv, x, d = can_solve_with_solution_interpolation_inner(M, b)
   return flag, x, d
end

@doc raw"""
    _solve_rational(M::MatElem{T}, b::MatElem{T}) where T <: RingElement

Given a non-singular $n\times n$ matrix over a ring and an $n\times m$
matrix over the same ring, return a tuple $x, d$ consisting of an
$n\times m$ matrix $x$ and a denominator $d$ such that $Ax = db$. The
denominator will be the determinant of $A$ up to sign. If $A$ is singular an
exception is raised.
"""
function _solve_rational(M::MatElem{T}, b::MatElem{T}) where T <: RingElement
   return _solve_ringelem(M, b)
end

function _solve_ringelem(M::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   return _solve_ff(M, b)
end

function _solve_rational(M::MatElem{T}, b::MatElem{T}) where {T <: PolyRingElem}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   flag = true
   try
      flag, x, d = _can_solve_with_solution_interpolation(M, b)
      !flag && error("No solution in _solve_rational")
      return x, d
   catch e
      if !isa(e, ErrorException)
         rethrow(e)
      end
      !flag && error("No solution in _solve_rational")
      return _solve_ff(M, b)
   end
end

# Find the pivot columns of an rref matrix
function find_pivot(A::MatElem{T}) where T <: RingElement
  p = Int[]
  j = 0
  for i = 1:nrows(A)
    j += 1
    if j > ncols(A)
      return p
    end
    while is_zero_entry(A, i, j)
      j += 1
      if j > ncols(A)
        return p
      end
    end
    push!(p, j)
  end
  return p
end

###############################################################################
#
#   Upper triangular solving
#
###############################################################################

@doc raw"""
    is_upper_triangular(A::MatrixElem)

Return `true` if $A$ is an upper triangular matrix, that is,
all entries below the main diagonal are zero. Note that this
definition also applies to non-square matrices.

Alias for `LinearAlgebra.istriu`.

# Examples
```jldoctest
julia> is_upper_triangular(QQ[1 2 ; 0 4])
true

julia> is_upper_triangular(QQ[1 0 ; 3 4])
false

julia> is_upper_triangular(QQ[1 2 ;])
true

julia> is_upper_triangular(QQ[1 ; 2])
false
```
"""
function is_upper_triangular(M::MatElem)
    m = ncols(M)
    for i = 2:nrows(M)
        for j = 1:min(i - 1, m)
            if !is_zero_entry(M, i, j)
                return false
            end
        end
    end
    return true
end

# Remove the following once Nemo is adjusted
_solve_triu_right(M, b; unipotent = false) = _solve_triu(M, b; unipotent, side = :right)

@doc raw"""
    _solve_triu(U::MatElem{T}, b::MatElem{T}; unipotent::Bool = false, side::Symbol = :left) where {T <: RingElement}

Let $U$ be a non-singular $n\times n$ upper triangular matrix $U$ over a field. If 
`side = :right`, let $b$ 
be an $n\times m$ matrix $b$ over the same field, return an
$n\times m$ matrix $x$ such that $Ux = b$. If this is not possible, an error
will be raised.

If `side = :left`, the default, $b$ has to be $m \times n$. In this case
$xU = b$ is solved - or an error raised.

See also [`AbstractAlgebra._solve_triu_left`](@ref) and [`Strassen`](@ref) for
  asymptotically fast versions.
"""
function _solve_triu(U::MatElem{T}, b::MatElem{T}; unipotent::Bool = false, side::Symbol = :left) where {T <: RingElement}
   if side == :left
     return _solve_triu_left(U, b; unipotent)
   end
   @assert side == :right
   n = nrows(U)
   m = ncols(b)
   R = base_ring(U)
   X = zero(b)
   tmp = Vector{elem_type(R)}(undef, n)
   t = R()
   for i = 1:m
      for j = 1:n
         tmp[j] = X[j, i]
      end
      for j = n:-1:1
         s = R(0)
         for k = j + 1:n
            s = addmul!(s, U[j, k], tmp[k], t)
#            s = s + U[j, k] * tmp[k]
         end
         s = sub!(s, b[j, i], s)
         if unipotent
           tmp[j] = s
         else
           tmp[j] = divexact(s, U[j,j])
         end
      end
      for j = 1:n
         X[j, i] = tmp[j]
      end
   end
   return X
end

@doc raw"""
    _solve_triu_left(U::MatElem{T}, b::MatElem{T}; unipotent::Bool = false) where {T <: RingElement}

Given a non-singular $n\times n$ matrix $U$ over a field which is upper
triangular, and an $n\times m$ matrix $b$ over the same ring, return an
$n\times m$ matrix $x$ such that $Ux = b$. If this is not possible, an error
will be raised.

See also [`_solve_triu`](@ref) and [`Strassen`](@ref) for asymptotically fast 
  versions.
"""
function _solve_triu_left(U::MatElem{T}, b::MatElem{T}; unipotent::Bool = false) where {T <: RingElement}
   n = ncols(U)
   m = nrows(b)
   R = base_ring(U)
   X = zero(b)
   tmp = Vector{elem_type(R)}(undef, n)
   t = R()
   for i = 1:m
      for j = 1:n
         tmp[j] = X[i, j]
      end
      for j = 1:n
         s = R()
         for k = 1:j-1
            s = addmul!(s, U[k, j], tmp[k], t)
         end
         s = sub!(s, b[i, j], s)
         if unipotent
           tmp[j] = s
         else
           tmp[j] = divexact(s, U[j,j])
         end
      end
      for j = 1:n
         X[i, j] = tmp[j]
      end
   end
   return X
end

#solves A x = B for A intended to be lower triangular
#only the lower part is used. if f is true, then the diagonal is assumed to be 1
#used to use lu!
#can be combined with Strassen._solve_tril!
function _solve_tril!(A::MatElem{T}, B::MatElem{T}, C::MatElem{T}, f::Int = 0) where T

  # a       x   u      ax = u
  # b c   * y = v      bx + cy = v
  # d e f   z   w      ....

  @assert ncols(A) == ncols(C)
  s = base_ring(A)(0)
  for i=1:ncols(A)
    for j = 1:nrows(A)
      t = C[j, i]
      for k = 1:j-1
        mul_red!(s, A[k, i], B[j, k], false)
        t = sub!(t, s)
      end
      reduce!(t)
      if f == 1
        A[j,i] = t
      else
        A[j,i] = divexact(t, B[j, j])
      end
    end
  end
end

###############################################################################
#
#
#
###############################################################################

@doc raw"""
    is_lower_triangular(A::MatrixElem)

Return `true` if $A$ is an lower triangular matrix, that is,
all entries above the main diagonal are zero. Note that this
definition also applies to non-square matrices.

Alias for `LinearAlgebra.istril`.

# Examples
```jldoctest
julia> is_lower_triangular(QQ[1 2 ; 0 4])
false

julia> is_lower_triangular(QQ[1 0 ; 3 4])
true

julia> is_lower_triangular(QQ[1 2 ;])
false

julia> is_lower_triangular(QQ[1 ; 2])
true
```
"""
function is_lower_triangular(M::MatElem)
    for i = 1:nrows(M)
        for j = i+1:ncols(M)
            if !is_zero_entry(M, i, j)
                return false
            end
        end
    end
    return true
end

@doc raw"""
    is_diagonal(A::MatrixElem)

Return `true` if $A$ is a diagonal matrix, that is,
all entries off the main diagonal are zero. Note that this
definition also applies to non-square matrices.

Alias for `LinearAlgebra.isdiag`.

# Examples
```jldoctest
julia> is_diagonal(QQ[1 0 ; 0 4])
true

julia> is_diagonal(QQ[1 2 ; 3 4])
false

julia> is_diagonal(QQ[1 0 ;])
true
```
"""
function is_diagonal(A::MatElem)
    for i = 1:ncols(A)
        for j = 1:nrows(A)
            if i != j && !is_zero_entry(A, j, i)
                return false
            end
        end
    end
    return true
end

###############################################################################
#
#   Inverse
#
###############################################################################

@doc raw"""
    pseudo_inv(M::MatrixElem{T}) where {T <: RingElement}

Given a non-singular $n\times n$ matrix $M$ over a ring return a tuple $X, d$
consisting of an $n\times n$ matrix $X$ and a denominator $d$ such that
$MX = dI_n$, where $I_n$ is the $n\times n$ identity matrix. The denominator
will be the determinant of $M$ up to sign. If $M$ is singular an exception
is raised.
"""
function pseudo_inv(M::MatrixElem{T}) where {T <: RingElement}
   is_square(M) || throw(DomainError(M, "Can not invert non-square Matrix"))
   flag, X, d = _can_solve_with_solution_fflu(M, identity_matrix(M))
   !flag && error("Singular matrix in pseudo_inv")
   return X, d
end

function Base.inv(M::MatrixElem{T}) where {T <: FieldElement}
   is_square(M) || throw(DomainError(M, "Can not invert non-square Matrix"))
   flag, A = can_solve_with_solution(M, identity_matrix(M))
   !flag && error("Singular matrix in inv")
   return A
end

@doc raw"""
    inv(M::MatrixElem{T}) where {T <: RingElement}

Given a non-singular $n\times n$ matrix over a ring, return an
$n\times n$ matrix $X$ such that $MX = I_n$, where $I_n$ is the $n\times n$
identity matrix. If $M$ is not invertible over the base ring an exception is
raised.
"""
function Base.inv(M::MatrixElem{T}) where {T <: RingElement}
   is_square(M) || throw(DomainError(M, "Can not invert non-square Matrix"))
   X, d = pseudo_inv(M)
   is_unit(d) || throw(DomainError(M, "Matrix is not invertible."))
   return divexact(X, d)
end

###############################################################################
#
#   Is invertible
#
###############################################################################

@doc raw"""
    is_invertible_with_inverse(A::MatrixElem{T}; side::Symbol = :left) where {T <: RingElement}

Given an $n \times m$ matrix $A$ over a ring, return a tuple `(flag, B)`. If
`side` is `:right` and `flag` is `true`, $B$ is a right inverse of $A$ i.e.
$A B$ is the $n \times n$ unit matrix. If `side` is `:left` and `flag` is
`true`, $B$ is a left inverse of $A$ i.e. $B A$ is the $m \times m$ unit matrix.
If `flag` is `false`, no right or left inverse exists.

To get the space of all inverses, note that if $B$ and $C$ are both right
inverses, then $A (B - C) = 0$, and similar for left inverses. Hence from one
inverse one can find all by making suitable use of [`kernel`](@ref).
"""
function is_invertible_with_inverse(A::MatrixElem{T}; side::Symbol = :left) where {T <: RingElement}
   if (side == :left && nrows(A) < ncols(A)) || (side == :right && ncols(A) < nrows(A))
      return (false, zero(A, 0, 0))
   end
   I = (side == :left) ? zero(A, ncols(A), ncols(A)) : zero(A, nrows(A), nrows(A))
   for i = 1:ncols(I)
      I[i, i] = one(base_ring(I))
   end
   return can_solve_with_solution(A, I; side = side)
end

@doc raw"""
    is_invertible(A::MatrixElem{T}) where {T <: RingElement}

Return true if a given square matrix is invertible, false otherwise. If
the inverse should also be computed, use `is_invertible_with_inverse`.
"""
is_invertible(A::MatrixElem{T}) where {T <: RingElement} = is_square(A) && is_unit(det(A))

is_invertible(A::MatrixElem{T}) where {T <: FieldElement} = nrows(A) == ncols(A) == rank(A)

###############################################################################
#
#   Nullspace
#
###############################################################################

@doc raw"""
    nullspace(M::MatElem{T}) where {T <: RingElement}

Return a tuple $(\nu, N)$ consisting of the nullity $\nu$ of $M$ and
a basis $N$ (consisting of column vectors) for the right nullspace of $M$,
i.e. such that $MN$ is the zero matrix. If $M$ is an $m\times n$ matrix
$N$ will be an $n\times \nu$ matrix. Note that the nullspace is taken to be
the vector space kernel over the fraction field of the base ring if the
latter is not a field. In AbstractAlgebra we use the name "kernel" for a
function to compute an integral kernel.

# Examples

```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> S = matrix_space(R, 4, 4)
Matrix space of 4 rows and 4 columns
  over univariate polynomial ring in x over integers

julia> M = S([-6*x^2+6*x+12 -12*x^2-21*x-15 -15*x^2+21*x+33 -21*x^2-9*x-9;
              -8*x^2+8*x+16 -16*x^2+38*x-20 90*x^2-82*x-44 60*x^2+54*x-34;
              -4*x^2+4*x+8 -8*x^2+13*x-10 35*x^2-31*x-14 22*x^2+21*x-15;
              -10*x^2+10*x+20 -20*x^2+70*x-25 150*x^2-140*x-85 105*x^2+90*x-50])
[  -6*x^2 + 6*x + 12   -12*x^2 - 21*x - 15    -15*x^2 + 21*x + 33     -21*x^2 - 9*x - 9]
[  -8*x^2 + 8*x + 16   -16*x^2 + 38*x - 20     90*x^2 - 82*x - 44    60*x^2 + 54*x - 34]
[   -4*x^2 + 4*x + 8    -8*x^2 + 13*x - 10     35*x^2 - 31*x - 14    22*x^2 + 21*x - 15]
[-10*x^2 + 10*x + 20   -20*x^2 + 70*x - 25   150*x^2 - 140*x - 85   105*x^2 + 90*x - 50]

julia> n, N = nullspace(M)
(2, [1320*x^4-330*x^2-1320*x-1320 1056*x^4+1254*x^3+1848*x^2-66*x-330; -660*x^4+1320*x^3+1188*x^2-1848*x-1056 -528*x^4+132*x^3+1584*x^2+660*x-264; 396*x^3-396*x^2-792*x 0; 0 396*x^3-396*x^2-792*x])
```
"""
function nullspace(M::MatElem{T}) where {T <: RingElement}
   n = ncols(M)
   rank, A, d = rref_rational(M)
   nullity = n - rank
   R = base_ring(M)
   U = zero(M, n, nullity)
   if rank == 0
      for i = 1:nullity
         U[i, i] = one(R)
      end
   elseif nullity != 0
      pivots = zeros(Int, rank)
      nonpivots = zeros(Int, nullity)
      j = k = 1
      for i = 1:rank
         while is_zero_entry(A, i, j)
            nonpivots[k] = j
            j += 1
            k += 1
         end
         pivots[i] = j
         j += 1
      end
      while k <= nullity
         nonpivots[k] = j
         j += 1
         k += 1
      end
      d = -A[1, pivots[1]]
      for i = 1:nullity
         for j = 1:rank
            U[pivots[j], i] = A[j, nonpivots[i]]
         end
         U[nonpivots[i], i] = d
      end
   end
   return nullity, U
end

@doc raw"""
    nullspace(M::MatElem{T}) where {T <: FieldElement}

Return a tuple $(\nu, N)$ consisting of the nullity $\nu$ of $M$ and
a basis $N$ (consisting of column vectors) for the right nullspace of $M$,
i.e. such that $MN$ is the zero matrix. If $M$ is an $m\times n$ matrix
$N$ will be an $n\times \nu$ matrix.
"""
function nullspace(M::MatElem{T}) where {T <: FieldElement}
   m = nrows(M)
   n = ncols(M)
   rank, A = rref(M)
   nullity = n - rank
   R = base_ring(M)
   X = zero(M, n, nullity)
   if rank == 0
      for i = 1:nullity
         X[i, i] = one(R)
      end
   elseif nullity != 0
      pivots = zeros(Int, max(m, n))
      np = rank
      j = k = 1
      for i = 1:rank
         while is_zero_entry(A, i, j)
            pivots[np + k] = j
            j += 1
            k += 1
         end
         pivots[i] = j
         j += 1
      end
      while k <= nullity
         pivots[np + k] = j
         j += 1
         k += 1
      end
      for i = 1:nullity
         for j = 1:rank
            X[pivots[j], i] = -A[j, pivots[np + i]]
         end
         X[pivots[np + i], i] = one(R)
      end
   end
   return nullity, X
end

###############################################################################
#
#   Nilpotency
#
###############################################################################

@doc raw"""
    is_nilpotent(A::MatrixElem{T}) where {T <: RingElement}

Return if `A` is nilpotent, i.e. if there exists a natural number $k$
such that $A^k = 0$. If `A` is not square an exception is raised.
"""
function is_nilpotent(A::MatrixElem{T}) where {T <: RingElement}
  is_domain_type(T) || error("Only supported over integral domains")
  !is_square(A) && error("Dimensions don't match in is_nilpotent")
  is_zero(tr(A)) || return false
  n = nrows(A)
  A = deepcopy(A)
  i = 1
  is_zero(A) && return true
  while i < n
    i *= 2
    A = mul!(A, A, A)
    is_zero(A) && return true
  end
  return false
end

###############################################################################
#
#   Hessenberg form
#
###############################################################################

function hessenberg!(A::MatrixElem{T}) where {T <: RingElement}
   !is_square(A) && error("Dimensions don't match in hessenberg")
   R = base_ring(A)
   n = nrows(A)
   u = R()
   t = R()
   for m = 2:n - 1
      i = m + 1
      while i <= n && is_zero_entry(A, i, m - 1)
         i += 1
      end
      if i != n + 1
         if !is_zero_entry(A, m, m - 1)
            i = m
         end
         h = -inv(A[i, m - 1])
         if i > m
            for j = m - 1:n
               A[i, j], A[m, j] = A[m, j], A[i, j]
            end
            for j = 1:n
               A[j, i], A[j, m] = A[j, m], A[j, i]
            end
         end
         for i = m + 1:n
            if !is_zero_entry(A, i, m - 1)
               u = mul!(u, A[i, m - 1], h)
               for j = m:n
                  t = mul!(t, u, A[m, j])
                  A[i, j] = add!(A[i, j], t)
               end
               u = -u
               for j = 1:n
                  t = mul!(t, u, A[j, i])
                  A[j, m] = add!(A[j, m], t)
               end
               A[i, m - 1] = R()
            end
         end
      end
   end
end

@doc raw"""
    hessenberg(A::MatrixElem{T}) where {T <: RingElement}

Return the Hessenberg form of $M$, i.e. an upper Hessenberg matrix
which is similar to $M$. The upper Hessenberg form has nonzero entries
above and on the diagonal and in the diagonal line immediately below the
diagonal.
"""
function hessenberg(A::MatrixElem{T}) where {T <: RingElement}
   !is_square(A) && error("Dimensions don't match in hessenberg")
   M = deepcopy(A)
   hessenberg!(M)
   return M
end

@doc raw"""
    is_hessenberg(A::MatrixElem{T}) where {T <: RingElement}

Return `true` if $M$ is in Hessenberg form, otherwise returns `false`.
"""
function is_hessenberg(A::MatrixElem{T}) where {T <: RingElement}
   if !is_square(A)
      return false
   end
   n = nrows(A)
   for i = 3:n
      for j = 1:i - 2
         if !is_zero_entry(A, i, j)
            return false
         end
      end
   end
   return true
end

###############################################################################
#
#   Characteristic polynomial
#
###############################################################################

function charpoly_hessenberg!(S::Ring, A::MatrixElem{T}) where {T <: RingElement}
   !is_square(A) && error("Dimensions don't match in charpoly")
   R = base_ring(A)
   base_ring(S) != base_ring(A) && error("Cannot coerce into polynomial ring")
   n = nrows(A)
   if n == 0
      return one(S)
   end
   if n == 1
      return gen(S) - A[1, 1]
   end
   hessenberg!(A)
   P = Vector{elem_type(S)}(undef, n + 1)
   P[1] = one(S)
   x = gen(S)
   for m = 1:n
      P[m + 1] = (x - A[m, m])*P[m]
      t = one(R)
      for i = 1:m - 1
         t = mul!(t, t, A[m - i + 1, m - i])
         P[m + 1] -= t*A[m - i, m]*P[m - i]
      end
   end
   return P[n + 1]
end

function charpoly_danilevsky_ff!(S::Ring, A::MatrixElem{T}) where {T <: RingElement}
   !is_square(A) && error("Dimensions don't match in charpoly")
   R = base_ring(A)
   base_ring(S) != base_ring(A) && error("Cannot coerce into polynomial ring")
   n = nrows(A)
   if n == 0
      return one(S)
   end
   if n == 1
      return gen(S) - A[1, 1]
   end
   d = one(R)
   t = R()
   V = Vector{T}(undef, n)
   W = Vector{T}(undef, n)
   pol = one(S)
   i = 1
   while i < n
      h = A[n - i + 1, n - i]
      while iszero(h)
         k = 1
         while k < n - i && is_zero_entry(A, n - i + 1, n - i - k)
            k += 1
         end
         if k == n - i
            b = S()
            fit!(b, i + 1)
            b = setcoeff!(b, i, one(R))
            for kk = 1:i
               b = setcoeff!(b, kk - 1, -A[n - i + 1, n - kk + 1]*d)
            end
            pol *= b
            n -= i
            i = 1
            if n == 1
               pol *= (gen(S) - A[1, 1]*d)
               return pol
            end
         else
            for j = 1:n
               A[n - i - k, j], A[n - i, j] = A[n - i, j], A[n - i - k, j]
            end
            for j = 1:n - i + 1
               A[j, n - i - k], A[j, n - i] = A[j, n - i], A[j, n - i - k]
            end
         end
         h = A[n - i + 1, n - i]
      end
      for j = 1:n
         V[j] = -A[n - i + 1, j]
         W[j] = deepcopy(A[n - i + 1, j])
      end
      for j = 1:n - i
         for kk = 1:n - i - 1
            t = mul_red!(t, A[j, n - i], V[kk], false)
            A[j, kk] = mul_red!(A[j, kk], A[j, kk], h, false)
            A[j, kk] = add!(A[j, kk], t)
            A[j, kk] = reduce!(A[j, kk])
         end
         for kk = n - i + 1:n
            t = mul_red!(t, A[j, n - i], V[kk], false)
            A[j, kk] = mul_red!(A[j, kk], A[j, kk], h, false)
            A[j, kk] = add!(A[j, kk], t)
            A[j, kk] = reduce!(A[j, kk])
         end
      end
      for kk = 1:n
         A[n - i + 1, kk] = R()
      end
      for j = 1:n - i
         for kk = 1:n - i - 1
            A[j, kk] = mul!(A[j, kk], A[j, kk], d)
         end
         for kk = n - i + 1:n
            A[j, kk] = mul!(A[j, kk], A[j, kk], d)
         end
      end
      A[n - i + 1, n - i] = deepcopy(h)
      for j = 1:n - i - 1
         s = R()
         for kk = 1:n - i
            s = addmul_delayed_reduction!(s, A[kk, j], W[kk], t)
         end
         s = reduce!(s)
         A[n - i, j] = s
      end
      for j = n - i:n - 1
         s = R()
         for kk = 1:n - i
            s = addmul_delayed_reduction!(s, A[kk, j], W[kk], t)
         end
         s = addmul_delayed_reduction!(s, h, W[j + 1], t)
         s = reduce!(s)
         A[n - i, j] = s
      end
      s = R()
      for kk = 1:n - i
         s = addmul_delayed_reduction!(s, A[kk, n], W[kk], t)
      end
      s = reduce!(s)
      A[n - i, n] = s
      for kk = 1:n
         A[n - i, kk] = mul!(A[n - i, kk], A[n - i, kk], d)
      end
      d = inv(h)
      parent(d) # To work around a bug in julia
      i += 1
   end
   b = S()
   fit!(b, n + 1)
   b = setcoeff!(b, n, one(R))
   for i = 1:n
      c = -A[1, n - i + 1]*d
      b = setcoeff!(b, i - 1, c)
   end
   return pol*b
end

function charpoly_danilevsky!(S::Ring, A::MatrixElem{T}) where {T <: RingElement}
   !is_square(A) && error("Dimensions don't match in charpoly")
   R = base_ring(A)
   base_ring(S) != base_ring(A) && error("Cannot coerce into polynomial ring")
   n = nrows(A)
   if n == 0
      return one(S)
   end
   if n == 1
      return gen(S) - A[1, 1]
   end
   t = R()
   V = Vector{T}(undef, n)
   W = Vector{T}(undef, n)
   pol = one(S)
   i = 1
   while i < n
      h = A[n - i + 1, n - i]
      while iszero(h)
         k = 1
         while k < n - i && is_zero_entry(A, n - i + 1, n - i - k)
            k += 1
         end
         if k == n - i
            b = S()
            fit!(b, i + 1)
            b = setcoeff!(b, i, one(R))
            for kk = 1:i
               b = setcoeff!(b, kk - 1, -A[n - i + 1, n - kk + 1])
            end
            pol *= b
            n -= i
            i = 1
            if n == 1
               pol *= (gen(S) - A[1, 1])
               return pol
            end
         else
            for j = 1:n
               A[n - i - k, j], A[n - i, j] = A[n - i, j], A[n - i - k, j]
            end
            for j = 1:n - i + 1
               A[j, n - i - k], A[j, n - i] = A[j, n - i], A[j, n - i - k]
            end
         end
         h = A[n - i + 1, n - i]
      end
      h = -inv(h)
      for j = 1:n
         V[j] = A[n - i + 1, j]*h
         W[j] = deepcopy(A[n - i + 1, j])
      end
      h = -h
      for j = 1:n - i
         for k = 1:n - i - 1
            t = mul!(t, A[j, n - i], V[k])
            A[j, k] = add!(A[j, k], t)
         end
         for k = n - i + 1:n
            t = mul!(t, A[j, n - i], V[k])
            A[j, k] = add!(A[j, k], t)
         end
         A[j, n - i] = mul!(A[j, n - i], A[j, n - i], h)
      end
      for j = 1:n - i - 1
         s = R()
         for k = 1:n - i
            s = addmul_delayed_reduction!(s, A[k, j], W[k], t)
         end
         s = reduce!(s)
         A[n - i, j] = s
      end
      for j = n - i:n - 1
         s = R()
         for k = 1:n - i
            s = addmul_delayed_reduction!(s, A[k, j], W[k], t)
         end
         s = add!(s, W[j + 1])
         s = reduce!(s)
         A[n - i, j] = s
      end
      s = R()
      for k = 1:n - i
         s = addmul_delayed_reduction!(s, A[k, n], W[k], t)
      end
      s = reduce!(s)
      A[n - i, n] = s
      i += 1
   end
   b = S()
   fit!(b, n + 1)
   b = setcoeff!(b, n, one(R))
   for i = 1:n
      b = setcoeff!(b, i - 1, -A[1, n - i + 1])
   end
   return pol*b
end

@doc raw"""
    charpoly(Y::MatrixElem{T}) where {T <: RingElement}
    charpoly(S::PolyRing{T}, Y::MatrixElem{T}) where {T <: RingElement}

Return the characteristic polynomial $p$ of the square matrix $Y$.
If a polynomial ring $S$ over the same base ring as $Y$ is supplied,
the resulting polynomial is an element of it.

# Examples

```jldoctest
julia> R, = residue_ring(ZZ, 7);

julia> S = matrix_space(R, 4, 4)
Matrix space of 4 rows and 4 columns
  over residue ring of integers modulo 7

julia> T, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

julia> M = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);
              R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])
[1   2   4   3]
[2   5   1   0]
[6   1   3   2]
[1   1   3   5]

julia> A = charpoly(T, M)
y^4 + 2*y^2 + 6*y + 2

julia> A = charpoly(M)
x^4 + 2*x^2 + 6*x + 2
```
"""
function charpoly(S::PolyRing{T}, Y::MatrixElem{T}) where {T <: RingElement}
   !is_square(Y) && error("Dimensions don't match in charpoly")
   R = base_ring(Y)
   base_ring(S) != base_ring(Y) && error("Cannot coerce into polynomial ring")
   n = nrows(Y)
   if n == 0
      return one(S)
   end
   F = Vector{elem_type(R)}(undef, n)
   A = Vector{elem_type(R)}(undef, n)
   M = Matrix{elem_type(R)}(undef, n - 1, n)
   F[1] = -Y[1, 1]
   for i = 2:n
      F[i] = R()
      for j = 1:i
         M[1, j] = Y[j, i]
      end
      A[1] = Y[i, i]
      p = R()
      for j = 2:i - 1
         for k = 1:i
            s = R()
            for l = 1:i
               s = addmul_delayed_reduction!(s, Y[k, l], M[j - 1, l], p)
            end
            s = reduce!(s)
            M[j, k] = s
         end
         A[j] = M[j, i]
      end
      s = R()
      for j = 1:i
         s = addmul_delayed_reduction!(s, Y[i, j], M[i - 1, j], p)
      end
      s = reduce!(s)
      A[i] = s
      for j = 1:i
         s = -F[j]
         for k = 1:j - 1
            s = addmul_delayed_reduction!(s, A[k], F[j - k], p)
         end
         s = reduce!(s)
         F[j] = -s - A[j]
     end
   end
   z = gen(S)
   f = z^n
   for i = 1:n
      f = setcoeff!(f, n - i, F[i])
   end
   return f
end

function charpoly(Y::MatrixElem)
   R = base_ring(Y)
   Rx, x = polynomial_ring(R; cached=false)
   return charpoly(Rx, Y)
end

###############################################################################
#
#   Minimal polynomial
#
###############################################################################

# We let the rows of A be the Krylov sequence v, Mv, M^2v, ... where v
# is a standard basis vector. We find a minimal polynomial P_v by finding
# a linear relation amongst the rows (rows are added one at a time until
# a relation is found). We then append the rows from A to B and row reduce
# by all the rows already in B. We then look for a new standard basis vector
# v not in the span of the rows in B and repeat the above. The arrays P1 and
# P2 are because we can't efficiently swap rows. The i-th entry encodes the
# row of A (in the case of P1) or B (in the case of P2) which has a pivot in
# column i. The arrays L1 and L2 encode the number of columns that contain
# nonzero entries for each row of A and B respectively. These are required
# by the reduce_row! function. If charpoly_only is set to true, the function
# will return just the polynomial obtained from the first Krylov subspace.
# If the charpoly is the minpoly, this returned polynomial will be the
# charpoly iff it has degree n. Otherwise it is meaningless (but it is
# extremely fast to compute over some fields).

function minpoly(S::PolyRing{T}, M::MatElem{T}, charpoly_only::Bool = false) where {T <: FieldElement}
   !is_square(M) && error("Not a square matrix in minpoly")
   base_ring(S) != base_ring(M) && error("Unable to coerce polynomial")
   n = nrows(M)
   if n == 0
      return one(S)
   end
   R = base_ring(M)
   p = one(S)
   A = similar(M, n + 1, 2n + 1)
   B = similar(M, n, n)
   L1 = Int[n + i for i in 1:n + 1]
   L2 = Int[n for i in 1:n]
   P2 = zeros(Int, n)
   P2[1] = 1
   c2 = 1
   r2 = 1
   first_poly = true
   while r2 <= n
      P1 = Int[0 for i in 1:2n + 1]
      v = zero(M, n, 1)
      for j = 1:n
         B[r2, j] = v[j, 1]
         A[1, j] = R()
      end
      P1[c2] = 1
      P2[c2] = r2
      v[c2, 1] = one(R)
      B[r2, c2] = v[c2, 1]
      A[1, c2] = one(R)
      A[1, n + 1] = one(R)
      indep = true
      r1 = 1
      c1 = 0
      while c1 <= n && r1 <= n
         r1 += 1
         r2 = indep ? r2 + 1 : r2
         v = M*v
         for j = 1:n
            A[r1, j] = deepcopy(v[j, 1])
         end
         for j = n + 1:n + r1 - 1
            A[r1, j] = zero(R)
         end
         A[r1, n + r1] = one(R)
         c1 = reduce_row!(A, P1, L1, r1)
         if indep && r2 <= n && !first_poly
            for j = 1:n
               B[r2, j] = deepcopy(v[j, 1])
            end
            c = reduce_row!(B, P2, L2, r2)
            indep = c != 0
         end
      end
      if first_poly
         for j = 1:n
            P2[j] = P1[j]
         end
         r2 = r1
      end
      c = 0
      for j = c2 + 1:n
         if P2[j] == 0
            c = j
            break
         end
      end
      c2 = c
      b = S()
      fit!(b, r1)
      h = inv(A[r1, n + r1])
      for i = 1:r1
         b = setcoeff!(b, i - 1, A[r1, n + i]*h)
      end
      p = lcm(p, b)
      if charpoly_only == true
         return p
      end
      if first_poly && r2 <= n
         for j = 1:r1 - 1
            for k = 1:n
               B[j, k] = A[j, k]
            end
         end
      end
      first_poly = false
   end
   return p
end

@doc raw"""
    minpoly(M::MatElem{T}) where {T <: RingElement}
    minpoly(S::PolyRing{T}, M::MatElem{T}) where {T <: RingElement}

Return the minimal polynomial $p$ of the square matrix $M$.
If a polynomial ring $S$ over the same base ring as $Y$ is supplied,
the resulting polynomial is an element of it.

# Examples

```jldoctest
julia> R = GF(13)
Finite field F_13

julia> S, y = polynomial_ring(R, :y)
(Univariate polynomial ring in y over R, y)

julia> M = R[7 6 1;
             7 7 5;
             8 12 5]
[7    6   1]
[7    7   5]
[8   12   5]

julia> A = minpoly(S, M)
y^2 + 10*y

julia> A = minpoly(M)
x^2 + 10*x

```
"""
function minpoly(S::PolyRing{T}, M::MatElem{T}, charpoly_only::Bool = false) where {T <: RingElement}
   !is_square(M) && error("Not a square matrix in minpoly")
   base_ring(S) != base_ring(M) && error("Unable to coerce polynomial")
   n = nrows(M)
   if n == 0
      return one(S)
   end
   R = base_ring(M)
   p = one(S)
   A = similar(M, n + 1, 2n + 1)
   B = similar(M, n, n)
   L1 = zeros(Int, n + 1)
   for i in 1:n + 1
      L1[i] = i + n
   end
   L2 = fill(n, n)
   P2 = zeros(Int, n)
   P2[1] = 1
   c2 = 1
   r2 = 1
   first_poly = true
   while r2 <= n
      P1 = [0 for i in 1:2n + 1]
      v = zero(M, n, 1)
      for j = 1:n
         B[r2, j] = v[j, 1]
         A[1, j] = R()
      end
      P1[c2] = 1
      P2[c2] = r2
      v[c2, 1] = one(R)
      B[r2, c2] = v[c2, 1]
      for s = 1:c2 - 1
         if P2[s] != 0
            B[r2, c2] *= B[P2[s], s]
         end
      end
      A[1, c2] = one(R)
      A[1, n + 1] = one(R)
      indep = true
      r1 = 1
      c1 = 0
      while c1 <= n && r1 <= n
         r1 += 1
         r2 = indep ? r2 + 1 : r2
         v = M*v
         for j = 1:n
            A[r1, j] = deepcopy(v[j, 1])
         end
         for j = n + 1:n + r1 - 1
            A[r1, j] = zero(R)
         end
         A[r1, n + r1] = one(R)
         c1 = reduce_row!(A, P1, L1, r1)
         if indep && r2 <= n && !first_poly
            for j = 1:n
               B[r2, j] = deepcopy(v[j, 1])
            end
            c = reduce_row!(B, P2, L2, r2)
            indep = c != 0
         end
      end
      if first_poly
         for j = 1:n
            P2[j] = P1[j]
         end
         r2 = r1
      end
      c = 0
      for j = c2 + 1:n
         if P2[j] == 0
            c = j
            break
         end
      end
      c2 = c
      b = S()
      fit!(b, r1)
      for i = 1:r1
         b = setcoeff!(b, i - 1, A[r1, n + i])
      end
      b = reverse(b, r1)
      b = primpart(b)
      b = reverse(b, r1)
      p = lcm(p, b)
      if charpoly_only == true
         return divexact(p, canonical_unit(p))
      end
      if first_poly && r2 <= n
         for j = 1:r1 - 1
            for k = 1:n
               B[j, k] = A[j, k]
            end
         end
      end
      first_poly = false
   end
   return divexact(p, canonical_unit(p))
end

function minpoly(M::MatElem{T}, charpoly_only::Bool = false) where {T <: RingElement}
   R = base_ring(M)
   Rx, x = polynomial_ring(R; cached=false)
   return minpoly(Rx, M, charpoly_only)
end

###############################################################################
#
#   Hermite Normal Form
#
###############################################################################

function hnf_cohen(A::MatrixElem{T}) where {T <: RingElement}
   H, U = hnf_cohen_with_transform(A)
   return H
end

function hnf_cohen_with_transform(A::MatrixElem{T}) where {T <: RingElement}
   H = deepcopy(A)
   m = nrows(H)
   U = identity_matrix(A, m)
   hnf_cohen!(H, U)
   return H, U
end

function hnf_cohen!(H::MatrixElem{T}, U::MatrixElem{T}) where {T <: RingElement}
   m = nrows(H)
   n = ncols(H)
   l = min(m, n)
   k = 1
   t = base_ring(H)()
   t1 = base_ring(H)()
   t2 = base_ring(H)()
   for i = 1:l
      for j = k + 1:m
         if is_zero_entry(H, j, i)
            continue
         end
         d, u, v = gcdx(H[k, i], H[j, i])
         a = divexact(H[k, i], d)
         b = -divexact(H[j, i], d)
         for c = i:n
            t = deepcopy(H[j, c])
            t1 = mul_red!(t1, a, H[j, c], false)
            t2 = mul_red!(t2, b, H[k, c], false)
            H[j, c] = t1 + t2
            H[j, c] = reduce!(H[j, c])
            t1 = mul_red!(t1, u, H[k, c], false)
            t2 = mul_red!(t2, v, t, false)
            H[k, c] = t1 + t2
            H[k, c] = reduce!(H[k, c])
         end
         for c = 1:m
            t = deepcopy(U[j,c])
            t1 = mul_red!(t1, a, U[j, c], false)
            t2 = mul_red!(t2, b, U[k, c], false)
            U[j, c] = t1 + t2
            U[j, c] = reduce!(U[j, c])
            t1 = mul_red!(t1, u, U[k, c], false)
            t2 = mul_red!(t2, v, t, false)
            U[k, c] = t1 + t2
            U[k, c] = reduce!(U[k, c])
         end
      end
      if is_zero_entry(H, k, i)
         continue
      end
      cu = canonical_unit(H[k, i])
      if !isone(cu)
         for c = i:n
            H[k, c] = divexact(H[k, c], cu)
        end
         for c = 1:m
            U[k, c] = divexact(U[k, c], cu)
         end
      end
      for j = 1:k-1
         q = -div(H[j,i], H[k, i])
         for c = i:n
            t = mul!(t, q, H[k, c])
            H[j, c] = add!(H[j, c], t)
         end
         for c = 1:m
            t = mul!(t, q, U[k, c])
            U[j, c] = add!(U[j, c], t)
         end
      end
      k += 1
   end
   return nothing
end

#  Hermite normal form via Kannan-Bachem for matrices of full column rank
#
#  This is the algorithm of Kannan, Bachem, "Polynomial algorithms for computing
#  the Smith and Hermite normal forms of an integer matrix", Siam J. Comput.,
#  Vol. 8, No. 4, pp. 499-507.

@doc raw"""
    hnf_minors(A::MatrixElem{T}) where {T <: RingElement}

Compute the upper right row Hermite normal form of $A$ using the algorithm of
Kannan-Bachem. The input must have full column rank.
"""
function hnf_minors(A::MatrixElem{T}) where {T <: RingElement}
   H = deepcopy(A)
   _hnf_minors!(H, similar(A, 0, 0), Val(false))
   return H
end

@doc raw"""
    hnf_minors_with_transform(A::MatrixElem{T}) where {T <: RingElement}

Compute the upper right row Hermite normal form $H$ of $A$ and an invertible
matrix $U$ with $UA = H$ using the algorithm of Kannan-Bachem. The input must
have full column rank.
"""
function hnf_minors_with_transform(A::MatrixElem{T}) where {T <: RingElement}
   H = deepcopy(A)
   U = similar(A, nrows(A), nrows(A))
   _hnf_minors!(H, U, Val(true))
   return H, U
end

function _hnf_minors!(H::MatrixElem{T}, U::MatrixElem{T}, ::Val{with_transform} = Val(false)) where {T <: RingElement, with_transform}
   m = nrows(H)
   n = ncols(H)

   l = m
   k = 0

   R = base_ring(H)

   if with_transform
      for i in 1:m
         for j in 1:m
            if j == i
               U[i, j] = one(R)
            else
               U[i, j] = zero(R)
            end
         end
      end
   end

   t = zero(R)
   b = zero(R)
   t2 = zero(R)
   r1d = zero(R)
   r2d = zero(R)
   q = zero(R)
   minus_one = base_ring(H)(-1)

   while k <= n - 1
      k = k + 1
      if k == 1 && !is_zero_entry(H, k, k)
         for j2 in 1:n
            H[k, j2] = deepcopy(H[k, j2])
         end
      end

      # We have a matrix of form ( * * * * )
      #                          ( 0 * * * )
      #                          ( 0 0 * * ) <- k - 1
      #                          ( * * * * ) <- k
      # We want to produce zeroes in the first k entries
      # of row k.

      for j in 1:k-1
         # Since we want to mutate the k-th row, first make a copy
         if j == 1
            for j2 in j:n
               H[k, j2] = deepcopy(H[k, j2])
            end
         end

         # Shortcuts
         if is_zero_entry(H, k, j)
            continue
         end

         di, q = divides(H[k, j], H[j, j])
         if di
            q = mul!(q, q, minus_one)
            for j2 in j:n
               t = mul!(t, q, H[j, j2])
               H[k, j2] = add!(H[k, j2], H[k, j2], t)
            end
            if with_transform
               for j2 in 1:m
                  t = mul!(t, q, U[j, j2])
                  U[k, j2] = add!(U[k, j2], U[k, j2], t)
               end
            end
            continue
         end

         # Generic case
         d, u, v = gcdx(H[j, j], H[k, j])

         r1d = divexact(H[j, j], d)
         r2d = divexact(H[k, j], d)

         r2d = mul!(r2d, r2d, minus_one)
         for j2 in j:n
            b = mul_red!(b, u, H[j, j2], false)
            t2 = mul_red!(t2, v, H[k, j2], false)
            H[k, j2] = mul_red!(H[k, j2], H[k, j2], r1d, false)
            t = mul_red!(t, r2d, H[j, j2], false)
            H[k, j2] = add!(H[k, j2], H[k, j2], t)
            H[j, j2] = add!(H[j, j2], b, t2)
            H[k, j2] = reduce!(H[k, j2])
            H[j, j2] = reduce!(H[j, j2])
         end
         if with_transform
            for j2 in 1:m
               b = mul_red!(b, u, U[j, j2], false)
               t2 = mul_red!(t2, v, U[k, j2], false)
               U[k, j2] = mul_red!(U[k, j2], U[k, j2], r1d, false)
               t = mul_red!(t, r2d, U[j, j2], false)
               U[k, j2] = add!(U[k, j2], U[k, j2], t)
               U[j, j2] = add!(U[j, j2], b, t2)
               U[k, j2] = reduce!(U[k, j2])
               U[j, j2] = reduce!(U[j, j2])
            end
         end
      end

      if is_zero_entry(H, k, k)
         swap_rows!(H, k, l)
         if with_transform
            swap_rows!(U, k, l)
         end
         l = l - 1
         k = k - 1
         continue
      end

      u = canonical_unit(H[k, k])
      if !isone(u)
         u = R(inv(u))
         for j in k:n
            H[k, j] = mul!(H[k, j], H[k, j], u)
         end
         if with_transform
            for j in 1:m
               U[k, j] = mul!(U[k, j], U[k, j], u)
            end
         end
      end

      for i in (k - 1):-1:1
         for j in (i + 1):k
            q = div(H[i, j], H[j, j])
            if iszero(q)
              continue
            end
            q = mul!(q, q, minus_one)
            for j2 in j:n
               t = mul!(t, q, H[j, j2])
               H[i, j2] = add!(H[i, j2], H[i, j2], t)
            end
            if with_transform
               for j2 in 1:m
                  t = mul!(t, q, U[j, j2])
                  U[i, j2] = add!(U[i, j2], U[i, j2], t)
               end
            end
         end
      end
      l = m
   end

   # Now the matrix is of form
   # ( * * * )
   # ( 0 * * )
   # ( 0 0 * )
   # ( * * * ) <- n + 1
   # ( * * * )

   for k in (n + 1):m
      for j in 1:n
         if j == 1
            for j2 in 1:n
               H[k, j2] = deepcopy(H[k, j2])
            end
         end

         # We do the same shortcuts as above
         if is_zero_entry(H, k, j)
            continue
         end

         di, q = divides(H[k, j], H[j, j])
         if di
            q = mul!(q, q, minus_one)
            for j2 in j:n
               t = mul!(t, q, H[j, j2])
               H[k, j2] = add!(H[k, j2], H[k, j2], t)
            end
            if with_transform
               for j2 in 1:m
                  t = mul!(t, q, U[j, j2])
                  U[k, j2] = add!(U[k, j2], U[k, j2], t)
               end
            end
            continue
         end

         d, u, v = gcdx(H[j, j], H[k, j])
         r1d = divexact(H[j, j], d)
         r2d = divexact(H[k, j], d)
         r2d = mul!(r2d, r2d, minus_one)
         for j2 in j:n
            b = mul_red!(b, u, H[j, j2], false)
            t2 = mul_red!(t2, v, H[k, j2], false)
            H[k, j2] = mul_red!(H[k, j2], H[k, j2], r1d, false)
            t = mul_red!(t, r2d, H[j, j2], false)
            H[k, j2] = add!(H[k, j2], H[k, j2], t)
            H[j, j2] = add!(H[j, j2], b, t2)
            H[k, j2] = reduce!(H[k, j2])
            H[j, j2] = reduce!(H[j, j2])
         end
         if with_transform
            for j2 in 1:m
               b = mul_red!(b, u, U[j, j2], false)
               t2 = mul_red!(t2, v, U[k, j2], false)
               U[k, j2] = mul_red!(U[k, j2], U[k, j2], r1d, false)
               t = mul_red!(t, r2d, U[j, j2], false)
               U[k, j2] = add!(U[k, j2], U[k, j2], t)
               U[j, j2] = add!(U[j, j2], b, t2)
               U[k, j2] = reduce!(U[k, j2])
               U[j, j2] = reduce!(U[j, j2])
            end
         end
      end
      for i in n:-1:1
         for j in (i + 1):n
            q = div(H[i, j], H[j, j])
            if iszero(q)
              continue
            end
            q = mul!(q, q, minus_one)
            for j2 in j:n
               t = mul!(t, q, H[j, j2])
               H[i, j2] = add!(H[i, j2], H[i, j2], t)
            end
            if with_transform
               for j2 in 1:m
                  t = mul!(t, q, U[j, j2])
                  U[i, j2] = add!(U[i, j2], U[i, j2], t)
               end
            end
         end
      end
   end
   return H
end

#  Hermite normal form for arbitrary matrices via a modification of the
#  Kannan-Bachem algorithm

@doc raw"""
    hnf_kb(A::MatrixElem{T}) where {T <: RingElement}

Compute the upper right row Hermite normal form of $A$ using a modification
of the algorithm of Kannan-Bachem.
"""
function hnf_kb(A::MatrixElem{T}) where {T <: RingElement}
   return _hnf_kb(A, Val(false))
end

@doc raw"""
    hnf_kb_with_transform(A::MatrixElem{T}) where {T <: RingElement}

Compute the upper right row Hermite normal form $H$ of $A$ and an invertible
matrix $U$ with $UA = H$ using a modification of the algorithm of
Kannan-Bachem.
"""
function hnf_kb_with_transform(A::MatrixElem{T}) where {T <: RingElement}
   return _hnf_kb(A, Val(true))
end

function _hnf_kb(A, ::Val{with_transform} = Val(false)) where {with_transform}
   H = deepcopy(A)
   m = nrows(H)
   if with_transform
      U = identity_matrix(A, m)
      hnf_kb!(H, U, true)
      return H, U
   else
      U = similar(A, 0, 0)
      hnf_kb!(H, U, false)
      return H
   end
end

function kb_search_first_pivot(H, start_element::Int = 1)
   for r = start_element:nrows(H)
      for c = start_element:ncols(H)
         if !is_zero_entry(H, r, c)
            return r, c
         end
      end
   end
   return 0, 0
end

# Reduces the entries above H[pivot[c], c]
function kb_reduce_column!(H::MatrixElem{T}, U::MatrixElem{T}, pivot::Vector{Int}, c::Int, with_trafo::Bool, start_element::Int = 1) where {T <: RingElement}

   # Let c = 4 and pivot[c] = 4. H could look like this:
   # ( 0 . * # * )
   # ( . * * # * )
   # ( 0 0 0 0 . )
   # ( 0 0 0 . * )
   # ( * * * * * )
   #
   # (. are pivots, we want to reduce the entries marked with #)
   # The #'s are in rows whose pivot is in a column left of column c.

   r = pivot[c]
   t = base_ring(H)()
   for i = start_element:c - 1
      p = pivot[i]
      if p == 0
         continue
      end
      # So, the pivot in row p is in a column left of c.
      if is_zero_entry(H, p, c)
         continue
      end
      q = -div(H[p, c], H[r, c])
      for j = c:ncols(H)
         t = mul!(t, q, H[r, j])
         H[p, j] += t
      end
      if with_trafo
         for j = 1:ncols(U)
            t = mul!(t, q, U[r, j])
            U[p, j] += t
         end
      end
   end
   return nothing
end

# Multiplies row r by a unit such that the entry H[r, c] is "canonical"
function kb_canonical_row!(H, U, r::Int, c::Int, with_trafo::Bool)
   cu = canonical_unit(H[r, c])
   if !isone(cu)
      for j = c:ncols(H)
         H[r, j] = divexact(H[r, j], cu)
      end
      if with_trafo
         for j = 1:ncols(U)
            U[r, j] = divexact(U[r, j], cu)
         end
      end
   end
   return nothing
end

function kb_sort_rows!(H::MatrixElem{T}, U::MatrixElem{T}, pivot::Vector{Int}, with_trafo::Bool, start_element::Int = 1) where {T <:RingElement}
   m = nrows(H)
   n = ncols(H)
   pivot2 = zeros(Int, m)
   for i = 1:n
      if pivot[i] == 0
         continue
      end
      pivot2[pivot[i]] = i
   end

   r1 = start_element
   for i = start_element:n
      r2 = pivot[i]
      if r2 == 0
         continue
      end
      if r1 != r2
         swap_rows!(H, r1, r2)
         with_trafo ? swap_rows!(U, r1, r2) : nothing
         p = pivot2[r1]
         pivot[i] = r1
         if p != 0
            pivot[p] = r2
         end
         pivot2[r1] = i
         pivot2[r2] = p
      end
      r1 += 1
      if r1 == m
         break
      end
   end
   return nothing
end

function hnf_kb!(H::MatrixElem{T}, U::MatrixElem{T}, with_trafo::Bool = false, start_element::Int = 1) where {T <: RingElement}
   m = nrows(H)
   n = ncols(H)
   pivot = zeros(Int, n) # pivot[j] == i if the pivot of column j is in row i

   # Find the first non-zero entry of H
   row1, col1 = kb_search_first_pivot(H, start_element)
   if row1 == 0
      return nothing
   end
   pivot[col1] = row1
   kb_canonical_row!(H, U, row1, col1, with_trafo)
   pivot_max = col1
   t = base_ring(H)()
   t1 = base_ring(H)()
   t2 = base_ring(H)()
   for i = row1 + 1:m
      new_pivot = false
      for j = start_element:n
         if is_zero_entry(H, i, j)
            continue
         end
         if pivot[j] == 0
            # We found a non-zero entry in a column without a pivot: This is a
            # new pivot
            pivot[j] = i
            pivot_max = max(pivot_max, j)
            new_pivot = true
         else
            # We have a pivot for this column: Use it to write 0 in H[i, j]
            p = pivot[j]
            d, u, v = gcdx(H[p, j], H[i, j])
            a = divexact(H[p, j], d)
            b = -divexact(H[i, j], d)
            for c = j:n
               t = deepcopy(H[i, c])
               t1 = mul_red!(t1, a, H[i, c], false)
               t2 = mul_red!(t2, b, H[p, c], false)
               H[i, c] = reduce!(t1 + t2)
               t1 = mul_red!(t1, u, H[p, c], false)
               t2 = mul_red!(t2, v, t, false)
               H[p, c] = reduce!(t1 + t2)
            end
            if with_trafo
               for c = 1:m
                  t = deepcopy(U[i, c])
                  t1 = mul_red!(t1, a, U[i, c], false)
                  t2 = mul_red!(t2, b, U[p, c], false)
                  U[i, c] = reduce!(t1 + t2)
                  t1 = mul_red!(t1, u, U[p, c], false)
                  t2 = mul_red!(t2, v, t, false)
                  U[p, c] = reduce!(t1 + t2)
               end
            end
         end

         # We changed the pivot of column j (or found a new one).
         # We have do reduce the entries marked with # in
         # ( 0 0 0 . * )
         # ( . # # * * )
         # ( 0 0 . * * )
         # ( 0 . # * * )
         # ( * * * * * )
         # where . are pivots and i = 4, j = 2. (This example is for the
         # "new pivot" case.)
         kb_canonical_row!(H, U, pivot[j], j, with_trafo)
         for c = j:pivot_max
            if pivot[c] == 0
               continue
            end
            kb_reduce_column!(H, U, pivot, c, with_trafo, start_element)
         end
         if new_pivot
            break
         end
      end
   end
   kb_sort_rows!(H, U, pivot, with_trafo, start_element)
   return nothing
end

@doc raw"""
    hnf(A::MatrixElem{T}) where {T <: RingElement}

Return the upper right row Hermite normal form of $A$.
"""
function hnf(A::MatrixElem{T}) where {T <: RingElement}
  return hnf_kb(A)
end

@doc raw"""
    hnf_with_transform(A)

Return the tuple $H, U$ consisting of the upper right row Hermite normal
form $H$ of $A$ together with invertible matrix $U$ such that $UA = H$.
"""
function hnf_with_transform(A)
  return hnf_kb_with_transform(A)
end

@doc raw"""
    is_hnf(M::MatrixElem{T}) where T <: RingElement

Return `true` if the matrix is in Hermite normal form.
"""
function is_hnf(M::MatrixElem{T}) where T <: RingElement
   r = nrows(M)
   c = ncols(M)
   row = 1
   col = 1
   pivots = zeros(Int, r)
   # first check the staircase, since it is cheap to do
   while row <= r
      while col <= c && M[row, col] == 0
         for i = row + 1:r
            if M[i, col] != 0
               return false
            end
         end
         col += 1
      end
      if col <= c # found pivot
         pivots[row] = col
         for i = row + 1:r
            if M[i, col] != 0
               return false
            end
         end
      end
      row += 1
      col += 1
   end
   # now check everything above pivots is reduced (expensive)
   row = 1
   while row <= r && pivots[row] != 0
      col = pivots[row]
      p = M[row, col]
      for i = 1:row - 1
         qq, rr = divrem(M[i, col], p)
         if rr != M[i, col]
            return false
         end
      end
      row += 1
   end
   return true
end

###############################################################################
#
#   Smith Normal Form
#
###############################################################################

@doc raw"""
    is_snf(A::MatrixElem{T}) where T <: RingElement

Return `true` if $A$ is in Smith Normal Form.
"""
function is_snf(A::MatrixElem{T}) where T <: RingElement
   m = nrows(A)
   n = ncols(A)
   a = A[1, 1]
   for i = 2:min(m, n)
      q, r = divrem(A[i, i], a)
      if !iszero(r)
         return false
      end
      a = A[i,i]
   end
   for i = 1:n
      for j = 1:m
         if i == j
            continue
         end
         if !is_zero_entry(A, j, i)
            return false
         end
      end
   end
   return true
end

function snf_kb(A::MatrixElem{T}) where {T <: RingElement}
   return _snf_kb(A, Val(false))
end

function snf_kb_with_transform(A::MatrixElem{T}) where {T <: RingElement}
   return _snf_kb(A, Val(true))
end

function _snf_kb(A::MatrixElem{T}, ::Val{with_transform} = Val(false)) where {T <: RingElement, with_transform}
   S = deepcopy(A)
   m = nrows(S)
   n = ncols(S)
   if with_transform
      U = identity_matrix(A, m)
      K = identity_matrix(A, n)
      snf_kb!(S, U, K, true)
      return S, U, K
   else
      U = similar(A, 0, 0)
      K = U
      snf_kb!(S, U, K, false)
      return S
   end
end

function kb_clear_row!(S::MatrixElem{T}, K::MatrixElem{T}, i::Int, with_trafo::Bool) where {T <: RingElement}
   m = nrows(S)
   n = ncols(S)
   t = base_ring(S)()
   t1 = base_ring(S)()
   t2 = base_ring(S)()
   for j = i+1:n
      if is_zero_entry(S, i, j)
         continue
      end
      d, u, v = gcdx(S[i, i], S[i, j])
      a = divexact(S[i ,i], d)
      b = -divexact(S[i, j], d)
      for r = i:m
         t = deepcopy(S[r, j])
         t1 = mul_red!(t1, a, S[r, j], false)
         t2 = mul_red!(t2, b, S[r, i], false)
         S[r, j] = reduce!(t1 + t2)
         t1 = mul_red!(t1, u, S[r, i], false)
         t2 = mul_red!(t2, v, t, false)
         S[r, i] = reduce!(t1 + t2)
      end
      if with_trafo
         for r = 1:n
            t = deepcopy(K[r,j])
            t1 = mul_red!(t1, a, K[r, j], false)
            t2 = mul_red!(t2, b, K[r, i], false)
            K[r, j] = reduce!(t1 + t2)
            t1 = mul_red!(t1, u, K[r, i], false)
            t2 = mul_red!(t2, v, t, false)
            K[r, i] = reduce!(t1 + t2)
         end
      end
   end
   return nothing
end

function snf_kb!(S::MatrixElem{T}, U::MatrixElem{T}, K::MatrixElem{T}, with_trafo::Bool = false) where {T <: RingElement}
   m = nrows(S)
   n = ncols(S)
   l = min(m, n)
   i = 1
   t = base_ring(S)()
   t1 = base_ring(S)()
   t2 = base_ring(S)()
   while i <= l
      kb_clear_row!(S, K, i, with_trafo)
      hnf_kb!(S, U, with_trafo, i)
      c = i + 1
      while c <= n && is_zero_entry(S, i, c)
         c += 1
      end
      if c != n + 1
         continue
      end
      i+=1
   end
   for i = 1:l-1
      for j = i + 1:l
         if isone(S[i, i])
           break
         end
         if is_zero_entry(S, i, i) && is_zero_entry(S, j, j)
            continue
         end
         d, u, v = gcdx(S[i, i], S[j, j])
         if with_trafo
            q = -divexact(S[j, j], d)
            t1 = mul!(t1, q, v)
            for c = 1:m
               t = deepcopy(U[i, c])
               U[i, c] += U[j, c]
               t2 = mul_red!(t2, t1, U[j, c], false)
               U[j, c] += t2
               t2 = mul_red!(t2, t1, t, false)
               U[j, c] = reduce!(U[j, c] + t2)
            end
            q1 = -divexact(S[j, j], d)
            q2 = divexact(S[i, i], d)
            for r = 1:n
               t = deepcopy(K[r, i])
               t1 = mul_red!(t1, K[r, i], u, false)
               t2 = mul_red!(t2, K[r, j], v, false)
               K[r, i] = reduce!(t1 + t2)
               t1 = mul_red!(t1, t, q1, false)
               t2 = mul_red!(t2, K[r, j], q2, false)
               K[r, j] = reduce!(t1 + t2)
            end
         end
         S[j, j] = S[i, i]*divexact(S[j, j], d)
         S[i, i] = d
      end
   end
   return nothing
end

@doc raw"""
    snf(A::MatrixElem{T}) where {T <: RingElement}

Return the Smith normal form of $A$.
"""
function snf(a::MatrixElem{T}) where {T <: RingElement}
  return snf_kb(a)
end

@doc raw"""
    snf_with_transform(A)

Return the tuple $S, T, U$ consisting of the Smith normal form $S$ of $A$
together with invertible matrices $T$ and $U$ such that $TAU = S$.
"""
function snf_with_transform(a::MatrixElem{T}) where {T <: RingElement}
  return snf_kb_with_transform(a)
end

################################################################################
#
#   Popov Form
#
################################################################################

@doc raw"""
    is_weak_popov(P::MatrixElem{T}, rank::Int) where T <: PolyRingElem

Return `true` if $P$ is a matrix in weak Popov form of the given rank.
"""
function is_weak_popov(P::MatrixElem{T}, rank::Int) where T <: PolyRingElem
   zero_rows = 0
   pivots = zeros(ncols(P))
   for r = 1:nrows(P)
      p = find_pivot_popov(P, r)
      if P[r, p] == 0
         zero_rows += 1
         continue
      end
      # There is already a pivot in this column
      if pivots[p] != 0
         return false
      end
      pivots[p] = r
   end
   if zero_rows != nrows(P) - rank
      return false
   end
   return true
end

@doc raw"""
    is_popov(P::MatrixElem{T}, rank::Int) where T <: PolyRingElem

Return `true` if $P$ is a matrix in Popov form with the given rank.
"""
function is_popov(P::MatrixElem{T}, rank::Int) where T <: PolyRingElem
   zero_rows = 0
   for r = 1:nrows(P)
      p = find_pivot_popov(P, r)
      if P[r, p] != 0
         break
      end
      zero_rows += 1
   end
   # The zero rows must all be on top.
   if zero_rows != nrows(P) - rank
      return false
   end
   pivotscr = zeros(Int, ncols(P)) # pivotscr[i] == j means the pivot of column i is in row j
   pivots = zeros(Int, nrows(P)) # the other way round
   for r = zero_rows + 1:nrows(P)
      p = find_pivot_popov(P, r)
      if P[r, p] == 0
         return false
      end
      if pivotscr[p] != 0
         # There is already a pivot in this column
         return false
      end
      pivotscr[p] = r
      pivots[r] = p
   end
   for r = zero_rows + 1:nrows(P)
      p = pivots[r]
      f = P[r, p]
      if !isone(leading_coefficient(f))
         return false
      end
      for i = 1:nrows(P)
         i == r ? continue : nothing
         if degree(P[i, p]) >= degree(f)
            return false
         end
      end
      if r == nrows(P)
         break
      end
      g = P[r + 1, pivots[r + 1]]
      if degree(f) >= degree(g)
         if pivots[r] >= pivots[r + 1]
            return false
         end
      end
   end
   return true
end

@doc raw"""
    weak_popov(A::MatElem{T}) where {T <: PolyRingElem}

Return the weak Popov form of $A$.
"""
function weak_popov(A::MatElem{T}) where {T <: PolyRingElem}
   return _weak_popov(A, Val(false))
end

@doc raw"""
    weak_popov_with_transform(A::MatElem{T}) where {T <: PolyRingElem}

Compute a tuple $(P, U)$ where $P$ is the weak Popov form of $A$ and $U$
is a transformation matrix so that $P = UA$.
"""
function weak_popov_with_transform(A::MatElem{T}) where {T <: PolyRingElem}
   return _weak_popov(A, Val(true))
end

function _weak_popov(A::MatElem{T}, ::Val{with_transform} = Val(false)) where {T <: PolyRingElem, with_transform}
   P = deepcopy(A)
   m = nrows(P)
   W = similar(A, 0, 0)
   if with_transform
      U = identity_matrix(A, m)
      weak_popov!(P, W, U, false, true)
      return P, U
   else
      U = similar(A, 0, 0)
      weak_popov!(P, W, U, false, false)
      return P
   end
end

@doc raw"""
    extended_weak_popov(A::MatElem{T}, V::MatElem{T}) where {T <: PolyRingElem}

Compute the weak Popov form $P$ of $A$ by applying simple row transformations
on $A$ and a vector $W$ by applying the same transformations on the vector $V$.
Return the tuple $(P, W)$.
"""
function extended_weak_popov(A::MatElem{T}, V::MatElem{T}) where {T <: PolyRingElem}
   return _extended_weak_popov(A, V, Val(false))
end

@doc raw"""
    extended_weak_popov_with_transform(A::MatElem{T}, V::MatElem{T}) where {T <: PolyRingElem}

Compute the weak Popov form $P$ of $A$ by applying simple row transformations
on $A$, a vector $W$ by applying the same transformations on the vector $V$,
and a transformation matrix $U$ so that $P = UA$.
Return the tuple $(P, W, U)$.
"""
function extended_weak_popov_with_transform(A::MatElem{T}, V::MatElem{T}) where {T <: PolyRingElem}
   return _extended_weak_popov(A, V, Val(true))
end

function _extended_weak_popov(A::MatElem{T}, V::MatElem{T}, ::Val{with_transform} = Val(false)) where {T <: PolyRingElem, with_transform}
   @assert nrows(V) == nrows(A) && ncols(V) == 1
   P = deepcopy(A)
   W = deepcopy(V)
   m = nrows(P)
   if with_transform
      U = identity_matrix(A)
      weak_popov!(P, W, U, true, true)
      return P, W, U
   else
      U = similar(A, 0, 0)
      weak_popov!(P, W, U, true, false)
      return P, W
   end
end

function find_pivot_popov(P::MatElem{T}, r::Int, last_col::Int = 0) where {T <: PolyRingElem}
   last_col == 0 ? n = ncols(P) : n = last_col
   pivot = n
   for c = n-1:-1:1
      if degree(P[r,c]) > degree(P[r,pivot])
         pivot = c
      end
   end
   return pivot
end

function init_pivots_popov(P::MatElem{T}, last_row::Int = 0, last_col::Int = 0) where {T <: PolyRingElem}
   last_row == 0 ? m = nrows(P) : m = last_row
   last_col == 0 ? n = ncols(P) : n = last_col
   pivots = Vector{Vector{Int}}(undef, n)
   for i = 1:n
      pivots[i] = zeros(Int, 0)
   end
   # pivots[i] contains the indices of the rows in which the pivot element is in the ith column.
   for r = 1:m
      pivot = find_pivot_popov(P, r, last_col)
      !is_zero_entry(P, r, pivot) ? push!(pivots[pivot], r) : nothing
   end
   return pivots
end

function weak_popov!(P::MatElem{T}, W::MatElem{T}, U::MatElem{T}, extended::Bool = false,
                                       with_trafo::Bool = false, last_row::Int = 0, last_col::Int = 0) where {T <: PolyRingElem}
   pivots = init_pivots_popov(P, last_row, last_col)
   weak_popov_with_pivots!(P, W, U, pivots, extended, with_trafo, last_row, last_col)
   return nothing
end

#=
The weak Popov form is defined by T. Mulders and A. Storjohann in
"On lattice reduction for polynomial matrices"
=#
function weak_popov_with_pivots!(P::MatElem{T}, W::MatElem{T}, U::MatElem{T}, pivots::Vector{Vector{Int}},
                                                   extended::Bool = false, with_trafo::Bool = false, last_row::Int = 0, last_col::Int = 0) where {T <: PolyRingElem}
   last_row == 0 ? m = nrows(P) : m = last_row
   last_col == 0 ? n = ncols(P) : n = last_col
   @assert length(pivots) >= n

   t = base_ring(P)()
   change = true
   while change
      change = false
      for i = 1:n
         if length(pivots[i]) <= 1
            continue
         end
         change = true
         # Reduce with the pivot of minimal degree.
         # TODO: Remove the list comprehension as soon as indmin(f, A) works
         #pivotInd = indmin(degree(P[j, i]) for j in pivots[i])
         pivotInd = argmin([degree(P[j, i]) for j in pivots[i]])
         pivot = pivots[i][pivotInd]
         for j = 1:length(pivots[i])
            if j == pivotInd
               continue
            end
            q = -div(P[pivots[i][j], i], P[pivot, i])
            for c = 1:n
               t = mul!(t, q, P[pivot, c])
               P[pivots[i][j], c] = add!(P[pivots[i][j], c], t)
            end
            if with_trafo
               for c = 1:ncols(U)
                  t = mul!(t, q, U[pivot, c])
                  U[pivots[i][j], c] = add!(U[pivots[i][j], c], t)
               end
            end
            if extended
               t = mul!(t, q, W[pivot,1])
               W[pivots[i][j], 1] = add!(W[pivots[i][j], 1], t)
            end
         end
         old_pivots = pivots[i]
         pivots[i] = [pivot]
         for j = 1:length(old_pivots)
            if j == pivotInd
               continue
            end
            p = find_pivot_popov(P, old_pivots[j], last_col)
            if !is_zero_entry(P, old_pivots[j], p)
               push!(pivots[p], old_pivots[j])
            end
         end
      end
   end
   return nothing
end

@doc raw"""
    rank_profile_popov(A::MatElem{T}) where {T <: PolyRingElem}

Return an array of $r$ row indices such that these rows of $A$ are linearly
independent, where $r$ is the rank of $A$.
"""
function rank_profile_popov(A::MatElem{T}) where {T <: PolyRingElem}
   B = deepcopy(A)
   m = nrows(A)
   n = ncols(A)
   U = similar(A, 0, 0)
   V = U
   r = 0
   rank_profile = Vector{Int}(undef, 0)
   pivots = Vector{Vector{Int}}(undef, n)
   for i = 1:n
      pivots[i] = zeros(Int, 0)
   end
   p = find_pivot_popov(B, 1)
   if !is_zero_entry(B, 1, p)
      push!(pivots[p], 1)
      r = 1
      push!(rank_profile, 1)
   end
   for i = 2:m
      p = find_pivot_popov(B, i)
      !is_zero_entry(B, i, p) ? push!(pivots[p], i) : nothing
      weak_popov_with_pivots!(B, V, U, pivots, false, false, i)
      s = 0
      for j = 1:n
         # length(pivots[j]) is either 1 or 0, since B is in weak
         # Popov form.
         s += length(pivots[j])
      end
      if s != r
         push!(rank_profile, i)
         r = s
      end
   end
   return rank_profile
end

function det_popov(A::MatElem{T}) where {T <: PolyRingElem}
   nrows(A) != ncols(A) && error("Not a square matrix in det_popov.")
   B = deepcopy(A)
   n = ncols(B)
   R = base_ring(B)
   det = one(R)
   U = similar(A, 0, 0)
   V = U
   t = R()
   pivots1 = init_pivots_popov(B)
   weak_popov_with_pivots!(B, V, U, pivots1, false, false)
   pivots = zeros(Int, n)
   diag_elems = zeros(Int, n)
   for i = 1:n
      if isassigned(pivots1, i) && length(pivots1[i]) > 0
         pivots[i] = pivots1[i][1]
      else
         # If there is no pivot in the ith column, A has not full rank.
         return zero(R)
      end
   end
   for i = n-1:-1:1
      # "Remove" the column i+1 and compute a weak Popov Form of the
      # remaining matrix.
      r1 = pivots[i+1]
      c = find_pivot_popov(B, r1, i)
      # If the pivot B[r1, c] is zero then the row is zero.
      while !is_zero_entry(B, r1, c)
         r2 = pivots[c]
         if degree(B[r2, c]) > degree(B[r1,c])
            r1, r2 = r2, r1
            pivots[c] = r2
         end
         q = -div(B[r1, c], B[r2, c])
         for j = 1:i + 1
            t = mul!(t, q, B[r2, j])
            B[r1, j] = add!(B[r1, j], t)
         end
         c = find_pivot_popov(B, r1, i)
      end
      if is_zero_entry(B, r1, i+1)
         return zero(R)
      end
      diag_elems[i+1] = r1
      det = mul!(det, det, B[r1,i+1])
   end
   det = mul!(det, det, B[pivots[1],1])
   diag_elems[1] = pivots[1]
   number_of_swaps = 0
   # Adjust the sign of det by sorting the diagonal elements.
   for i = 1:n
      while diag_elems[i] != i
         r = diag_elems[i]
         diag_elems[i] = diag_elems[r]
         diag_elems[r] = r
         number_of_swaps += 1
      end
   end
   if number_of_swaps%2 == 1
      det = mul!(det, det, R(-1))
   end
   return det
end

@doc raw"""
    popov(A::MatElem{T}) where {T <: PolyRingElem}

Return the Popov form of $A$.
"""
function popov(A::MatElem{T}) where {T <: PolyRingElem}
   return _popov(A, Val(false))
end

@doc raw"""
    popov_with_transform(A::MatElem{T}) where {T <: PolyRingElem}

Compute a tuple $(P, U)$ where $P$ is the Popov form of $A$ and $U$
is a transformation matrix so that $P = UA$.
"""
function popov_with_transform(A::MatElem{T}) where {T <: PolyRingElem}
   return _popov(A, Val(true))
end

function _popov(A::MatElem{T}, ::Val{with_transform} = Val(false)) where {T <: PolyRingElem, with_transform}
   P = deepcopy(A)
   m = nrows(P)
   if with_transform
      U = identity_matrix(A, m)
      popov!(P, U, true)
      return P, U
   else
      U = similar(A, 0, 0)
      popov!(P, U, false)
      return P
   end
end

function asc_order_popov!(P::MatElem{T}, U::MatElem{T}, pivots::Vector{Vector{Int}}, with_trafo::Bool) where {T <: PolyRingElem}
   m = nrows(P)
   n = ncols(P)
   pivots2 = Vector{NTuple{3,Int}}(undef, m)
   for r = 1:m
      pivots2[r] = (r,n,-1)
   end
   for c = 1:n
      if length(pivots[c]) == 0
         continue
      end
      r = pivots[c][1]
      pivots2[r] = (r, c, degree(P[r,c]))
   end
   sort!(pivots2, lt = (x,y) -> ( x[3] < y[3] || ( x[3] == y[3] && x[2] <= y[2] ) ))
   row_nums = [ i for i = 1:m ]
   for r = 1:m
      if pivots2[r][3] != -1
         c = pivots2[r][2]
         pivots[c] = [r]
      end
      i = pivots2[r][1]
      r2 = row_nums[i]
      if r == r2
         continue
      end
      swap_rows!(P, r, r2)
      with_trafo ? swap_rows!(U, r, r2) : nothing
      j = findfirst(isequal(r), row_nums)
      row_nums[i] = r
      row_nums[j] = r2
   end
   return nothing
end

# Mulders, Storjohann: "On lattice reduction for polynomial matrices", Section 7
function popov!(P::MatElem{T}, U::MatElem{T}, with_trafo::Bool = false) where {T <: PolyRingElem}
   m = nrows(P)
   n = ncols(P)
   W = similar(U, 0, 0)
   pivots = init_pivots_popov(P)
   weak_popov_with_pivots!(P, W, U, pivots, false, with_trafo)
   asc_order_popov!(P, U, pivots, with_trafo)
   pivotColumns = zeros(Int, m)
   for c = 1:n
      if length(pivots[c]) == 0
         continue
      end
      # P is in weak Popov form, so any non-zero row contains exactly one pivot
      pivotColumns[pivots[c][1]] = c
   end

   t = base_ring(P)()
   for r = 1:m
      # Reduce row r assuming that rows 1, ..., r - 1 are reduced
      if iszero(pivotColumns[r])
         continue
      end
      r2 = 1
      while !iszero(r2)
         r2 = 0
         maxDegreeDiff = -1
         # Find the column c2 for which the difference
         # degree(P[r, c2]) - degree(P[r2, c2]) is maximal, where r2 is the row
         # in which the pivot is in column c2.
         for i = 1:r - 1
            ci = pivotColumns[i]
            if iszero(ci)
               continue
            end
            d = degree(P[r, ci]) - degree(P[i, ci])
            if d <= maxDegreeDiff
               continue
            end
            maxDegreeDiff = d
            r2 = i
         end
         iszero(r2) ? break : nothing
         c2 = pivotColumns[r2]
         q = -div(P[r, c2], P[r2, c2])
         for c = 1:n
            t = mul!(t, q, P[r2, c])
            P[r, c] = add!(P[r, c], t)
         end
         if with_trafo
            for c = 1:ncols(U)
               t = mul!(t, q, U[r2, c])
               U[r, c] = add!(U[r, c], t)
            end
         end
      end
   end
   # Make the pivots monic
   for i = 1:n
      if length(pivots[i]) == 0
         continue
      end
      r = pivots[i][1]
      cu = canonical_unit(P[r, i])
      if !isone(cu)
         for j = 1:n
            P[r, j] = divexact(P[r, j], cu)
         end
         if with_trafo
            for j = 1:ncols(U)
               U[r, j] = divexact(U[r, j], cu)
            end
         end
      end
   end
   return nothing
end

function hnf_via_popov(A::MatElem{T}) where {T <: PolyRingElem}
   return _hnf_via_popov(A, Val(false))
end

function hnf_via_popov_with_transform(A::MatElem{T}) where {T <: PolyRingElem}
   return _hnf_via_popov(A, Val(true))
end

function _hnf_via_popov(A::MatElem{T}, ::Val{with_transform} = Val(false)) where {T <: PolyRingElem, with_transform}
   H = deepcopy(A)
   m = nrows(H)
   if with_transform
      U = identity_matrix(A, m)
      hnf_via_popov!(H, U, true)
      return H, U
   else
      U = similar(A, 0, 0)
      hnf_via_popov!(H, U, false)
      return H
   end
end

function hnf_via_popov_reduce_row!(H::MatElem{T}, U::MatElem{T}, pivots_hermite::Vector{Int}, r::Int, with_trafo::Bool) where {T <: PolyRingElem}
   n = ncols(H)
   t = base_ring(H)()
   for c = 1:n
      if pivots_hermite[c] == 0
         continue
      end
      pivot = pivots_hermite[c]
      q = -div(H[r, c], H[pivot, c])
      for j = c:n
         t = mul!(t, q, H[pivot, j])
         H[r, j] = add!(H[r, j], t)
      end
      if with_trafo
         for j = 1:ncols(U)
            t = mul!(t, q, U[pivot, j])
            U[r, j] = add!(U[r, j], t)
         end
      end
   end
   return nothing
end

function hnf_via_popov_reduce_column!(H::MatElem{T}, U::MatElem{T}, pivots_hermite::Vector{Int}, c::Int, with_trafo::Bool) where {T <: PolyRingElem}
   m = nrows(H)
   n = ncols(H)
   t = base_ring(H)()
   r = pivots_hermite[c]
   for i = 1:m
      if i == r
         continue
      end
      if degree(H[i, c]) < degree(H[r, c])
         continue
      end
      q = -div(H[i, c], H[r, c])
      for j = 1:n
         t = mul!(t, q, H[r, j])
         H[i, j] = add!(H[i, j], t)
      end
      if with_trafo
         for j = 1:ncols(U)
            t = mul!(t, q, U[r, j])
            U[i, j] = add!(U[i, j], t)
         end
      end
   end
   return nothing
end

function hnf_via_popov!(H::MatElem{T}, U::MatElem{T}, with_trafo::Bool = false) where {T <: PolyRingElem}
   m = nrows(H)
   n = ncols(H)
   R = base_ring(H)
   W = similar(H, 0, 0)
   t = R()
   pivots = init_pivots_popov(H)
   weak_popov_with_pivots!(H, W, U, pivots, false, with_trafo)
   pivots_popov = zeros(Int, n)
   for j = 1:n
      if isassigned(pivots, j) && length(pivots[j]) > 0
         pivots_popov[j] = pivots[j][1]
      else
         error("The rank must be equal to the number of columns.")
      end
   end
   pivots_hermite = zeros(Int, n)
   for i = n-1:-1:1
      # "Remove" the column i+1 and compute a weak Popov Form of the
      # remaining matrix.
      r1 = pivots_popov[i + 1]
      c = find_pivot_popov(H, r1, i)
      new_pivot = true
      # If the pivot H[r1, c] is zero then the row is zero.
      while !is_zero_entry(H, r1, c)
         r2 = pivots_popov[c]
         if degree(H[r2, c]) > degree(H[r1,c])
            r1, r2 = r2, r1
            pivots_popov[c] = r2
         end
         q = -div(H[r1, c], H[r2, c])
         for j = 1:n
            t = mul!(t, q, H[r2, j])
            H[r1, j] = add!(H[r1, j], t)
         end
         if with_trafo
            for j = 1:ncols(U)
               t = mul!(t, q, U[r2, j])
               U[r1, j] = add!(U[r1, j], t)
            end
         end
         hnf_via_popov_reduce_row!(H, U, pivots_hermite, r1, with_trafo)
         c = find_pivot_popov(H, r1, i)
      end
      new_pivot ? nothing : continue
      pivots_hermite[i+1] = r1
      hnf_via_popov_reduce_column!(H, U, pivots_hermite, i+1, with_trafo)
      l = pivots_popov[i]
      hnf_via_popov_reduce_row!(H, U, pivots_hermite, l, with_trafo)
   end
   pivots_hermite[1] = pivots_popov[1]
   kb_sort_rows!(H, U, pivots_hermite, with_trafo)
   for c = 1:n
      kb_canonical_row!(H, U, pivots_hermite[c], c, with_trafo)
   end
   return nothing
end

###############################################################################
#
#   Transforms
#
###############################################################################

@doc raw"""
    similarity!(A::MatrixElem{T}, r::Int, d::T) where {T <: RingElement}

Applies a similarity transform to the $n\times n$ matrix $M$ in-place. Let
$P$ be the $n\times n$ identity matrix that has had all zero entries of row
$r$ replaced with $d$, then the transform applied is equivalent to
$M = P^{-1}MP$. We require $M$ to be a square matrix. A similarity transform
preserves the minimal and characteristic polynomials of a matrix.

# Examples

```jldoctest
julia> R, = residue_ring(ZZ, 7);

julia> S = matrix_space(R, 4, 4)
Matrix space of 4 rows and 4 columns
  over residue ring of integers modulo 7

julia> M = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);
              R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])
[1   2   4   3]
[2   5   1   0]
[6   1   3   2]
[1   1   3   5]

julia> similarity!(M, 1, R(3))

```
"""
function similarity!(A::MatrixElem{T}, r::Int, d::T) where {T <: RingElement}
   n = nrows(A)
   t = base_ring(A)()
   for i = 1:n
      for j = 1:r - 1
         t = mul!(t, A[i, r], d)
         A[i, j] = add!(A[i, j], t)
      end
      for j = r + 1:n
         t = mul!(t, A[i, r], d)
         A[i, j] = add!(A[i, j], t)
      end
   end
   d = -d
   for i = 1:n
      for j = 1:r - 1
         A[r, i] = addmul_delayed_reduction!(A[r, i], A[j, i], d, t)
      end
      for j = r + 1:n
         A[r, i] = addmul_delayed_reduction!(A[r, i], A[j, i], d, t)
      end
      A[r, i] = reduce!(A[r, i])
   end
end

###############################################################################
#
#   Row & column permutations
#
###############################################################################

@doc raw"""
    swap_rows(a::MatrixElem{T}, i::Int, j::Int) where T <: NCRingElement

Return a matrix $b$ with the entries of $a$, where the $i$th and $j$th
row are swapped.

# Examples
```jldoctest
julia> M = identity_matrix(ZZ, 3)
[1   0   0]
[0   1   0]
[0   0   1]

julia> swap_rows(M, 1, 2)
[0   1   0]
[1   0   0]
[0   0   1]

julia> M  # was not modified
[1   0   0]
[0   1   0]
[0   0   1]
```
"""
function swap_rows(a::MatrixElem{T}, i::Int, j::Int) where T <: NCRingElement
   (1 <= i <= nrows(a) && 1 <= j <= nrows(a)) || throw(BoundsError())
   b = deepcopy(a)
   swap_rows!(b, i, j)
   return b
end

@doc raw"""
    swap_rows!(a::MatrixElem{T}, i::Int, j::Int) where T <: NCRingElement

Swap the $i$th and $j$th row of $a$ in place. The function returns the mutated
matrix (since matrices are assumed to be mutable in AbstractAlgebra.jl).

# Examples
```jldoctest
julia> M = identity_matrix(ZZ, 3)
[1   0   0]
[0   1   0]
[0   0   1]

julia> swap_rows!(M, 1, 2)
[0   1   0]
[1   0   0]
[0   0   1]

julia> M  # was modified
[0   1   0]
[1   0   0]
[0   0   1]
```
"""
function swap_rows!(a::MatrixElem{T}, i::Int, j::Int) where T <: NCRingElement
   if i != j
      for k = 1:ncols(a)
         a[i, k], a[j, k] = a[j, k], a[i, k]
      end
   end
   return a
end

@doc raw"""
    swap_cols(a::MatrixElem{T}, i::Int, j::Int) where T <: NCRingElement

Return a matrix $b$ with the entries of $a$, where the $i$th and $j$th
row are swapped.
"""
function swap_cols(a::MatrixElem{T}, i::Int, j::Int) where T <: NCRingElement
   (1 <= i <= ncols(a) && 1 <= j <= ncols(a)) || throw(BoundsError())
   b = deepcopy(a)
   swap_cols!(b, i, j)
   return b
end

@doc raw"""
    swap_cols!(a::MatrixElem{T}, i::Int, j::Int) where T <: NCRingElement

Swap the $i$th and $j$th column of $a$ in place. The function returns the mutated
matrix (since matrices are assumed to be mutable in AbstractAlgebra.jl).
"""
function swap_cols!(a::MatrixElem{T}, i::Int, j::Int) where T <: NCRingElement
   if i != j
      for k = 1:nrows(a)
         a[k, i], a[k, j] = a[k, j], a[k, i]
      end
   end
   return a
end

@doc raw"""
    reverse_rows!(a::MatrixElem{T}) where T <: NCRingElement

Swap the $i$th and $r - i$th row of $a$ for $1 \leq i \leq r/2$,
where $r$ is the number of rows of $a$.
"""
function reverse_rows!(a::MatrixElem{T}) where T <: NCRingElement
   k = div(nrows(a), 2)
   for i in 1:k
      swap_rows!(a, i, nrows(a) - i + 1)
   end
   return a
end

@doc raw"""
    reverse_rows(a::MatrixElem{T}) where T <: NCRingElement

Return a matrix $b$ with the entries of $a$, where the $i$th and $r - i$th
row is swapped for $1 \leq i \leq r/2$. Here $r$ is the number of rows of
$a$.
"""
function reverse_rows(a::MatrixElem{T}) where T <: NCRingElement
   b = deepcopy(a)
   return reverse_rows!(b)
end

@doc raw"""
    reverse_cols!(a::MatrixElem{T}) where T <: NCRingElement

Swap the $i$th and $r - i$th column of $a$ for $1 \leq i \leq c/2$,
where $c$ is the number of columns of $a$.
"""
function reverse_cols!(a::MatrixElem{T}) where T <: NCRingElement
   k = div(ncols(a), 2)
   for i in 1:k
      swap_cols!(a, i, ncols(a) - i + 1)
   end
   return a
end

@doc raw"""
    reverse_cols(a::MatrixElem{T}) where T <: NCRingElement

Return a matrix $b$ with the entries of $a$, where the $i$th and $r - i$th
column is swapped for $1 \leq i \leq c/2$. Here $c$ is the number of columns
of$a$.
"""
function reverse_cols(a::MatrixElem{T}) where T <: NCRingElement
   b = deepcopy(a)
   return reverse_cols!(b)
end

################################################################################
#
#  Elementary row/column transformations
#
################################################################################

@doc raw"""
    add_column!(a::MatrixElem{T}, s::RingElement, i::Int, j::Int, rows = 1:nrows(a)) where T <: RingElement

Add $s$ times the $i$-th row to the $j$-th row of $a$.

By default, the transformation is applied to all rows of $a$. This can be
changed using the optional `rows` argument.
"""
function add_column!(a::MatrixElem{T}, s::RingElement, i::Int, j::Int, rows = 1:nrows(a)) where T <: RingElement
   v = base_ring(a)(s)
   nc = ncols(a)
   !_checkbounds(nc, i) && error("Column index ($i) must be between 1 and $nc")
   !_checkbounds(nc, j) && error("Column index ($j) must be between 1 and $nc")
   temp = base_ring(a)()
   for r in rows
      temp = mul!(temp, v, a[r, i])
      a[r, j] += temp # cannot mutate matrix entries
   end
   return a
end

@doc raw"""
    add_column(a::MatrixElem{T}, s::RingElement, i::Int, j::Int, rows = 1:nrows(a)) where T <: RingElement

Create a copy of $a$ and add $s$ times the $i$-th row to the $j$-th row of $a$.

By default, the transformation is applied to all rows of $a$. This can be
changed using the optional `rows` argument.

"""
function add_column(a::MatrixElem{T}, s::RingElement, i::Int, j::Int, rows = 1:nrows(a)) where T <: RingElement
   b = deepcopy(a)
   return add_column!(b, s, i, j, rows)
end

@doc raw"""
    add_row!(a::MatrixElem{T}, s::RingElement, i::Int, j::Int, cols = 1:ncols(a)) where T <: RingElement

Add $s$ times the $i$-th row to the $j$-th row of $a$.

By default, the transformation is applied to all columns of $a$. This can be
changed using the optional `cols` argument.
"""
function add_row!(a::MatrixElem{T}, s::RingElement, i::Int, j::Int, cols = 1:ncols(a)) where T <: RingElement
   v = base_ring(a)(s)
   nr = nrows(a)
   !_checkbounds(nr, i) && error("Row index ($i) must be between 1 and $nr")
   !_checkbounds(nr, j) && error("Row index ($j) must be between 1 and $nr")
   temp = base_ring(a)()
   for c in cols
      temp = mul!(temp, v, a[i, c])
      a[j, c] += temp # cannot mutate matrix entries
   end
   return a
end

@doc raw"""
    add_row(a::MatrixElem{T}, s::RingElement, i::Int, j::Int, cols = 1:ncols(a)) where T <: RingElement

Create a copy of $a$ and add $s$ times the $i$-th row to the $j$-th row of $a$.

By default, the transformation is applied to all columns of $a$. This can be
changed using the optional `cols` argument.
"""
function add_row(a::MatrixElem{T}, s::RingElement, i::Int, j::Int, cols = 1:ncols(a)) where T <: RingElement
   b = deepcopy(a)
   return add_row!(b, s, i, j, cols)
end

# Multiply column

@doc raw"""
    multiply_column!(a::MatrixElem{T}, s::RingElement, i::Int, rows = 1:nrows(a)) where T <: RingElement

Multiply the $i$th column of $a$ with $s$.

By default, the transformation is applied to all rows of $a$. This can be
changed using the optional `rows` argument.
"""
function multiply_column!(a::MatrixElem{T}, s::RingElement, i::Int, rows = 1:nrows(a)) where T <: RingElement
   c = base_ring(a)(s)
   nc = ncols(a)
   !_checkbounds(nc, i) && error("Row index ($i) must be between 1 and $nc")
   temp = base_ring(a)()
   for r in rows
      a[r, i] = c*a[r, i] # cannot mutate matrix entries
   end
   return a
end

@doc raw"""
    multiply_column(a::MatrixElem{T}, s::RingElement, i::Int, rows = 1:nrows(a)) where T <: RingElement

Create a copy of $a$ and multiply the $i$th column of $a$ with $s$.

By default, the transformation is applied to all rows of $a$. This can be
changed using the optional `rows` argument.
"""
function multiply_column(a::MatrixElem{T}, s::RingElement, i::Int, rows = 1:nrows(a)) where T <: RingElement
   b = deepcopy(a)
   return multiply_column!(b, s, i, rows)
end

# Multiply row

@doc raw"""
    multiply_row!(a::MatrixElem{T}, s::RingElement, i::Int, cols = 1:ncols(a)) where T <: RingElement

Multiply the $i$th row of $a$ with $s$.

By default, the transformation is applied to all columns of $a$. This can be
changed using the optional `cols` argument.
"""
function multiply_row!(a::MatrixElem{T}, s::RingElement, i::Int, cols = 1:ncols(a)) where T <: RingElement
   c = base_ring(a)(s)
   nr = nrows(a)
   !_checkbounds(nr, i) && error("Row index ($i) must be between 1 and $nr")
   temp = base_ring(a)()
   for r in cols
      a[i, r] = c*a[i, r] # cannot mutate matrix entries
   end
   return a
end

@doc raw"""
    multiply_row(a::MatrixElem{T}, s::RingElement, i::Int, cols = 1:ncols(a)) where T <: RingElement

Create a copy of $a$ and multiply  the $i$th row of $a$ with $s$.

By default, the transformation is applied to all columns of $a$. This can be
changed using the optional `cols` argument.
"""
function multiply_row(a::MatrixElem{T}, s::RingElement, i::Int, cols = 1:ncols(a)) where T <: RingElement
   b = deepcopy(a)
   return multiply_row!(b, s, i, cols)
end

###############################################################################
#
#   Concatenation
#
###############################################################################

@doc raw"""
    vcat(A::MatrixElem{T}...) where T <: NCRingElement -> MatrixElem

Return the horizontal concatenation of the matrices $A$.
All component matrices need to have the same base ring and number of columns.
"""
function Base.vcat(A::MatrixElem...)
  # We don't add a type parameter T <: NCRingElement, so that this function is
  # called for e.g. vcat(QQ[1 0; 0 1], ZZ[1 0; 0 1]) and ERRORS instead of
  # producing an array of the arguments.
  return _vcat(A)
end

# this leads to an ambiguity when calling `reduce(hcat, Union{}[])`, but we don't have a better solution right now
Base.reduce(::typeof(vcat), A::AbstractVector{<:MatrixElem}) = _vcat(A)

function _vcat(A)
  if length(A) == 0
    error("Number of matrices to concatenate must be positive")
  end

  if any(x -> ncols(x) != ncols(A[1]), A)
    error("Matrices must have the same number of columns")
  end

  if any(x -> base_ring(x) !== base_ring(A[1]), A)
    error("Matrices must have the same base ring")
  end

  M = similar(A[1], sum(nrows, A), ncols(A[1]))
  i = 1
  for j in 1:length(A)
    M[i:i + nrows(A[j]) - 1, :] = A[j]
    i += nrows(A[j])
  end
  return M
end

@doc raw"""
    hcat(A::MatrixElem{T}...) where T <: NCRingElement -> MatrixElem

Return the horizontal concatenating of the matrices $A$.
All component matrices need to have the same base ring and number of rows.
"""
function Base.hcat(A::MatrixElem...)
  # We don't add a type parameter T <: NCRingElement, so that this function is
  # called for e.g. vcat(QQ[1 0; 0 1], ZZ[1 0; 0 1]) and ERRORS instead of
  # producing an array of the arguments.
  return _hcat(A)
end

# this leads to an ambiguity when calling `reduce(hcat, Union{}[])`, but we don't have a better solution right now
Base.reduce(::typeof(hcat), A::AbstractVector{<:MatrixElem}) = _hcat(A)

function _hcat(A)
  if length(A) == 0
    error("Number of matrices to concatenate must be positive")
  end

  if any(x -> nrows(x) != nrows(A[1]), A)
    error("Matrices must have the same number of rows")
  end

  if any(x -> base_ring(x) !== base_ring(A[1]), A)
    error("Matrices must have the same base ring")
  end

  M = similar(A[1], nrows(A[1]), sum(ncols, A))
  i = 1
  for j in 1:length(A)
    M[:, i:i + ncols(A[j]) - 1] = A[j]
    i += ncols(A[j])
  end
  return M
end

function Base.cat(A::MatrixElem, As::MatrixElem...; dims)
  @assert dims == (1,2) || isa(dims, Int)

  if isa(dims, Int)
    if dims == 1
      return hcat(A, As...)
    elseif dims == 2
      return vcat(A, As...)
    else
      error("dims must be 1, 2, or (1,2)")
    end
  end

  X = hcat(A, zero(A, nrows(A), sum(Int[ncols(As[j]) for j=1:length(As)])))
  for i in 1:length(As)
    X = vcat(X, hcat(zero(A, nrows(As[i]), ncols(A) + sum(Int[ncols(As[j]) for j=1:i-1])), As[i], zero(A, nrows(As[i]), sum(Int[ncols(As[j]) for j in i+1:length(As)]))))
  end
  return X
end

function Base.hvcat(rows::Tuple{Vararg{Int}}, A::MatrixElem...)
  if any(x -> base_ring(x) !== base_ring(A[1]), A)
    error("Matrices must have the same base ring")
  end

  nr = 0
  k = 1
  for i in 1:length(rows)
    nr += nrows(A[k])
    k += rows[i]
  end

  nc = sum(ncols(A[i]) for i in 1:rows[1])

  M = similar(A[1], nr, nc)
  mat_offset = 0
  row_offset = 0
  for j in 1:length(rows)
    s = 0
    for i in 1:rows[j]
      N = A[mat_offset + i]
      M[row_offset + 1:row_offset + nrows(N), s + 1:s + ncols(N)] = N
      s += ncols(N)
    end
    row_offset += nrows(A[1 + mat_offset])
    mat_offset += rows[j]
  end

  return M
end

###############################################################################
#
#   Change Base Ring
#
###############################################################################

# like change_base_ring, but without initializing the entries
# this function exists until a better API is implemented
_change_base_ring(R::NCRing, a::MatElem) = dense_matrix_type(R)(R, undef, nrows(a), ncols(a))
_change_base_ring(R::NCRing, a::MatRingElem) = matrix_ring(R, nrows(a))()

@doc raw"""
    change_base_ring(R::NCRing, M::MatrixElem{T}) where T <: NCRingElement

Return the matrix obtained by coercing each entry into `R`.
"""
function change_base_ring(R::NCRing, M::MatrixElem{T}) where T <: NCRingElement
   N = _change_base_ring(R, M)
   for i = 1:nrows(M), j = 1:ncols(M)
      N[i,j] = R(M[i,j])
   end
   return N
end

###############################################################################
#
#   Map
#
###############################################################################

@doc raw"""
    map_entries!(f, dst::MatrixElem{T}, src::MatrixElem{U}) where {T <: NCRingElement, U <: NCRingElement}

Like `map_entries`, but stores the result in `dst` rather than a new matrix.
"""
function map_entries!(f::S, dst::MatrixElem{T}, src::MatrixElem{U}) where {S, T <: NCRingElement, U <: NCRingElement}
   for i = 1:nrows(src), j = 1:ncols(src)
      dst[i, j] = f(src[i, j])
   end
   dst
end

@doc raw"""
    map!(f, dst::MatrixElem{T}, src::MatrixElem{U}) where {T <: NCRingElement, U <: NCRingElement}

Like `map`, but stores the result in `dst` rather than a new matrix.
This is equivalent to `map_entries!(f, dst, src)`.
"""
Base.map!(f::S, dst::MatrixElem{T}, src::MatrixElem{U}) where {S, T <: NCRingElement, U <: NCRingElement} = map_entries!(f, dst, src)

@doc raw"""
    map_entries(f, a::MatrixElem{T}) where T <: NCRingElement

Transform matrix `a` by applying `f` on each element.
"""
function map_entries(f::S, a::MatrixElem{T}) where {S, T <: NCRingElement}
   isempty(a) && return _change_base_ring(parent(f(zero(base_ring(a)))), a)
   b11 = f(a[1, 1])
   b = _change_base_ring(parent(b11), a)
   b[1, 1] = b11
   for i = 1:nrows(a), j = 1:ncols(a)
      i == j == 1 && continue
      b[i, j] = f(a[i, j])
   end
   b
end

@doc raw"""
    map(f, a::MatrixElem{T}) where T <: NCRingElement

Transform matrix `a` by applying `f` on each element.
This is equivalent to `map_entries(f, a)`.
"""
Base.map(f::S, a::MatrixElem{T}) where {S, T <: NCRingElement} = map_entries(f, a)

###############################################################################
#
#   Random generation
#
###############################################################################

RandomExtensions.maketype(S::MatSpace, _) = elem_type(S)

function RandomExtensions.make(S::MatSpace, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1]) # forward to default Make constructor
   else
      Make(S, make(R, vs...))
   end
end

# Sampler for a MatSpace not needing arguments (passed via make)
# this allows to obtain the Sampler in simple cases without having to know about make
# (when one can do `rand(M)`, one can expect to be able to do `rand(Sampler(rng, M))`)
Random.Sampler(::Type{RNG}, S::MatSpace, n::Random.Repetition
               ) where {RNG<:AbstractRNG} =
   Random.Sampler(RNG, make(S), n)

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{<:MatElem,
                                         <:MatSpace}})
   S, v = sp[][1:end]
   M = S()
   R = base_ring(S)
   for i = 1:nrows(M)
      for j = 1:ncols(M)
         M[i, j] = rand(rng, v)
      end
   end
   return M
end

rand(rng::AbstractRNG, S::MatSpace, v...) = rand(rng, make(S, v...))

rand(S::MatSpace, v...) = rand(Random.default_rng(), S, v...)

# resolve ambiguities
rand(rng::AbstractRNG, S::MatSpace, dims::Integer...) =
   rand(rng, make(S), dims...)

rand(S::MatSpace, dims::Integer...) = rand(Random.default_rng(), S, dims...)

function randmat_triu(rng::AbstractRNG, S::MatSpace, v...)
   M = S()
   R = base_ring(S)
   for i = 1:nrows(M)
      for j = 1:i - 1
         M[i, j] = R()
      end
      for j = i:ncols(M)
         M[i, j] = rand(rng, R, v...)
      end
      while is_zero_entry(M, i, i)
         M[i, i] = rand(rng, R, v...)
      end
   end
   return M
end

randmat_triu(S::MatSpace, v...) = randmat_triu(Random.default_rng(), S, v...)

function randmat_with_rank(rng::AbstractRNG, S::MatSpace{T}, rank::Int, v...) where {T <: RingElement}
   if !is_domain_type(T) && !(T <: ResElem)
      error("Not implemented")
   end
   M = S()
   R = base_ring(S)
   for i = 1:rank
      for j = 1:i - 1
         M[i, j] = R()
      end
      M[i, i] = rand(rng, R, v...)
      while is_zero_entry(M, i, i)
         M[i, i] = rand(rng, R, v...)
      end
      for j = i + 1:ncols(M)
         M[i, j] = rand(rng, R, v...)
      end
   end
   for i = rank + 1:nrows(M)
      for j = 1:ncols(M)
         M[i, j] = R()
      end
   end
   m = nrows(M)
   if m > 1
      for i = 1:4*m
         r1 = rand(rng, 1:m)
         r2 = rand(rng, 1:m - 1)
         r2 = r2 >= r1 ? r2 + 1 : r2
         d = rand(rng, -5:5)
         for j = 1:ncols(M)
            M[r1, j] = M[r1, j] + d*M[r2, j]
         end
      end
   end
   return M
end

randmat_with_rank(S::MatSpace{T}, rank::Int, v...) where {T <: RingElement} =
   randmat_with_rank(Random.default_rng(), S, rank, v...)

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(S::MatSpace)
  R = base_ring(S)
  return S(elem_type(R)[ConformanceTests.generate_element(R) for i in 1:nrows(S), j in 1:ncols(S)])
end

################################################################################
#
#   Matrix constructors
#
################################################################################

@doc raw"""
    matrix(R::Ring, arr::AbstractMatrix{T}) where {T}

Constructs the matrix over $R$ with entries as in `arr`.

# Examples

```jldoctest
julia> matrix(GF(3), [1 2 ; 3 4])
[1   2]
[0   1]

julia> using LinearAlgebra ; matrix(GF(5), I(2))
[1   0]
[0   1]
```
"""
function matrix(R::NCRing, arr::AbstractMatrix{T}) where {T}
   Base.require_one_based_indexing(arr)
   if elem_type(R) === T && all(e -> parent(e) === R, arr)
      z = Generic.MatSpaceElem{elem_type(R)}(R, arr)
      return z
   else
      mat = (arr isa Matrix{T}) ? arr : convert(Matrix{T}, arr)
      arr_coerce = convert(Matrix{elem_type(R)}, map(R, mat))::Matrix{elem_type(R)}
      return matrix(R, arr_coerce)
   end
end

function matrix(R::NCRing, arr::MatElem)
   return map_entries(R, arr)
end

function matrix(R::NCRing, arr::MatRingElem)
   M = dense_matrix_type(R)(R, undef, nrows(arr), ncols(arr))
   for i in 1:nrows(arr), j in 1:ncols(arr)
      M[i, j] = arr[i, j]
   end
   return M
end

function matrix(mat::MatrixElem{T}) where {T<:NCRingElement}
   return matrix(base_ring(mat), mat)
end

function matrix(arr::AbstractMatrix{T}) where {T<:NCRingElement}
   Base.require_one_based_indexing(arr)
   r, c = size(arr)
   (r < 0 || c < 0) && error("Array must be non-empty")
   R = parent(arr[1, 1])
   all(e -> parent(e) === R, arr) || error("Non-compatible elements")
   return matrix(R, arr)
end

function matrix(arr::AbstractVector{T}) where {T<:NCRingElement}
   return matrix(reshape(arr, length(arr), 1))
end

function matrix(arr::AbstractVector{<:AbstractVector{T}}) where {T<:NCRingElement}
    return matrix(permutedims(reduce(hcat, arr)))
end

function matrix(R::NCRing, arr::AbstractVector{<:AbstractVector})
   return matrix(R, permutedims(reduce(hcat, arr)))
end

@doc raw"""
    matrix(R::Ring, r::Int, c::Int, arr::AbstractVector{T}) where {T}

Constructs the $r \times c$ matrix over $R$, where the entries are taken
row-wise from `arr`.
"""
function matrix(R::NCRing, r::Int, c::Int, arr::AbstractVecOrMat{T}) where T
   _check_dim(r, c, arr)
   ndims(arr) == 2 && return matrix(R, arr)
   if elem_type(R) === T && all(e -> parent(e) === R, arr)
     z = Generic.MatSpaceElem{elem_type(R)}(R, r, c, arr)
     return z
   else
     arr_coerce = convert(Vector{elem_type(R)}, map(R, arr))::Vector{elem_type(R)}
     return matrix(R, r, c, arr_coerce)
   end
end

################################################################################
#
#   Zero matrix
#
################################################################################

@doc raw"""
    zero_matrix(R::Ring, r::Int, c::Int)

Return the $r \times c$ zero matrix over $R$.
"""
function zero_matrix(R::NCRing, r::Int, c::Int)
  (r < 0 || c < 0) && error("Dimensions must be non-negative")
  mat = dense_matrix_type(R)(R, undef, r, c)
  return is_zero_initialized(mat) ? mat : zero!(mat)
end

zero_matrix(::Type{MatElem}, R::Ring, n::Int, m::Int) = zero_matrix(R, n, m)

################################################################################
#
#   Ones matrix
#
################################################################################

@doc raw"""
    ones_matrix(R::Ring, r::Int, c::Int)

Return the $r \times c$ ones matrix over $R$.
"""
function ones_matrix(R::NCRing, r::Int, c::Int)
   z = dense_matrix_type(R)(R, undef, r, c)
   for i in 1:r, j in 1:c
      z[i, j] = one(R)
   end
   return z
end

################################################################################
#
#   Identity matrix
#
################################################################################

@doc raw"""
    identity_matrix(R::NCRing, n::Int)

Return the $n \times n$ identity matrix over $R$.
"""
identity_matrix(R::NCRing, n::Int) = diagonal_matrix(one(R), n)

@doc raw"""
    identity_matrix(M::MatElem{T}) where T <: NCRingElement

Construct the identity matrix in the same matrix space as `M`, i.e.
with ones down the diagonal and zeroes elsewhere. `M` must be square.
This is an alias for `one(M)`.
"""
function identity_matrix(M::MatElem{T}) where T <: NCRingElement
   identity_matrix(check_square(M), nrows(M))
end

function identity_matrix(M::MatElem{T}, n::Int) where T <: NCRingElement
   z = zero(M, n, n)
   R = base_ring(M)
   for i = 1:n
      z[i, i] = one(R)
   end
   z
end

identity_matrix(::Type{MatElem}, R::Ring, n::Int) = identity_matrix(R, n)

################################################################################
#
#   Scalar matrix
#
################################################################################

@doc raw"""
    scalar_matrix(R::NCRing, n::Int, a::NCRingElement)
    scalar_matrix(n::Int, a::NCRingElement)

Return the $n \times n$ matrix over `R` with `a` along the main diagonal and
zeroes elsewhere. If `R` is not specified, it defaults to `parent(a)`.
"""
function scalar_matrix(R::NCRing, n::Int, a::NCRingElement)
   return diagonal_matrix(R(a), n)
end

function scalar_matrix(n::Int, a::NCRingElement)
   return diagonal_matrix(a, n)
end

################################################################################
#
#   Diagonal matrix
#
################################################################################

@doc raw"""
    diagonal_matrix(x::NCRingElement, m::Int, [n::Int])

Return the $m \times n$ matrix over $R$ with `x` along the main diagonal and
zeroes elsewhere. If `n` is not specified, it defaults to `m`.

# Examples
```jldoctest
julia> diagonal_matrix(ZZ(2), 2, 3)
[2   0   0]
[0   2   0]

julia> diagonal_matrix(QQ(-1), 3)
[-1//1    0//1    0//1]
[ 0//1   -1//1    0//1]
[ 0//1    0//1   -1//1]
```
"""
function diagonal_matrix(x::NCRingElement, m::Int, n::Int)
   z = zero_matrix(parent(x), m, n)
   for i in 1:min(m, n)
      z[i, i] = x
   end
   return z
end

diagonal_matrix(x::NCRingElement, m::Int) = diagonal_matrix(x, m, m)

@doc raw"""
    diagonal_matrix(x::T...) where T <: NCRingElement -> MatElem{T}
    diagonal_matrix(x::AbstractVector{T}) where T <: NCRingElement -> MatElem{T}
    diagonal_matrix(R::NCRing, x::AbstractVector{T}) where T <: NCRingElement -> MatElem{T}

Returns a diagonal matrix whose diagonal entries are the elements of $x$.
If a ring $R$ is given then it is used a parent for the entries of the created
matrix. Otherwise the parent is inferred from the vector $x$.

# Examples

```jldoctest
julia> diagonal_matrix(ZZ(1), ZZ(2))
[1   0]
[0   2]

julia> diagonal_matrix([ZZ(3), ZZ(4)])
[3   0]
[0   4]

julia> diagonal_matrix(ZZ, [5, 6])
[5   0]
[0   6]
```
"""
function diagonal_matrix(R::NCRing, x::AbstractVector{<:NCRingElement})
    Base.require_one_based_indexing(x)
    x = R.(x)
    M = zero_matrix(R, length(x), length(x))
    for i = 1:length(x)
        M[i, i] = x[i]
    end
    return M
end

function diagonal_matrix(x::T, xs::T...) where {T<:NCRingElement}
    return diagonal_matrix([x, xs...])
end

function diagonal_matrix(x::AbstractVector{<:NCRingElement})
   @req !isempty(x) "Cannot infer base ring from empty vector; consider passing the desired base ring as first argument to `diagonal_matrix`"
   return diagonal_matrix(parent(first(x)), x)
end

@doc raw"""
    diagonal_matrix(V::Vector{T}) where T <: MatElem -> MatElem

Returns a block diagonal matrix whose diagonal blocks are the matrices in $x$.
"""
function diagonal_matrix(V::Vector{T}) where {T<:MatElem}
    return block_diagonal_matrix(V)
end

function diagonal_matrix(x::T, xs::T...) where {T<:MatElem}
    return block_diagonal_matrix([x, xs...])
end

function diagonal_matrix(R::NCRing, V::Vector{<:MatElem})
    if length(V) == 0
        return zero_matrix(R, 0, 0)
    else
        return block_diagonal_matrix(map(x -> change_base_ring(R, x), V))
    end
end


###############################################################################
#
#   Lower triangular matrix
#
###############################################################################

@doc raw"""
    lower_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}

Return the $n$ by $n$ matrix whose entries on and below the main diagonal are
the elements of `L`, and which has zeroes elsewhere.
The value of $n$ is determined by the condition that `L` has length
$n(n+1)/2$.

An exception is thrown if there is no integer $n$ with this property.

# Examples
```jldoctest
julia> lower_triangular_matrix([1, 2, 3])
[1   0]
[2   3]
```
"""
function lower_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
   l = length(L)
   @req l > 0 "Input vector must be nonempty"
   (flag, m) = is_square_with_sqrt(8*l+1)
   @req flag "Input vector of invalid length"
   n = div(m-1, 2)
   R = parent(L[1])
   M = zero_matrix(R, n, n)
   pos = 1
   for i in 1:n, j in 1:i
      M[i,j] = L[pos]
      pos += 1
   end
   return M
end

###############################################################################
#
#   Upper triangular matrix
#
###############################################################################

@doc raw"""
    upper_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}

Return the $n$ by $n$ matrix whose entries on and above the main diagonal are
the elements of `L`, and which has zeroes elsewhere.
The value of $n$ is determined by the condition that `L` has length
$n(n+1)/2$.

An exception is thrown if there is no integer $n$ with this property.

# Examples
```jldoctest
julia> upper_triangular_matrix([1, 2, 3])
[1   2]
[0   3]
```
"""
function upper_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
   l = length(L)
   @req l > 0 "Input vector must be nonempty"
   (flag, m) = is_square_with_sqrt(8*l+1)
   @req flag "Input vector of invalid length"
   n = div(m-1, 2)
   R = parent(L[1])
   M = zero_matrix(R, n, n)
   pos = 1
   for i in 1:n, j in i:n
      M[i,j] = L[pos]
      pos += 1
   end
   return M
end

###############################################################################
#
#   Strictly lower triangular matrix
#
###############################################################################

@doc raw"""
    strictly_lower_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}

Return the $n$ by $n$ matrix whose entries below the main diagonal are
the elements of `L`, and which has zeroes elsewhere.
The value of $n$ is determined by the condition that `L` has length
$(n-1)n/2$.

An exception is thrown if there is no integer $n$ with this property.

# Examples
```jldoctest
julia> strictly_lower_triangular_matrix([1, 2, 3])
[0   0   0]
[1   0   0]
[2   3   0]
```
"""
function strictly_lower_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
   l = length(L)
   @req l > 0 "Input vector must be nonempty"
   (flag, m) = is_square_with_sqrt(8*l+1)
   @req flag "Input vector of invalid length"
   n = div(m+1, 2)
   R = parent(L[1])
   M = zero_matrix(R, n, n)
   pos = 1
   for i in 2:n, j in 1:(i-1)
      M[i,j] = L[pos]
      pos += 1
   end
   return M
end

###############################################################################
#
#   Strictly upper triangular matrix
#
###############################################################################

@doc raw"""
    strictly_upper_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}

Return the $n$ by $n$ matrix whose entries above the main diagonal are
the elements of `L`, and which has zeroes elsewhere.
The value of $n$ is determined by the condition that `L` has length
$(n-1)n/2$.

An exception is thrown if there is no integer $n$ with this property.

# Examples
```jldoctest
julia> strictly_upper_triangular_matrix([1, 2, 3])
[0   1   2]
[0   0   3]
[0   0   0]
```
"""
function strictly_upper_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
   l = length(L)
   @req l > 0 "Input vector must be nonempty"
   (flag, m) = is_square_with_sqrt(8*l+1)
   @req flag "Input vector of invalid length"
   n = div(m+1, 2)
   R = parent(L[1])
   M = zero_matrix(R, n, n)
   pos = 1
   for i in 1:(n-1), j in (i+1):n
      M[i,j] = L[pos]
      pos += 1
   end
   return M
end

###############################################################################
#
#   matrix_space constructor
#
###############################################################################

@doc raw"""
    matrix_space(R::NCRing, r::Int, c::Int)

Return parent object corresponding to the space of $r\times c$ matrices over
the ring $R$.
"""
function matrix_space(R::NCRing, r::Int, c::Int; cached::Bool = true)
  # TODO: the 'cached' argument is ignored and mainly here for backwards compatibility
  # (and perhaps future compatibility, in case we need it again)
  (r < 0 || c < 0) && error("Dimensions must be non-negative")
  T = elem_type(R)
  return MatSpace{T}(R, r, c)
end

###############################################################################
#
#   Matrix M = R[...] syntax
#
################################################################################

Base.typed_hvncat(R::NCRing, args...) = _matrix(R, hvncat(args...))
Base.typed_hvcat(R::NCRing, args...) = _matrix(R, hvcat(args...))
Base.typed_hcat(R::NCRing, args...) = _matrix(R, hcat(args...))
Base.typed_vcat(R::NCRing, args...) = _matrix(R, vcat(args...))
_matrix(R::NCRing, a::AbstractVector) = matrix(R, length(a), isempty(a) ? 0 : 1, a)
_matrix(R::NCRing, a::AbstractMatrix) = matrix(R, a)
