export scalar_matrix, diagonal_matrix, iszero_row

################################################################################
#
#  Matrix from Hecke
#
################################################################################
"""
This file contains functions transported from Hecke. It will act as a kind of
quarantine zone until the AbstractAlgebra package maintainers are happy with
the state of the code.
"""

################################################################################
#
#  Dense matrix types
#
################################################################################

dense_matrix_type(::Type{T}) where {T} = Generic.MatSpaceElem{T}

coefficient_type(::Type{Generic.Mat{T}}) where {T} = T

################################################################################
#
#  Matrix constructors
#
################################################################################

@doc Markdown.doc"""
    zero_matrix(::Type{MatElem}, R::Ring, n::Int)
    zero_matrix(::Type{MatElem}, R::Ring, n::Int, m::Int)
Return the zero `n x n` (resp. `n x m`) matrix over the ring `R`.
"""
function zero_matrix(::Type{MatElem}, R::Ring, n::Int)
  return zero_matrix(R, n)
end

function zero_matrix(::Type{MatElem}, R::Ring, n::Int, m::Int)
  return zero_matrix(R, n, m)
end

@doc Markdown.doc"""
    matrix(A::Array{T, 2}) where T <: RingElem
Construct an AbstractAlgebra matrix whose entries are `A`. Note that `A`
must be non-empty.
"""
function matrix(A::Array{T, 2}) where T <: RingElem
  r, c = size(A)
  isempty(A) && throw(DomainError(A, "Array must be non-empty."))
  m = matrix(parent(A[1, 1]), A)
  return m
end

@doc Markdown.doc"""
    matrix(A::Array{T, 1}) where T <: RingElem
Construct an `n x 1` AbstractAlgebra matrix `M` such that `M[j,1] = A[j]`.
Note that `A` must be non-empty.
"""
function matrix(A::Array{T, 1}) where T <: RingElem
  return matrix(reshape(A,length(A),1))
end

@doc Markdown.doc"""
    scalar_matrix(R::Ring, n::Int, a::RingElement)
Construct the AbstractAlgebra matrix `a*I`, with `I` the `n x n` identity matrix.
"""
function scalar_matrix(R::Ring, n::Int, a::RingElement)
  b = R(a)
  z = zero_matrix(R, n, n)
  for i in 1:n
    z[i, i] = b
  end
  return z
end

################################################################################
#
#  Zero checking.
#
################################################################################

@doc Markdown.doc"""
    iszero_row(M::MatElem{T}, i::Int) where T
Check if the `i`-th row of the matrix `M` is zero. Note this function is much
faster than either `iszero(M[1,:])` or `iszero(A.entries[i,:])`.
"""

function iszero_row(M::Union{MatElem{T}, Array{T,2}}, i::Int) where T
  for j in 1:ncols(M)
    if !iszero(M[i,j])
      return false
    end
  end
  return true
end


################################################################################
#
#  Diagonal (block) matrix
#
################################################################################

@doc Markdown.doc"""
    diagonal_matrix(x::Vector{T}) where T <: NCRingElem -> MatElem{T}
    diagonal_matrix(x::T...) where T <: NCRingElem -> MatElem{T}
Returns a diagonal matrix whose diagonal entries are the element of `x`.
"""
function diagonal_matrix(x::Vector{T}) where T <: NCRingElem
  isempty(x) && throw(DomainError(x, "Array must be non-empty."))
  M = zero_matrix(parent(x[1]), length(x), length(x))
  for i = 1:length(x)
    M[i, i] = x[i]
  end
  return M
end

function diagonal_matrix(x::T...) where T <: NCRingElem
  return diagonal_matrix(collect(x))
end

@doc Markdown.doc"""
    diagonal_matrix(x::Vector{T}) where T <: MatElem -> MatElem
    diagonal_matrix(x::T...) where T <: MatElem -> MatElem

Returns a block diagonal matrix whose diagonal blocks are the matrices in $x$.
"""
function diagonal_matrix(x::Vector{T}) where T <: MatElem
  return diagonal_matrix(x...)
end

function diagonal_matrix(x::T...) where T <: MatElem
    (isempty(x) || all([isempty(M) for M in x])) && (
        throw(DomainError(x, "Not all matrices can be empty")))
    
    return cat(x..., dims = (1, 2))
end


@doc Markdown.doc"""
    isdiagonal(A::Mat)

Tests if `A` is diagonal; i.e, if `A[i,j] == 0` for all `i != j`, which is well-defined
for rectangular matrices. We consider the empty matrix to be diagonal as well.
"""
function isdiagonal(A::MatElem)
    for i = 1:ncols(A)
        #Pre-diagonal.
        for j = 1:i-1
            !iszero(A[j, i]) && return false
        end
        #Post-diagonal.
        for j = i+1:nrows(A)
            !iszero(A[j, i]) && return false
        end
    end
    return true
end


################################################################################
#
#  Unsafe arithmetic functions for generic matrices
#
################################################################################

#=
function mul!(c::MatElem, a::MatElem, b::MatElem)

    if any(ncols(a) != nrows(b), nrows(c) != nrows(a), ncols(c) != ncols(b))
        @info "Matrix dimensions:" size(a) size(b) size(c)
        throw(DomainError((a,b,c), "Incompatible matrix dimensions"))
    end

    if c === a || c === b
        # If the output is either of the inputs, memory must be allocated anyway.
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

function add!(c::MatElem, a::MatElem, b::MatElem)
    !(parent(a) == parent(b) == parent(c)) && throw(DomainError(
        (parent(a), parent(b), parent(c)), "Parents don't match."))

    for i = 1:nrows(c)
        for j = 1:ncols(c)
            c[i, j] = add!(c[i, j], a[i, j], b[i, j])
        end
    end
    return c
end

function sub!(c::MatElem, a::MatElem, b::MatElem)
    !(parent(a) == parent(b) == parent(c)) && throw(DomainError(
        (parent(a), parent(b), parent(c)), "Parents don't match."))

    for i = 1:nrows(c)
        for j = 1:ncols(c)
            c[i, j] = sub!(c[i, j], a[i, j], b[i, j])
        end
    end
    return c
end
=#
