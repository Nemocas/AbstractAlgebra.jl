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

# TODO LIST:
#= Pure abstract algebra.

T = tested.

DONE:
XX Unsafe functions for generic matrices [basically, cannot implement.]
T- Zero matrix constructors
T- iszero_row
T- diagonal_matrix, isdiagonal, diagonal
T- concatination [Already there]
-- array interface [Some already there]
-- reduce_mod
-- find_pivot
-- can_solve and friends... (Big...)
-- minpoly/charpoly
-- basic eigenvector

NOT DONE:
-- Kernel function.
-- Kernel with base ring
-- is upper/lower triangular
-- triangular solving

-- where the hell is the Array/Matrix interface promised in the Documentation?

=#

#= Also involves NEMO

NOT DONE:
-- Dense matrix types
-- Saturation
-- Zero matrix constructors
-- iszero_row [FLINT call]
-- hnf
-- divexact
-- is_lll_reduced
-- maximum / minimum
-- lifting to over-rings.
-- copy
-- powering
-- round_scale
-- shift
-- mod, mod_sym
-- Smith normal form

=#


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
    diagonal_matrix(R::Ring, x::Vector)
    diagonal_matrix(x::Vector{T}) where T <: NCRingElem -> MatElem{T}
    diagonal_matrix(x::T...) where T <: NCRingElem -> MatElem{T}
Returns a diagonal matrix whose diagonal entries are the element of `x`. In the
first method, the elements of `x` are assigned via `setindex!`, so conversion is
automatic.
"""
function diagonal_matrix(R::Ring, x::Vector)
  M = zero_matrix(R, length(x), length(x))
  for i = 1:length(x)
    M[i, i] = x[i]
  end
  return M
end

function diagonal_matrix(x::Vector{T}) where T <: NCRingElem
    isempty(x) && throw(DomainError(x, "Array must be non-empty."))
    return diagonal_matrix(parent(x[1]), x)
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

@doc Markdown.doc"""
    diagonal(A::Mat{T}) -> Vector{T}

Returns the diagonal of `A` is as array.
"""
diagonal(A::Generic.Mat{T}) where {T} = T[A[i, i] for i in 1:nrows(A)]


################################################################################
#
#  hvcat
#
################################################################################

# Via grep, this method is completely unused in Hecke. We add it here to annoy
# people until they remove it.

function Base.hvcat(rows::Tuple{Vararg{Int}}, A::MatElem...)

    @warn "This function is not actually used for anything." 
    B = hcat([A[i] for i=1:rows[1]]...)
  o = rows[1]
  for j=2:length(rows)
    C = hcat([A[i+o] for i=1:rows[j]]...)
    o += rows[j]
    B = vcat(B, C)
  end
  return B
end

################################################################################
#
#  Array interface (extensions)
#
################################################################################

function Base.keys(A::MatElem)
    return keys(A.entries)
end

function Base.getindex(A::MatElem, I::CartesianIndex{2})
    return A[I[1], I[1]]
end

function Base.getindex(A::MatElem, n::Int)
    1 <= n <= length(A) || throw(BoundsError(A,n))
    return A[1 + ((n-1) % nrows(A)), 1 + div((n-1), nrows(A))]
end

function setindex!(A::MatElem{T}, n::Int, s::T) where T <: RingElem
    1 <= n <= length(A) || throw(BoundsError(A,n))
    A[1 + ((n-1) % nrows(A)), 1 + div((n-1), nrows(A))] = s
    return s
end

function Base.stride(A::MatElem, n::Int)
    n <= 1 && return 1
    n == 2 && return nrows(A)
    return length(A)
end

##


function iterate(A::MatElem, state::Int = 0)
    
    # Annoy Hecke devs until they change code.
    if state == 0
        @warn "`iterate` output shape has changed for matrices. Please adjust accordingly."
    end
    
    if state < length(A)
        state += 1
        return A[state], state
    end
    return nothing
end

Base.IteratorSize(M::MatElem) = Base.HasShape{2}()
Base.IteratorEltype(M::MatElem) = Base.HasEltype()
Base.eltype(M::MatElem) = elem_type(base_ring(M))


################################################################################
#
#  Minpoly and Charpoly
#
################################################################################

function minpoly(M::MatElem)
  k = base_ring(M)
  kx, x = PolynomialRing(k, cached = false)
  return minpoly(kx, M)
end

function charpoly(M::MatElem)
  k = base_ring(M)
  kx, x = PolynomialRing(k, cached = false)
  return charpoly(kx, M)
end
