###############################################################################
#
#   Matrix.jl : Additional AbstractAlgebra functionality for Julia Matrices
#
###############################################################################

number_of_rows(A::Matrix{T}) where {T} = size(A)[1]
number_of_columns(A::Matrix{T}) where {T} = size(A)[2]

zero_matrix(::Type{Int}, r, c) = zeros(Int, r, c)

###############################################################################
#
#   Conversion from MatrixElem
#
###############################################################################

"""
    Matrix(A::MatrixElem{T}) where {T<:NCRingElement}
    Matrix{U}(A::MatrixElem{T}) where {U<:NCRingElement, T<:NCRingElement}

Convert `A` to a Julia `Matrix{U}` of the same dimensions with the same elements.
If `U` is omitted then `eltype(M)` is used in its place.

# Examples
```jldoctest
julia> A = ZZ[1 2 3; 4 5 6]
[1   2   3]
[4   5   6]

julia> Matrix(A)
2×3 Matrix{BigInt}:
 1  2  3
 4  5  6

julia> Matrix{Int}(A)
2×3 Matrix{Int64}:
 1  2  3
 4  5  6
```
"""
Matrix(M::MatrixElem{T}) where {T<:NCRingElement} = Matrix{eltype(M)}(M)
Matrix{U}(M::MatrixElem{T}) where {U<:NCRingElement, T<:NCRingElement} = U[M[i, j] for i = 1:nrows(M), j = 1:ncols(M)]


"""
    Array(A::MatrixElem{T}) where T <: NCRingElement

Convert `A` to a Julia `Matrix` of the same dimensions with the same elements.

# Examples
```jldoctest
julia> R, x = ZZ[:x]; A = R[x^0 x^1; x^2 x^3]
[  1     x]
[x^2   x^3]

julia> Array(A)
2×2 Matrix{AbstractAlgebra.Generic.Poly{BigInt}}:
 1    x
 x^2  x^3
```
"""
Array(M::MatrixElem{T}) where {T<:NCRingElement} = Matrix(M)


###############################################################################
#
#   Array creation functions
#
###############################################################################

Array(R::NCRing, r::Int...) = Array{elem_type(R)}(undef, r)

function zeros(R::NCRing, r::Int...)
   T = elem_type(R)
   A = Array{T}(undef, r)
   for i in eachindex(A)
      A[i] = R()
   end
   return A
end

###############################################################################
#
#   Row & column permutations
#
###############################################################################

function swap_rows(a::Matrix{T}, i::Int, j::Int) where T <: NCRingElement
   (1 <= i <= nrows(a) && 1 <= j <= nrows(a)) || throw(BoundsError())
   b = deepcopy(a)
   swap_rows!(b, i, j)
   return b
end

function swap_rows!(a::Matrix{T}, i::Int, j::Int) where T <: NCRingElement
   if i != j
      for k = 1:ncols(a)
         a[i, k], a[j, k] = a[j, k], a[i, k]
      end
   end
   return a
end

function swap_cols(a::Matrix{T}, i::Int, j::Int) where T <: NCRingElement
   (1 <= i <= ncols(a) && 1 <= j <= ncols(a)) || throw(BoundsError())
   b = deepcopy(a)
   swap_cols!(b, i, j)
   return b
end

function swap_cols!(a::Matrix{T}, i::Int, j::Int) where T <: NCRingElement
   if i != j
      for k = 1:nrows(a)
         a[k, i], a[k, j] = a[k, j], a[k, i]
      end
   end
   return a
end
