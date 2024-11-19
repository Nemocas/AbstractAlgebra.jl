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
