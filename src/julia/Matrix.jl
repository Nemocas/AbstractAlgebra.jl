###############################################################################
#
#   Matrix.jl : Additional AbstractAlgebra functionality for Julia Matrices
#
###############################################################################

nrows(A::Matrix{T}) where {T} = size(A)[1]
ncols(A::Matrix{T}) where {T} = size(A)[2]

function is_zero_row(M::Matrix, i::Int)
    for j = 1:ncols(M)
        if !iszero(M[i, j])
            return false
        end
    end
    return true
end

# TODO: add `is_zero_column(M::Matrix{T}, i::Int) where {T<:Integer}` and specialized functions in Nemo

zero_matrix(::Type{Int}, r, c) = zeros(Int, r, c)
