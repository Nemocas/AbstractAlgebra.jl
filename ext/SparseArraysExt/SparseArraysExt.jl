module SparseArraysExt

using AbstractAlgebra
using SparseArrays: SparseMatrixCSC, sparse

import AbstractAlgebra: matrix_repr

# Documented in src/generic/PermGroups.jl
function matrix_repr(::Type{M}, a::Perm{T}) where {M<:SparseMatrixCSC, T}
   n = length(a.d)
   return convert(M, sparse(1:n, a.d, ones(T, n)))
end

end # module
