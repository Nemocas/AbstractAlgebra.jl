module SparseArraysExt

using AbstractAlgebra
using AbstractAlgebra.Generic: YoungTableau, SkewDiagram
using SparseArrays: sparse, spzeros

import AbstractAlgebra: matrix_repr

# Docstrings for these methods are in src/generic/PermGroups.jl and
# src/generic/YoungTabs.jl.

matrix_repr(a::Perm{T}) where {T<:Integer} = sparse(collect(T, 1:length(a.d)), a.d, ones(T, length(a.d)))

function matrix_repr(Y::YoungTableau{T}) where {T<:Integer}
   tab = spzeros(T, length(Y.part), Y.part[1])
   k = 1
   for (idx, p) in enumerate(Y.part)
      tab[idx, 1:p] = Y.fill[k:k+p-1]
      k += p
   end
   return tab
end

function matrix_repr(xi::SkewDiagram)
   skdiag = spzeros(eltype(xi), size(xi)...)
   for i in 1:length(xi.mu)
      skdiag[i, xi.mu[i]+1:xi.lam[i]] .= 1
   end
   for i in length(xi.mu)+1:length(xi.lam)
      skdiag[i, 1:xi.lam[i]] .= 1
   end
   return skdiag
end

end # module
