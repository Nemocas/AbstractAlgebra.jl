###############################################################################
#
#   Module.jl : Module functionality for modules over Euclidean domains
#
###############################################################################

export iscompatible

###############################################################################
#
#   Basic manipulation
#
###############################################################################

@doc Markdown.doc"""
    zero(M::AbstractAlgebra.FPModule{T}) where T <: RingElement
> Return the zero element of the module $M$.
"""
function zero(M::AbstractAlgebra.FPModule{T}) where T <: RingElement
   R = base_ring(M)
   return M([zero(R) for i in 1:ngens(M)])
end

@doc Markdown.doc"""
    relations(M::AbstractAlgebra.FPModule{T}) where T <: RingElement
> Return all relations between generators of the given module, where each
> relation is given as an element of `M`. The relation matrix whose rows
> are the returned relations will be in reduced form (hnf/rref).
"""
relations(M::AbstractAlgebra.FPModule{T}) where T <: RingElement = M.rels::Vector{dense_matrix_type(T)}

@doc Markdown.doc"""
    iscompatible(M::AbstractAlgebra.FPModule{T}, N::AbstractAlgebra.FPModule{T}) where T <: RingElement
> Return `true, P` if the given modules are compatible, i.e. that they are
> (transitively) submodules of the same module, P. Otherwise return `false, M`.
"""
function iscompatible(M::AbstractAlgebra.FPModule{T}, N::AbstractAlgebra.FPModule{T}) where T <: RingElement
   M1 = M
   M2 = N
   while isa(M1, Submodule)
      M2 = N
      while isa(M2, Submodule)
         if M1 === M2
            return true, M1
         end
         M2 = supermodule(M2)
      end
      M1 = supermodule(M1)
   end
   while isa(M2, Submodule)
      M2 = supermodule(M2)
   end
   if M1 === M2
      return true, M1
   end
   return false, M
end

function Base.intersect(M::AbstractAlgebra.FPModule{T}, N::AbstractAlgebra.FPModule{T}) where T <: RingElement
   # Compute the common supermodule P of M and N
   flag, P = iscompatible(M, N)
   !flag && error("Modules not compatible")
   # Compute the generators of M as elements of P
   G1 = gens(M)
   M1 = M
   while M1 !== P
      G1 = [M1.map(v) for v in G1]
      M1 = supermodule(M1)
   end
   # Compute the generators of N as elements of P
   G2 = gens(N)
   M2 = N
   while M2 !== P
      G2 = [M2.map(v) for v in G2]
      M2 = supermodule(M2)
   end
   # Make matrix containing all generators and relations as rows
   r1 = ngens(M)
   r2 = ngens(N)
   rels = relations(P)
   r3 = length(rels)
   c = ngens(P)
   mat = matrix(base_ring(M), r1 + r2 + r3, c, [0 for i in 1:(r1 + r2 + r3)*c])
   # We flip the rows of the matrix so the input to Submodule is in upper
   # triangular form
   rn = r1 + r2 + r3
   for i = 1:r1
      for j = 1:c
         mat[rn - i + 1, j] = G1[i].v[1, j]
      end
   end
   for i = 1:r2
      for j = 1:c
         mat[rn - i - r1 + 1, j] = G2[i].v[1, j]
      end
   end
   for i = 1:r3
      for j = 1:c
         mat[rn - i - r1 - r2 + 1, j] = rels[i][1, j]
      end
   end
   # Find the left kernel space of the matrix
   nc, K = left_kernel(mat)
   # Last r1 elements of a row correspond to a generators of intersection
   # We flip the rows of K so the input to Submodule is upper triangular
   # and the columns so that they correspond to the original order before
   # flipping above
   I = [M([K[nc - j + 1, rn - i + 1] for i in 1:r1]) for j in 1:nc]
   return Submodule(M, I)
end

function ==(M::AbstractAlgebra.FPModule{T}, N::AbstractAlgebra.FPModule{T}) where T <: RingElement
   # Compute the common supermodule P of M and N
   flag, P = iscompatible(M, N)
   !flag && error("Modules not compatible")
   # Compute the generators of M as elements of P
   G1 = gens(M)
   M1 = M
   while M1 !== P
      G1 = [M1.map(v) for v in G1]
      M1 = supermodule(M1)
   end
   # Compute the generators of N as elements of P
   G2 = gens(N)
   M2 = N
   while M2 !== P
      G2 = [M2.map(v) for v in G2]
      M2 = supermodule(M2)
   end
   # Put (rewritten) gens of M and N into matrices with relations of P
   rels = relations(P)
   c = ngens(P)
   r1 = ngens(M)
   r2 = ngens(N)
   mat1 = matrix(base_ring(M), r1 + length(rels), c,
                 [0 for i in 1:(r1 + length(rels))*c])
   for i = 1:r1
      for j = 1:c
         mat1[i, j] = G1[i].v[1, j]
      end
   end
   mat2 = matrix(base_ring(M), r2 + length(rels), c,
                 [0 for i in 1:(r2 + length(rels))*c])
   for i = 1:r2
      for j = 1:c
         mat2[i, j] = G2[i].v[1, j]
      end
   end
   for i = 1:length(rels)
      for j = 1:c
         mat1[i + r1, j] = rels[i][1, j]
         mat2[i + r2, j] = rels[i][1, j]
      end
   end
   # Put the matrices into reduced form
   mat1 = reduced_form(mat1)
   mat2 = reduced_form(mat2)
   # Check containment of rewritten gens of M in row space of mat2
   for v in G1
      flag, r = can_solve_left_reduced_triu(v.v, mat2)
      if !flag
         return false
      end
   end
   # Check containment of rewritten gens of N in row space of mat1
   for v in G2
      flag, r = can_solve_left_reduced_triu(v.v, mat1)
      if !flag
         return false
      end
   end
   return true
end

###############################################################################
#
#   Random generation
#
###############################################################################

function rand(M::AbstractAlgebra.FPModule{T}, vals...) where T <: RingElement
   R = base_ring(M)
   v = [rand(R, vals...) for i in 1:ngens(M)]
   return M(v)
end

###############################################################################
#
#   Helper functions
#
###############################################################################

# Assumes M is in reduced form (hnf/rref). Removes zero rows. Returns a tuple
# gen_cols, culled, pivots where all rows and columns corresponding to unit
# pivots have been removed, gen_cols is a list of columns without unit pivots,
# culled is an array of row (indices) that have not been removed and pivots[i]
# is the pivot column of the $i$-th row of the culled matrix
function cull_matrix(M::AbstractAlgebra.MatElem{T}) where T <: RingElement
   # count the nonzero rows
   nrels = nrows(M)
   while nrels > 0 && iszero_row(M, nrels)
      nrels -= 1
   end
   # find relations with non-unit pivot
   gen_cols = Vector{Int}(undef, 0)
   culled = Vector{Int}(undef, 0)
   pivots = Vector{Int}(undef, 0)
   col = 1
   row = 1
   new_col = 1
   for i in 1:nrels
      while M[i, col] == 0
         push!(gen_cols, col)
         col += 1
         new_col += 1
      end
      if !isunit(M[i, col])
         push!(culled, row)
         push!(gen_cols, col)
         push!(pivots, new_col)
         new_col += 1
      end
      col += 1
      row += 1
   end
   while col <= ncols(M)
      push!(gen_cols, col)
      col += 1
   end
   return gen_cols, culled, pivots
end

# As per the previous function but operates on a vector of rows.
# The duplication is unavoidable until we get views. The parameter
# ngens is just the number of columns of the matrix (which cannot
# otherwise be determined if the vector M is empty).
function cull_matrix_rows(M::Vector{<:AbstractAlgebra.MatElem{T}}, ngens::Int) where T <: RingElement
   # count the nonzero rows
   nrels = length(M)
   while nrels > 0 && iszero(M[nrels])
      nrels -= 1
   end
   # find relations with non-unit pivot
   gen_cols = Vector{Int}(undef, 0)
   culled = Vector{Int}(undef, 0)
   pivots = Vector{Int}(undef, 0)
   col = 1
   row = 1
   new_col = 1
   for i in 1:nrels
      while M[i][col] == 0
         push!(gen_cols, col)
         col += 1
         new_col += 1
      end
      if !isunit(M[i][col])
         push!(culled, row)
         push!(gen_cols, col)
         push!(pivots, new_col)
         new_col += 1
      end
      col += 1
      row += 1
   end
   while col <= ngens
      push!(gen_cols, col)
      col += 1
   end
   return gen_cols, culled, pivots
end
