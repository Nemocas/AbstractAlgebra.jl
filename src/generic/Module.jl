###############################################################################
#
#   Module.jl : Generic modules over Euclidean domains
#
###############################################################################

###############################################################################
#
#   Basic manipulation
#
###############################################################################

# The .v field in FPModuleElem has an abstract type. We introduce an internal
# type stable accessor function to make the compiler happy.
@inline function _matrix(v::AbstractAlgebra.FPModuleElem{T}) where T
  return (v.v)::dense_matrix_type(T)
end

@doc Markdown.doc"""
    rels(M::AbstractAlgebra.FPModule{T}) where T <: RingElement

Return a vector of all the relations between generators of the given
module, where each relation is given as row matrix. The relation matrix
whose rows are the returned relations will be in reduced form (hnf/rref).
"""
rels(M::AbstractAlgebra.FPModule{T}) where T <: RingElement = M.rels::Vector{dense_matrix_type(T)}

@doc Markdown.doc"""
    iscompatible(M::AbstractAlgebra.FPModule{T}, N::AbstractAlgebra.FPModule{T}) where T <: RingElement

Return `true, P` if the given modules are compatible, i.e. that they are
(transitively) submodules of the same module, P. Otherwise return `false, M`.
"""
function iscompatible(M::AbstractAlgebra.FPModule{T}, N::AbstractAlgebra.FPModule{T}) where T <: RingElement
   check_parent(M, N)
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

@doc Markdown.doc"""
    issubmodule(M::AbstractAlgebra.FPModule{T}, N::AbstractAlgebra.FPModule{T}) where T <: RingElement

Return `true` if $N$ was constructed as a submodule of $M$. The relation
is taken transitively (i.e. subsubmodules are submodules for the purposes
of this relation, etc). The module $M$ is also considered a submodule of
itself for this relation.
"""
function issubmodule(M::AbstractAlgebra.FPModule{T}, N::AbstractAlgebra.FPModule{T}) where T <: RingElement
   check_parent(M, N)
   if M === N
      return true
   end
   while isa(N, Submodule)
      N = supermodule(N)
      if M === N
         return true
      end
   end
   return false
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
# is the pivot column of the $i$-th row of the culled matrix.
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
   new_col = 1
   for i in 1:nrels
      while M[i, col] == 0
         push!(gen_cols, col)
         col += 1
         new_col += 1
      end
      if !isunit(M[i, col])
         push!(culled, i)
         push!(gen_cols, col)
         push!(pivots, new_col)
         new_col += 1
      end
      col += 1
   end
   while col <= ncols(M)
      push!(gen_cols, col)
      col += 1
   end
   return gen_cols, culled, pivots
end
