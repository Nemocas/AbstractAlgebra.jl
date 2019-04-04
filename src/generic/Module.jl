###############################################################################
#
#   Module.jl : Module functionality for modules over Euclidean domains
#
###############################################################################

export iscompatible

@doc Markdown.doc"""
    iscompatible(M::AbstractAlgebra.Module{T}, N::AbstractAlgebra.Module{T}) where T <: RingElement
> Return `true, P` if the given modules are compatible, i.e. that they are
> (transitively) submodules of the same module, P. Otherwise return `false, M`.
"""
function iscompatible(M::AbstractAlgebra.Module{T}, N::AbstractAlgebra.Module{T}) where T <: RingElement
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

function Base.intersect(M::AbstractAlgebra.Module{T}, N::AbstractAlgebra.Module{T}) where T <: RingElement
   # Compute the common supermodule P of M and N
   flag, P = iscompatible(M, N)
   !flag && error("Modules not compatible")
   # Compute the generators of M as elements of P
   G1 = gens(M)
   M1 = M
   while M1 != P
      G1 = [M1.map(v) for v in G1]
      M1 = supermodule(M1)
   end
   # Compute the generators of N as elements of P
   G2 = gens(N)
   M2 = N
   while M2 != P
      G2 = [M2.map(v) for v in G2]
      M2 = supermodule(M2)
   end
   # Make matrix containing all generators as columns
   r1 = ngens(M)
   r2 = ngens(N)
   c = ngens(P)
   mat = matrix(base_ring(M), c, r1 + r2, [0 for i in 1:c*(r1 + r2)])
   for i = 1:r1
      for j = 1:c
         mat[j, i] = G1[i].v[1, j]
      end
   end
   for i = 1:r2
      for j = 1:c
         mat[j, i + r1] = G2[i].v[1, j]
      end
   end
   # Find the left kernel space of the matrix
   # z = matrix(base_ring(M), r1 + r2, 1, [0 for i in 1:r1 + r2])
   # K = solve(mat, z) # generators of kernel space are columns
   nc, K = nullspace(mat)
   # First r1 elements of a column correspond to a generators of intersection
   I = [M([K[i, j] for i in 1:r1]) for j in 1:nc]
   return Submodule(M, I)
end

