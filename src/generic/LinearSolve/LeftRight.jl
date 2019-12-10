
function solve_left(a::AbstractAlgebra.MatElem{S}, b::AbstractAlgebra.MatElem{S}) where S <: RingElement
    @warn "Function has been depreciated. Use `solve_hnf(a,b,side=:left)` instead." maxlog=1
    return solve_hnf(a,b,side=:left)
end

function solve_left(a::AbstractAlgebra.MatElem{S}, b::AbstractAlgebra.MatElem{S}) where S <: FieldElement
    @warn "Function has been depreciated. Use `solve_lu(a,b,side=:left)` instead." maxlog=1
    return solve_lu(a,b,side=:left)
end


@doc Markdown.doc"""
    solve_left(a::AbstractAlgebra.MatElem{S}, b::AbstractAlgebra.MatElem{S}) where S <: RingElement
> Given an $r\times n$ matrix $a$ over a ring and an $m\times n$ matrix $b$
> over the same ring, return an $m\times r$ matrix $x$ such that $xa = b$. If
> no such matrix exists, an exception is raised.
"""
function old_solve_left(a::AbstractAlgebra.MatElem{S}, b::AbstractAlgebra.MatElem{S}) where S <: RingElement
   @assert ncols(a) == ncols(b)
   H, T = hnf_with_transform(a)
   b = deepcopy(b)
   z = zero(a, nrows(b), nrows(a))
   l = min(ncols(a), nrows(a))
   t = base_ring(a)()
   for i = 1:nrows(b)
      for j = 1:l
         k = 1
         while k <= ncols(H) && iszero(H[j, k])
            k += 1
         end
         if k > ncols(H)
            continue
         end
         q, r = divrem(b[i, k], H[j, k])
         r != 0 && error("Unable to solve linear system")
         z[i, j] = q
         q = -q
         for h = k:ncols(H)
            t = mul!(t, q, H[j, h])
            b[i, h] = addeq!(b[i, h], t)
         end
      end
   end
   b != 0 && error("Unable to solve linear system")
   return z*T
end

#=
function solve_left(A::AbstractAlgebra.MatElem{T}, B::AbstractAlgebra.MatElem{T}) where T <: FieldElement
  R = base_ring(A)
  ncols(A) != ncols(B) && error("Incompatible matrices")
  mu = zero_matrix(R, ncols(A), nrows(A) + nrows(B))
  for i = 1:ncols(A)
     for j = 1:nrows(A)
        mu[i, j] = A[j, i]
     end
     for j = 1:nrows(B)
        mu[i, nrows(A) + j] = B[j, i]
     end
  end
  rk, mu = rref(mu)
  p = find_pivot(mu)
  if any(i -> i > nrows(A), p)
    error("Unable to solve linear system")
  end
  sol = zero_matrix(R, nrows(B), nrows(A))
  for i = 1:length(p)
    for j = 1:nrows(B)
      sol[j, p[i]] = mu[i, nrows(A) + j]
    end
  end
  return sol
end

=#

# Find the pivot columns of an rref matrix
function find_pivot(A::AbstractAlgebra.MatElem{T}) where T <: RingElement
  p = Int[]
  j = 0
  for i = 1:nrows(A)
    j += 1
    if j > ncols(A)
      return p
    end
    while iszero(A[i, j])
      j += 1
      if j > ncols(A)
        return p
      end
    end
    push!(p, j)
  end
  return p
end
