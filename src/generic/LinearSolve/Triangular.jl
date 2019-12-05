

# These functions are honestly just triangular solving in disguise.

function solve_fflu_precomp(p::Generic.Perm, FFLU::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   x = p * b
   n = nrows(x)
   m = ncols(x)
   R = base_ring(FFLU)

   t = base_ring(b)()
   s = base_ring(b)()
   minus_one = R(-1)

    # For each column of x
    for k in 1:m
        # Backsolve the lower triangular part. However, things have been rearranged
        # so that a single column is accessed in the inner (j) loop, rather than a row.
      for i in 1:(n - 1)
         t = mul!(t, x[i, k], minus_one)
         for j in (i + 1):n
            if i == 1
              x[j, k] = mul_red!(R(), x[j, k], FFLU[i, i], false)
            else
              x[j, k] = mul_red!(x[j, k], x[j, k], FFLU[i, i], false)
            end
            s = mul_red!(s, FFLU[j, i], t, false)
            x[j, k] = addeq!(x[j, k], s)
            x[j, k] = reduce!(x[j, k])
            if i > 1
                x[j, k] = divexact(x[j, k], FFLU[i - 1, i - 1])
            end
         end
      end

        # Backsolve the upper triangular part. However, things have been rearranged
        # so that a single column of `x` is accessed in the inner (j) loop, rather than a row.
      for i in (n - 1):-1:1
         if i > 1
            x[i, k] = mul!(x[i, k], x[i, k], FFLU[n, n])
         else
            x[i, k] = x[i, k] * FFLU[n, n]
         end
         for j in (i + 1):n
            t = mul!(t, x[j, k], FFLU[i, j])
            t = mul!(t, t, minus_one)
            x[i, k] = addeq!(x[i, k], t)
         end
         x[i, k] = divexact(x[i, k], FFLU[i, i])
      end
   end
   return x
end

function solve_lu_precomp(p::Generic.Perm, LU::MatElem{T}, b::MatrixElem{T}) where {T <: FieldElement}
    
   x = p * b
   n = nrows(x)
   m = ncols(x)
   R = base_ring(LU)

   t = base_ring(b)()
   s = base_ring(b)()

    # For each column of x.
   for k in 1:m
       x[1, k] = deepcopy(x[1, k])

       # The back-substitution, along rows.
      for i in 2:n
         for j in 1:(i - 1)
            # x[i, k] = x[i, k] - LU[i, j] * x[j, k]
            t = mul_red!(t, -LU[i, j], x[j, k], false)
            if j == 1
               x[i, k] = x[i, k] + t # LU[i, j] * x[j, k]
            else
               x[i, k] = addeq!(x[i, k], t)
            end
         end
         x[i, k] = reduce!(x[i, k])
      end

      # Now every entry of x is a proper copy, so we can change the entries
      # as much as we want.

      x[n, k] = divexact(x[n, k], LU[n, n])

      for i in (n - 1):-1:1
         for j in (i + 1):n
            # x[i, k] = x[i, k] - x[j, k] * LU[i, j]
            t = mul_red!(t, x[j, k], -LU[i, j], false)
            x[i, k] = addeq!(x[i, k], t)
         end
         x[i, k] = reduce!(x[i, k])
         x[i, k] = divexact(x[i, k], LU[i, i])
      end
   end
   return x
end



###############################################################################
#
#   Upper triangular solving
#
###############################################################################

@doc Markdown.doc"""
    solve_triu(U::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}, unit::Bool = false) where {T <: FieldElement}
> Given a non-singular $n\times n$ matrix over a field which is upper
> triangular, and an $n\times m$ matrix over the same field, return an
> $n\times m$ matrix $x$ such that $Ax = b$. If $A$ is singular an exception
> is raised. If unit is true then $U$ is assumed to have ones on its
> diagonal, and the diagonal will not be read.
"""
function solve_triu(U::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}, unit::Bool = false) where {T <: FieldElement}
   n = nrows(U)
   m = ncols(b)
   R = base_ring(U)
   X = zero(b)
   Tinv = Array{elem_type(R)}(undef, n)
   tmp = Array{elem_type(R)}(undef, n)
   if unit == false
      for i = 1:n
         Tinv[i] = inv(U[i, i])
      end
   end
   t = R()
   for i = 1:m
      for j = 1:n
         tmp[j] = X[j, i]
      end
      for j = n:-1:1
         s = R()
         for k = j + 1:n
            s = addmul_delayed_reduction!(s, U[j, k], tmp[k], t)
         end
         s = reduce!(s)
         s = b[j, i] - s
         if unit == false
            s = mul!(s, s, Tinv[j])
         end
         tmp[j] = s
      end
      for j = 1:n
         X[j, i] = tmp[j]
      end
   end
   return X
end

###############################################################################
#
#   Can solve
#
###############################################################################

@doc Markdown.doc"""
    can_solve_left_reduced_triu(r::AbstractAlgebra.MatElem{T},
                          M::AbstractAlgebra.MatElem{T}) where T <: RingElement
> Return a tuple `flag, x` where `flag` is set to true if $xM = r$ has a
> solution, where $M$ is an $m\times n$ matrix in (upper triangular) Hermite
> normal form or reduced row echelon form and $r$ and $x$ are row vectors with
> $m$ columns. If there is no solution, flag is set to `false` and $x$ is set
> to the zero row.
"""
function can_solve_left_reduced_triu(r::AbstractAlgebra.MatElem{T},
                          M::AbstractAlgebra.MatElem{T}) where T <: RingElement
   ncols(r) != ncols(M) && error("Incompatible matrices")
   r = deepcopy(r) # do not destroy input
   m = ncols(r)
   n = nrows(M)
   if n == 0
      return true, r
   end
   R = base_ring(r)
   x = zero_matrix(R, 1, n)
   j = 1 # row in M
   k = 1 # column in M
   t = R()
   for i = 1:m # column in r
      if iszero(r[1, i])
         continue
      end
      while k <= i && j <= n
         if iszero(M[j, k])
            k += 1
         elseif k < i
            j += 1
         else
            break
         end
      end
      if k != i
         return false, x
      end
      x[1, j], r[1, i] = divrem(r[1, i], M[j, k])
      if !iszero(r[1, i])
         return false, x
      end
      q = -x[1, j]
      for l = i + 1:m
         t = mul!(t, q, M[j, l])
         r[1, l] = addeq!(r[1, l], t)
      end
   end
   return true, x
end
