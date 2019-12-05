
function solve_ff(M::MatrixElem{T}, b::MatrixElem{T}) where {T <: FieldElement}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   !issquare(M) && error("Non-square matrix in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   m = nrows(M)
   x, d = solve_fflu(M, b)
   for i in 1:nrows(x)
      for j in 1:ncols(x)
         x[i, j] = divexact(x[i, j], d)
      end
   end
   return x
end

function solve_ff(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: RingElement}
   m = nrows(M)
   n = ncols(M)
   if m == 0 || n == 0
      return b, base_ring(M)()
   end
   return solve_fflu(M, b)
end


function solve_with_det(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: RingElement}
   # We cannot use solve_fflu directly, since it forgot about the (parity of
   # the) permutation.
   nrows(M) != ncols(M) && error("Non-square matrix")
   R = base_ring(M)
   FFLU = deepcopy(M)
   p = PermGroup(nrows(M))()
   r, d = fflu!(p, FFLU)
   if r < nrows(M)
      error("Singular matrix in solve_with_det")
   end
   x = solve_fflu_precomp(p, FFLU, b)
   # Now M*x = d*b, but d is only sign(P) * det(M)
   if parity(p) != 0
      minus_one = R(-1)
      for k in 1:ncols(x)
         for i in 1:nrows(x)
            # We are allowed to modify x in-place.
            x[i, k] = mul!(x[i, k], x[i, k], minus_one)
         end
      end
      d = mul!(d, d, minus_one)
   end
   return x, d
end

function solve_with_det(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: PolyElem}
   x, d = solve_interpolation(M, b)
   return x, d
end


function solve_interpolation(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: PolyElem}
   m = nrows(M)
   h = ncols(b)
   if m == 0
      return b, base_ring(M)()
   end
   R = base_ring(M)
   maxlen = 0
   for i = 1:m
      for j = 1:m
         maxlen = max(maxlen, length(M[i, j]))
      end
   end
   maxlen == 0 && error("Matrix is singular in solve")
   maxlenb = 0
   for i = 1:m
      for j = 1:h
         maxlenb = max(maxlenb, length(b[i, j]))
      end
   end
   # bound from xd = (M*)b where d is the det
   bound = (maxlen - 1)*(m - 1) + max(maxlenb, maxlen)
   tmat = matrix(base_ring(R), 0, 0, elem_type(base_ring(R))[])
   V = Array{typeof(tmat)}(undef, bound)
   d = Array{elem_type(base_ring(R))}(undef, bound)
   y = Array{elem_type(base_ring(R))}(undef, bound)
   bj = Array{elem_type(base_ring(R))}(undef, bound)
   X = similar(tmat, m, m)
   Y = similar(tmat, m, h)
   x = similar(b)
   b2 = div(bound, 2)
   pt1 = base_ring(R)(1 - b2)
   l = 1
   i = 1
   pt = 1
   while l <= bound
      y[l] = base_ring(R)(pt - b2)
      (y[l] == pt1 && pt != 1) && error("Not enough interpolation points in ring")
      for j = 1:m
         for k = 1:m
            X[j, k] = evaluate(M[j, k], y[l])
         end
         for k = 1:h
            Y[j, k] = evaluate(b[j, k], y[l])
         end
      end
      try
         V[l], d[l] = solve_with_det(X, Y)
         l = l + 1
      catch e
         if !(e isa ErrorException)
            rethrow(e)
         end
         i = i + 1
      end

      # We tested bound evaluation points and either an impossible
      # inverse was encountered, or the matrix was singular for all
      # the values.

      if i > bound && l == 1
         error("Impossible inverse or too many singular matrices in solve_interpolation")
      end

      pt = pt + 1
   end
   for k = 1:h
      for i = 1:m
         for j = 1:bound
            bj[j] = V[j][i, k]
         end
         x[i, k] = interpolate(R, y, bj)
      end
   end
   return x, interpolate(R, y, d)
end

@doc Markdown.doc"""
    solve(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: FieldElement}
> Given a non-singular $n\times n$ matrix over a field and an $n\times m$
> matrix over the same field, return $x$ an
> $n\times m$ matrix $x$ such that $Ax = b$.
> If $A$ is singular an exception is raised.
"""
function solve(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: FieldElement}
    return solve_ringelem(M, b)
end

@doc Markdown.doc"""
    solve_rational(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where T <: RingElement
> Given a non-singular $n\times n$ matrix over a ring and an $n\times m$
> matrix over the same ring, return a tuple $x, d$ consisting of an
> $n\times m$ matrix $x$ and a denominator $d$ such that $Ax = db$. The
> denominator will be the determinant of $A$ up to sign. If $A$ is singular an
> exception is raised.
"""
function solve_rational(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where T <: RingElement
   return solve_ringelem(M, b)
end

function solve_ringelem(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: RingElement}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   nrows(M) != ncols(M) && error("Non-square matrix in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   return solve_ff(M, b)
end

function solve_rational(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: PolyElem}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   nrows(M) != ncols(M) && error("Non-square matrix in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   try
      return solve_interpolation(M, b)
   catch e
      if !isa(e, ErrorException)
         rethrow(e)
      end
      return solve_ff(M, b)
   end
end
