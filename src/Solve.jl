module Solve

using AbstractAlgebra

@doc raw"""
    solve(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement
    solve(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T <: FieldElement

Return $x$ of same type as $b$ solving the linear system $Ax = b$, if `side == :right`
(default), or $xA = b$, if `side == :left`.

If no solution exists, an error is raised.

See also [`can_solve_with_solution`](@ref).
"""
function solve(A::MatElem{T}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :right) where T <: FieldElement
  fl, x = can_solve_with_solution(A, b, side = side)
  fl || throw(ArgumentError("Unable to solve linear system"))
  return x
end

@doc raw"""
    can_solve(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement
    can_solve(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T <: FieldElement

Return `true` if the linear system $Ax = b$ or $xA = b$ with `side == :right`
(default) or `side == :left`, respectively, has a solution and `false` otherwise.

See also [`can_solve_with_solution`](@ref).
"""
function can_solve(A::MatElem{T}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :right) where T <: FieldElement
  return _can_solve_internal(A, b, :only_check; side = side)[1]
end

@doc raw"""
    can_solve_with_solution(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement
    can_solve_with_solution(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T <: FieldElement

Return `true` and $x$ of same type as $b$ solving the linear system $Ax = b$, if
such a solution exists. Return `false` and an empty vector or matrix, if the
system has no solution.

If `side == :left`, the system $xA = b$ is solved.

See also [`solve`](@ref).
"""
function can_solve_with_solution(A::MatElem{T}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :right) where T <: FieldElement
  return _can_solve_internal(A, b, :with_solution; side = side)[1:2]
end

@doc raw"""
    can_solve_with_solution_and_kernel(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T <: FieldElement
    can_solve_with_solution_and_kernel(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T <: FieldElement

Return `true`, $x$ of same type as $b$ solving the linear system $Ax = b$,
together with a matrix $K$ giving the kernel of $A$ (i.e. $AK = 0$), if such
a solution exists. Return `false`, an empty vector or matrix and an empty matrix,
if the system has no solution.

If `side == :left`, the system $xA = b$ is solved.

See also [`solve`](@ref) and [`kernel`](@ref).
"""
function can_solve_with_solution_and_kernel(A::MatElem{T}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :right) where T <: FieldElement
  return _can_solve_internal(A, b, :with_kernel; side = side)
end

# Tries to solve Ax = b (side == :right) or xA = b (side == :left) possibly with kernel.
# Always returns a tuple (Bool, MatElem, MatElem).
# task may be:
# * :only_check -> It is only tested whether there is a solution, the two MatElem's are
#   "dummies"
# * :with_solution -> A solution is computed, the last MatElem is a "dummy"
# * :with_kernel -> A solution and the kernel is computed
function _can_solve_internal(A::MatElem{T}, b::MatElem{T}, task::Symbol; side::Symbol = :right) where T <: FieldElement
  if task !== :only_check && task !== :with_solution && task !== :with_kernel
    error("task $(task) not recognized")
  end

  if side !== :right && side !== :left
    throw(ArgumentError("Unsupported argument :$side for side: Must be :left or :right."))
  end

  isright = side === :right

  R = base_ring(A)

  # Build [ A | b ] and compute rref
  if isright
    nrows(A) != nrows(b) && error("Incompatible matrices")
    mu = hcat(A, b)
    ncolsA = ncols(A)
    ncolsb = ncols(b)
  else
    ncols(A) != ncols(b) && error("Incompatible matrices")
    mu = _hcat_transposed(A, b)
    # For side == :left we pretend that A and b are transposed!
    ncolsA = nrows(A)
    ncolsb = nrows(b)
  end

  rk = rref!(mu)
  p = pivot_and_non_pivot_cols(mu, rk)
  if any(i -> i > ncolsA, p[1:rk])
    return false, zero_matrix(R, 0, 0), zero_matrix(R, 0, 0)
  end
  if task === :only_check
    return true, zero_matrix(R, 0, 0), zero_matrix(R, 0, 0)
  end

  # Compute a solution
  sol = isright ? zero_matrix(R, ncols(A), ncols(b)) : zero_matrix(R, nrows(b), nrows(A))
  for i = 1:rk
    for j = 1:ncolsb
      if isright
        sol[p[i], j] = mu[i, ncols(A) + j]
      else
        sol[j, p[i]] = mu[i, nrows(A) + j]
      end
    end
  end
  if task === :with_solution
    return true, sol, zero_matrix(R, 0, 0)
  end

  # Build the kernel
  X = isright ? zero_matrix(R, ncols(A), nullity) : zero_matrix(R, nullity, nrows(A))
  for i = 1:nullity
    for j = 1:rk
      if isright
        X[p[j], i] = -mu[j, p[rk + i]]
      else
        X[i, p[j]] = -mu[j, p[rk + i]]
      end
    end
    if isright
      X[p[rk + i], i] = one(R)
    else
      X[i, p[rk + i]] = one(R)
    end
  end
  return true, sol, X
end

function _can_solve_internal(A::MatElem{T}, b::Vector{T}, task::Symbol; side::Symbol = :right) where T <: FieldElement
  if side !== :right && side !== :left
    throw(ArgumentError("Unsupported argument :$side for side: Must be :left or :right."))
  end

  isright = side === :right

  if isright
    nrows(A) != length(b) && error("Incompatible matrices")
    B = matrix(base_ring(A), nrows(A), 1, b)
  else # side == :left
    ncols(A) != length(b) && error("Incompatible matrices")
    B = matrix(base_ring(A), 1, ncols(A), b)
  end
  fl, sol, K = _can_solve_internal(A, B, task, side = side)
  if isright
    x = eltype(b)[ sol[i, 1] for i in 1:nrows(sol) ]
  else # side == :left
    x = eltype(b)[ sol[1, i] for i in 1:ncols(sol) ]
  end
  return fl, x, K
end

# Builds the matrix hcat(transpose(A), transpose(B))
function _hcat_transposed(A::MatElem{T}, B::MatElem{T}) where T
  @assert ncols(A) == ncols(B)
  C = zero_matrix(base_ring(A), ncols(A), nrows(A) + nrows(B))
  for i = 1:ncols(A)
    for j = 1:nrows(A)
      C[i, j] = A[j, i]
    end
    for j = 1:nrows(B)
      C[i, nrows(A) + j] = B[j, i]
    end
  end
  return C
end

# A is supposed to be in rref of rank r
# Return a Vector of length ncols(A) with the first r entries the pivot columns
# of A and the following entries the non-pivot columns (in ascending order).
function pivot_and_non_pivot_cols(A::MatElem, r::Int)
  p = zeros(Int, ncols(A))
  j = 1
  k = 1
  for i = 1:r
    while is_zero_entry(A, i, j)
      p[r + k] = j
      j += 1
      k += 1
    end
    p[i] = j
    j += 1
  end
  while k <= ncols(A) - r
    p[r + k] = j
    j += 1
    k += 1
  end

  return p
end

end
