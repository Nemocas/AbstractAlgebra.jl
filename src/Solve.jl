module Solve

using AbstractAlgebra

import AbstractAlgebra: base_ring, nrows, ncols

################################################################################
#
#  "Lazy" transpose of a matrix (very experimental)
#
################################################################################

mutable struct LazyTransposeMatElem{T} <: MatElem{T}
  M::MatElem{T}
end

data(M::LazyTransposeMatElem) = M.M

# The entries of M and the result are SHARED, so e.g. a setindex! will modify
# 'both' matrices. But this is the point: we don't want to actually transpose
# the matrix.
lazy_transpose(M::MatElem) = LazyTransposeMatElem(M)
lazy_transpose(M::LazyTransposeMatElem) = data(M)

# Change the order of rows and columns in nrows, ncols, getindex and setindex!
AbstractAlgebra.nrows(M::LazyTransposeMatElem) = ncols(data(M))
AbstractAlgebra.ncols(M::LazyTransposeMatElem) = nrows(data(M))

Base.getindex(M::LazyTransposeMatElem, r::Int, c::Int) = data(M)[c, r]
function Base.setindex!(M::LazyTransposeMatElem{T}, d::T, r::Int, c::Int) where T
  setindex!(M.M, d, c, r)
  return M
end

AbstractAlgebra.base_ring(M::LazyTransposeMatElem) = base_ring(data(M))

################################################################################
#
#  User facing functions for linear solving
#
################################################################################

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

################################################################################
#
#  Internal functionality
#
################################################################################

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

  isleft = side === :left

  R = base_ring(A)

  if isleft
    # For side == :left, we pretend that A and b are transposed
    A = lazy_transpose(A)
    b = lazy_transpose(b)
  end

  nrows(A) != nrows(b) && error("Incompatible matrices")
  mu = hcat(A, b)

  rk = rref!(mu)
  p = pivot_and_non_pivot_cols(mu, rk)
  if any(i -> i > ncols(A), p[1:rk])
    return false, zero_matrix(R, 0, 0), zero_matrix(R, 0, 0)
  end
  if task === :only_check
    return true, zero_matrix(R, 0, 0), zero_matrix(R, 0, 0)
  end

  # Compute a solution
  if !isleft
    sol = zero_matrix(R, ncols(A), ncols(b))
  else
    # Pretend again that sol is transposed, so that we can use the same code to
    # fill it
    sol = lazy_transpose(zero_matrix(R, ncols(b), ncols(A)))
  end
  for i = 1:rk
    for j = 1:ncols(b)
      sol[p[i], j] = mu[i, ncols(A) + j]
    end
  end
  if isleft
    # Transpose back
    sol = data(sol)
  end
  if task === :with_solution
    return true, sol, zero_matrix(R, 0, 0)
  end

  # Build the kernel
  nullity = ncols(A) - rk
  if !isleft
    X = zero_matrix(R, ncols(A), nullity)
  else
    X = lazy_transpose(zero_matrix(R, nullity, ncols(A)))
  end
  for i = 1:nullity
    for j = 1:rk
      X[p[j], i] = -mu[j, p[rk + i]]
    end
    X[p[rk + i], i] = one(R)
  end
  if isleft
    # Transpose back
    X = data(X)
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
