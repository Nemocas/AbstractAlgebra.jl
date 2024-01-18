module Solve

using AbstractAlgebra

import AbstractAlgebra: base_ring, nrows, ncols, kernel, matrix, rank

################################################################################
#
#  "Lazy" transpose of a matrix
#
################################################################################

mutable struct LazyTransposeMatElem{T, MatT} <: MatElem{T} where {MatT <: MatElem{T}}
  M::MatT
end

data(M::LazyTransposeMatElem) = M.M

# The entries of M and the result are SHARED, so e.g. a setindex! will modify
# 'both' matrices. But this is the point: we don't want to actually transpose
# the matrix.
lazy_transpose(M::MatElem{T}) where T = LazyTransposeMatElem{T, typeof(M)}(M)
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

Base.zero(M::LazyTransposeMatElem) = lazy_transpose(zero(data(M)))
Base.zero(M::LazyTransposeMatElem, i::Int, j::Int) = lazy_transpose(zero(data(M), j, i))

Base.similar(M::LazyTransposeMatElem) = lazy_transpose(similar(data(M)))
Base.similar(M::LazyTransposeMatElem, i::Int, j::Int) = lazy_transpose(similar(data(M), j, i))

################################################################################
#
#  Linear solving context object
#
################################################################################

mutable struct SolveCtx{T, MatT}
  A::MatT # matrix giving the linear system
  red::MatT # rref or HNF of A
  red_transp::LazyTransposeMatElem{T, MatT} # rref or HNF of transpose(A)
  trafo::MatT # transformation: trafo*A == red
  trafo_transp::LazyTransposeMatElem{T, MatT} # transformation: trafo_transp*transpose(A) == red_transp

  rank::Int # rank of A
  pivots::Vector{Int} # pivot and non-pivot columns of red
  pivots_transp::Vector{Int} # pivot and non-pivot columns of red_transp

  function SolveCtx(A::MatElem{T}) where T
    z = new{T, typeof(A)}()
    z.A = A
    z.rank = -1 # not known yet
    return z
  end
end

@doc raw"""
    solve_init(A::MatElem)

Return a context object `C` that allows to efficiently solve linear systems
$Ax = b$ or $xA = b$ for different $b$.
"""
function solve_init(A::MatElem)
  return SolveCtx(A)
end

matrix(C::SolveCtx) = C.A

function _init_reduce(C::SolveCtx{<:FieldElement})
  if isdefined(C, :red) && isdefined(C, :trafo)
    return nothing
  end

  r, R, U = _rref_with_transformation(matrix(C))
  set_rank!(C, r)
  C.red = R
  C.trafo = U
  return nothing
end

function _init_reduce(C::SolveCtx{<:RingElement})
  if isdefined(C, :red) && isdefined(C, :trafo)
    return nothing
  end

  R, U = hnf_with_transform(matrix(C))
  C.red = R
  C.trafo = U
  return nothing
end

function _init_reduce_transpose(C::SolveCtx{<:FieldElement})
  if isdefined(C, :red_transp) && isdefined(C, :trafo_transp)
    return nothing
  end

  r, R, U = _rref_with_transformation(lazy_transpose(matrix(C)))
  set_rank!(C, r)
  C.red_transp = R
  C.trafo_transp = U
  return nothing
end

function _init_reduce_transpose(C::SolveCtx{<:RingElement})
  if isdefined(C, :red_transp) && isdefined(C, :trafo_transp)
    return nothing
  end

  R, U = hnf_with_transform(lazy_transpose(matrix(C)))
  C.red_transp = R
  C.trafo_transp = U
  return nothing
end

function reduced_matrix(C::SolveCtx)
  _init_reduce(C)
  return C.red
end

function reduced_matrix_of_transpose(C::SolveCtx)
  _init_reduce_transpose(C)
  return C.red_transp
end

function transformation_matrix(C::SolveCtx)
  _init_reduce(C)
  return C.trafo
end

function transformation_matrix_of_transpose(C::SolveCtx)
  _init_reduce_transpose(C)
  return C.trafo_transp
end

function set_rank!(C::SolveCtx, r::Int)
  if C.rank >= 0
    @assert C.rank == r
  end
  C.rank = r
  return nothing
end

function AbstractAlgebra.rank(C::SolveCtx{<:FieldElement})
  if C.rank < 0
    _init_reduce(C)
  end
  return C.rank
end

AbstractAlgebra.nrows(C::SolveCtx) = nrows(matrix(C))
AbstractAlgebra.ncols(C::SolveCtx) = ncols(matrix(C))
AbstractAlgebra.base_ring(C::SolveCtx) = base_ring(matrix(C))

function pivot_and_non_pivot_cols(C::SolveCtx{<:FieldElement})
  if !isdefined(C, :pivots)
    R = reduced_matrix(C)
    r = rank(C)
    C.pivots = pivot_and_non_pivot_cols(R, r)
  end
  return C.pivots
end

function pivot_and_non_pivot_cols_of_transpose(C::SolveCtx{<:FieldElement})
  if !isdefined(C, :pivots_transp)
    R = reduced_matrix_of_transpose(C)
    r = rank(C)
    C.pivots_transp = pivot_and_non_pivot_cols(R, r)
  end
  return C.pivots_transp
end

################################################################################
#
#  User facing functions for linear solving
#
################################################################################

@doc raw"""
    solve(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T
    solve(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T
    solve(C::SolveCtx{T}, b::Vector{T}; side::Symbol = :right) where T
    solve(C::SolveCtx{T}, b::MatElem{T}; side::Symbol = :right) where T

Return $x$ of same type as $b$ solving the linear system $Ax = b$, if `side == :right`
(default), or $xA = b$, if `side == :left`.

If no solution exists, an error is raised.

If a context object `C` is supplied, then the above applies for `A = matrix(C)`.

See also [`can_solve_with_solution`](@ref).
"""
function solve(A::Union{MatElem{T}, SolveCtx{T}}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :right) where T
  fl, x = can_solve_with_solution(A, b, side = side)
  fl || throw(ArgumentError("Unable to solve linear system"))
  return x
end

@doc raw"""
    can_solve(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T
    can_solve(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T
    can_solve(A::SolveCtx{T}, b::Vector{T}; side::Symbol = :right) where T
    can_solve(A::SolveCtx{T}, b::MatElem{T}; side::Symbol = :right) where T

Return `true` if the linear system $Ax = b$ or $xA = b$ with `side == :right`
(default) or `side == :left`, respectively, has a solution and `false` otherwise.

If a context object `C` is supplied, then the above applies for `A = matrix(C)`.

See also [`can_solve_with_solution`](@ref).
"""
function can_solve(A::Union{MatElem{T}, SolveCtx{T}}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :right) where T
  return _can_solve_internal(A, b, :only_check; side = side)[1]
end

@doc raw"""
    can_solve_with_solution(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T
    can_solve_with_solution(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T
    can_solve_with_solution(A::SolveCtx{T}, b::Vector{T}; side::Symbol = :right) where T
    can_solve_with_solution(A::SolveCtx{T}, b::MatElem{T}; side::Symbol = :right) where T

Return `true` and $x$ of same type as $b$ solving the linear system $Ax = b$, if
such a solution exists. Return `false` and an empty vector or matrix, if the
system has no solution.

If `side == :left`, the system $xA = b$ is solved.

If a context object `C` is supplied, then the above applies for `A = matrix(C)`.

See also [`solve`](@ref).
"""
function can_solve_with_solution(A::Union{MatElem{T}, SolveCtx{T}}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :right) where T
  return _can_solve_internal(A, b, :with_solution; side = side)[1:2]
end

@doc raw"""
    can_solve_with_solution_and_kernel(A::MatElem{T}, b::Vector{T}; side::Symbol = :right) where T
    can_solve_with_solution_and_kernel(A::MatElem{T}, b::MatElem{T}; side::Symbol = :right) where T
    can_solve_with_solution_and_kernel(A::SolveCtx{T}, b::Vector{T}; side::Symbol = :right) where T
    can_solve_with_solution_and_kernel(A::SolveCtx{T}, b::MatElem{T}; side::Symbol = :right) where T

Return `true`, $x$ of same type as $b$ solving the linear system $Ax = b$,
together with a matrix $K$ giving the kernel of $A$ (i.e. $AK = 0$), if such
a solution exists. Return `false`, an empty vector or matrix and an empty matrix,
if the system has no solution.

If `side == :left`, the system $xA = b$ is solved.

If a context object `C` is supplied, then the above applies for `A = matrix(C)`.

See also [`solve`](@ref) and [`kernel`](@ref).
"""
function can_solve_with_solution_and_kernel(A::Union{MatElem{T}, SolveCtx{T}}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :right) where T
  return _can_solve_internal(A, b, :with_kernel; side = side)
end

function AbstractAlgebra.kernel(C::SolveCtx{<:FieldElement}; side::Symbol = :right)
  if side !== :right && side !== :left
    throw(ArgumentError("Unsupported argument :$side for side: Must be :left or :right."))
  end

  if side === :right
    return _kernel_of_rref(reduced_matrix(C), rank(C), pivot_and_non_pivot_cols(C))
  else
    nullity, X = _kernel_of_rref(reduced_matrix_of_transpose(C), rank(C), pivot_and_non_pivot_cols_of_transpose(C))
    # X is of type LazyTransposeMatElem
    return nullity, data(X)
  end
end

function AbstractAlgebra.kernel(C::SolveCtx{<:RingElement}; side::Symbol = :right)
  if side !== :right && side !== :left
    throw(ArgumentError("Unsupported argument :$side for side: Must be :left or :right."))
  end

  if side === :right
    return _kernel_of_hnf(matrix(C), reduced_matrix_of_transpose(C), transformation_matrix_of_transpose(C))
  else
    nullity, X = _kernel_of_hnf(lazy_transpose(matrix(C)), reduced_matrix(C), transformation_matrix(C))
    # X is of type LazyTransposeMatElem
    return nullity, data(X)
  end
end

################################################################################
#
#  Internal functionality
#
################################################################################

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

# Transform a right hand side of type Vector into a MatElem
function _can_solve_internal(A::Union{MatElem{T}, SolveCtx{T}}, b::Vector{T}, task::Symbol; side::Symbol = :right) where T
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

## _can_solve_internal over FIELDS
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

  R = base_ring(A)

  if side === :left
    # For side == :left, we pretend that A and b are transposed
    fl, sol, K = _can_solve_internal(lazy_transpose(A), lazy_transpose(b), task, side = :right)
    return fl, data(sol), data(K)
  end

  nrows(A) != nrows(b) && error("Incompatible matrices")
  mu = hcat(A, b)

  rk = rref!(mu)
  p = pivot_and_non_pivot_cols(mu, rk)
  if any(i -> i > ncols(A), p[1:rk])
    return false, zero(A, 0, 0), zero(A, 0, 0)
  end
  if task === :only_check
    return true, zero(A, 0, 0), zero(A, 0, 0)
  end

  # Compute a solution
  sol = zero(A, ncols(A), ncols(b))
  for i = 1:rk
    for j = 1:ncols(b)
      sol[p[i], j] = mu[i, ncols(A) + j]
    end
  end
  if task === :with_solution
    return true, sol, zero(A, 0, 0)
  end

  # Build the kernel
  nullity = ncols(A) - rk
  X = zero(A, ncols(A), nullity)
  for i = 1:nullity
    for j = 1:rk
      X[p[j], i] = -mu[j, p[rk + i]]
    end
    X[p[rk + i], i] = one(R)
  end

  return true, sol, X
end

## _can_solve_internal over RINGS
# See _can_solve_internal over fields for more explanation on the arguments and output
function _can_solve_internal(A::MatElem{T}, b::MatElem{T}, task::Symbol; side::Symbol = :right) where T <: RingElement
  if task !== :only_check && task !== :with_solution && task !== :with_kernel
    error("task $(task) not recognized")
  end

  if side !== :right && side !== :left
    throw(ArgumentError("Unsupported argument :$side for side: Must be :left or :right."))
  end

  R = base_ring(A)

  if side === :left
    # For side == :left, we pretend that A and b are transposed
    fl, sol, K = _can_solve_internal(lazy_transpose(A), lazy_transpose(b), task, side = :right)
    return fl, data(sol), data(K)
  end

  nrows(A) != nrows(b) && error("Incompatible matrices")

  H, S = hnf_with_transform(lazy_transpose(A))
  fl, sol = _can_solve_with_hnf(b, H, S, task)
  if !fl || task !== :with_kernel
    return fl, sol, zero(A, 0, 0)
  end

  n, N = _kernel_of_hnf(A, H, S)
  return true, sol, N
end

## _can_solve_internal over FIELDS with SOLVE CONTEXT
# See _can_solve_internal over fields for more explanation on the arguments and output
function _can_solve_internal(C::SolveCtx{T}, b::MatElem{T}, task::Symbol; side::Symbol = :right) where T <: FieldElement
  if task !== :only_check && task !== :with_solution && task !== :with_kernel
    error("task $(task) not recognized")
  end

  if side !== :right && side !== :left
    throw(ArgumentError("Unsupported argument :$side for side: Must be :left or :right."))
  end

  if side === :right
    fl, sol = _can_solve_with_rref(b, transformation_matrix(C), rank(C), pivot_and_non_pivot_cols(C), task)
  else
    fl, sol = _can_solve_with_rref(lazy_transpose(b), transformation_matrix_of_transpose(C), rank(C), pivot_and_non_pivot_cols_of_transpose(C), task)
    sol = data(sol)
  end
  if !fl || task !== :with_kernel
    return fl, sol, zero(b, 0, 0)
  end

  return true, sol, kernel(C, side = side)[2]
end

## _can_solve_internal over RINGS with SOLVE CONTEXT
# See _can_solve_internal over fields for more explanation on the arguments and output
function _can_solve_internal(C::SolveCtx{T}, b::MatElem{T}, task::Symbol; side::Symbol = :right) where T <: RingElement
  if task !== :only_check && task !== :with_solution && task !== :with_kernel
    error("task $(task) not recognized")
  end

  if side !== :right && side !== :left
    throw(ArgumentError("Unsupported argument :$side for side: Must be :left or :right."))
  end

  if side === :right
    fl, sol = _can_solve_with_hnf(b, reduced_matrix_of_transpose(C), transformation_matrix_of_transpose(C), task)
  else
    fl, sol = _can_solve_with_hnf(lazy_transpose(b), reduced_matrix(C), transformation_matrix(C), task)
    sol = data(sol)
  end
  if !fl || task !== :with_kernel
    return fl, sol, zero(b, 0, 0)
  end

  return true, sol, kernel(C, side = side)[2]
end

################################################################################
#
#  Internals for solving of row reduced matrices
#
################################################################################

# Solve Ax = b with U*A in rref of rank r.
# pivots must be of length ncols(A) and contain the pivot columns of U*A in the
# first r entries.
# Takes same options for `task` as _can_solve_internal but only returns (flag, solution)
# and no kernel.
function _can_solve_with_rref(b::MatElem{T}, U::MatElem{T}, r::Int, pivots::Vector{Int}, task::Symbol) where T <: FieldElement
  if task !== :only_check && task !== :with_solution && task !== :with_kernel
    error("task $(task) not recognized")
  end

  nrows(U) != nrows(b) && error("Incompatible matrices")

  bU = U*b
  if any(i -> !is_zero_row(bU, i), r + 1:nrows(bU))
    return false, zero(b, 0, 0)
  end
  if task === :only_check
    return true, zero(b, 0, 0)
  end

  # Compute a solution
  sol = zero(b, length(pivots), ncols(b))
  for i = 1:r
    for j = 1:ncols(b)
      sol[pivots[i], j] = bU[i, j]
    end
  end
  return true, sol
end

# Compute a matrix N with RN == 0 where the columns of N give a basis for the kernel.
# R must be in rref of rank r and pivots must be of length ncols(R) with the pivot
# columns in the first r entries and the non-pivot columns in the remaining entries.
function _kernel_of_rref(R::MatElem{T}, r::Int, pivots::Vector{Int}) where T <: FieldElement
  @assert length(pivots) == ncols(R)
  nullity = ncols(R) - r
  X = zero(R, ncols(R), nullity)
  for i = 1:nullity
    for j = 1:r
      X[pivots[j], i] = -R[j, pivots[r + i]]
    end
    X[pivots[r + i], i] = one(base_ring(R))
  end
  return nullity, X
end

# Solve Ax = b where H = U*transpose(A) is in HNF.
# Takes same options for `task` as _can_solve_internal but only returns (flag, solution)
# and no kernel.
function _can_solve_with_hnf(b::MatElem{T}, H::MatElem{T}, U::MatElem{T}, task::Symbol) where T <: RingElement
  if task !== :only_check && task !== :with_solution && task !== :with_kernel
    error("task $(task) not recognized")
  end

  ncols(H) != nrows(b) && error("Incompatible matrices")

  sol = lazy_transpose(zero(b, nrows(H), ncols(b)))
  l = min(nrows(H), ncols(H))
  b = deepcopy(b)
  for i = 1:ncols(b)
    for j = 1:l
      k = 1
      while k <= ncols(H) && is_zero_entry(H, j, k)
        k += 1
      end
      if k > ncols(H)
        continue
      end
      q, r = divrem(b[k, i], H[j, k])
      if !iszero(r)
        return false, zero(b, 0, 0)
      end
      for h = k:ncols(H)
        b[h, i] -= q*H[j, h]
      end
      sol[i, j] = q
    end
  end
  if !is_zero(b)
    return false, zero(b, 0, 0)
  end
  if task === :only_check
    return true, zero(b, 0, 0)
  end
  return true, lazy_transpose(U)*lazy_transpose(sol)
end

# Compute a matrix N with AN == 0 where the columns of N give a basis for the kernel
# and H = U*transpose(A) is in HNF.
# The matrix A is only needed to get the return type right (MatElem vs LazyTransposeMatElem)
function _kernel_of_hnf(A::MatElem{T}, H::MatElem{T}, U::MatElem{T}) where T <: RingElement
  nullity = nrows(H)
  for i = nrows(H):-1:1
    if !is_zero_row(H, i)
      nullity = nrows(H) - i
      break
    end
  end
  N = zero(A, nrows(H), nullity)
  for i = 1:nrows(N)
    for j = 1:ncols(N)
      N[i, j] = U[nrows(U) - j + 1, i]
    end
  end
  return nullity, N
end

# Copied from Hecke, to be replaced with echelon_form_with_transformation eventually
function _rref_with_transformation(M::MatElem{T}) where T <: FieldElement
  n = hcat(M, identity_matrix(base_ring(M), nrows(M)))
  rref!(n)
  s = nrows(n)
  while s > 0 && iszero(sub(n, s:s, 1:ncols(M)))
    s -= 1
  end
  return s, sub(n, 1:nrows(M), 1:ncols(M)), sub(n, 1:nrows(M), ncols(M)+1:ncols(n))
end

end
