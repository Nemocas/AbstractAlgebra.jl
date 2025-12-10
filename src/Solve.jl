module Solve

using AbstractAlgebra

using AbstractAlgebra.PrettyPrinting: @show_name
using AbstractAlgebra.PrettyPrinting: Dedent
using AbstractAlgebra.PrettyPrinting: Indent
using AbstractAlgebra.PrettyPrinting: pretty
import AbstractAlgebra: can_solve_with_solution
using AbstractAlgebra: _can_solve_with_solution_fflu
using AbstractAlgebra: _can_solve_with_solution_interpolation
using AbstractAlgebra: _solve_fflu_precomp
using AbstractAlgebra: echelon_form!
using AbstractAlgebra: howell_form!
import AbstractAlgebra: kernel
import AbstractAlgebra: matrix

import Base: show

export solve
export solve_init
export solve_context_type
export can_solve
export can_solve_with_solution_and_kernel

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

function AbstractAlgebra.view(M::LazyTransposeMatElem, r::Union{Colon, AbstractVector{Int}}, c::Union{Colon, AbstractVector{Int}})
  return lazy_transpose(view(data(M), c, r))
end

# If one of the arguments is an Int, then the result is 1-dimensional so wrapping
# it in a lazy_transpose does not make sense (and does not work, if we don't add
# something like lazy_transpose(::AbstractVector))
function AbstractAlgebra.view(M::LazyTransposeMatElem, r::Int, c::Union{Colon, AbstractVector{Int}})
  return view(data(M), c, r)
end
function AbstractAlgebra.view(M::LazyTransposeMatElem, r::Union{Colon, AbstractVector{Int}}, c::Int)
  return view(data(M), c, r)
end

AbstractAlgebra.base_ring(M::LazyTransposeMatElem) = base_ring(data(M))

Base.zero(M::LazyTransposeMatElem) = lazy_transpose(zero(data(M)))
Base.zero(M::LazyTransposeMatElem, i::Int, j::Int) = lazy_transpose(zero(data(M), j, i))

Base.similar(M::LazyTransposeMatElem) = lazy_transpose(similar(data(M)))
Base.similar(M::LazyTransposeMatElem, i::Int, j::Int) = lazy_transpose(similar(data(M), j, i))

################################################################################
#
#  Type traits for the matrix normal forms
#
################################################################################

abstract type MatrixNormalFormTrait end

struct HowellFormTrait <: MatrixNormalFormTrait end # Howell form for PIR
struct HermiteFormTrait <: MatrixNormalFormTrait end # Hermite form for PID
struct RREFTrait <: MatrixNormalFormTrait end # Row-reduced echelon form for fields
struct LUTrait <: MatrixNormalFormTrait end # LU factoring for fields
struct FFLUTrait <: MatrixNormalFormTrait end # "fraction free" LU factoring for fraction fields
struct MatrixInterpolateTrait <: MatrixNormalFormTrait end # interpolate in fraction fields of polynomial rings

function matrix_normal_form_type(R::Ring)
  if is_domain_type(elem_type(R))
    return HermiteFormTrait()
  else
    return HowellFormTrait()
  end
end

matrix_normal_form_type(::Field) = RREFTrait()

# The fflu approach is the fastest over a fraction field (see benchmarks on PR 661)
matrix_normal_form_type(::FracField) = FFLUTrait()
matrix_normal_form_type(::AbstractAlgebra.Rationals{BigInt}) = FFLUTrait()
matrix_normal_form_type(::FracField{T}) where {T <: PolyRingElem} = MatrixInterpolateTrait()

matrix_normal_form_type(A::MatElem) = matrix_normal_form_type(base_ring(A))

################################################################################
#
#  Linear solving context object
#
################################################################################

# T: element type of the base ring
# NFTrait: the used MatrixNormalFormTrait
# MatT: type of the matrix
# RedMatT: type of the reduced/canonical form of the matrix
#          Usually RedMatT == MatT, but not for fraction fields where we do fflu
# TranspRedMatT: type of the reduced/canonical form of the matrix
#          Usually TranspRedMatT == LazyTransposeMatElem{T, RedMatT}, but not
#          for, e.g. Flint types
@attributes mutable struct SolveCtx{T, NFTrait, MatT, RedMatT, TranspRedMatT}
  A::MatT # matrix giving the linear system
  red::RedMatT # reduced/canonical form of A (rref, hnf, lu)
  red_transp::TranspRedMatT # reduced/canonical form of transpose(A)

  # Used for rref / hnf
  trafo::RedMatT # transformation: trafo*A == red
  trafo_transp::TranspRedMatT # transformation: trafo_transp*transpose(A) == red_transp

  # Used for rref
  pivots::Vector{Int} # pivot and non-pivot columns of red
  pivots_transp::Vector{Int} # pivot and non-pivot columns of red_transp

  # Used for lu / fflu
  lu_perm::Generic.Perm{Int} # permutation used for the lu factorization of A
  lu_perm_transp::Generic.Perm{Int} # permutation used for the lu factorization of transpose(A)
  permuted_matrix::MatT # precomputed data for solving
  permuted_matrix_transp::MatT # precomputed data for solving

  # Used for fflu
  scaling_factor::T # factor by which any solution needs to be multiplied
  scaling_factor_transp::T # factor by which any solution needs to be multiplied

  rank::Int # rank of A

  kernel_left::MatT
  kernel_right::MatT

  function SolveCtx{T, NFTrait, MatT, RedMatT, TranspRedMatT}(A::MatT) where {T, NFTrait <: MatrixNormalFormTrait, MatT <: MatElem{T}, RedMatT <: MatElem, TranspRedMatT <: MatElem}
    z = new{T, NFTrait, MatT, RedMatT, TranspRedMatT}()
    z.A = A
    z.rank = -1 # not known yet
    return z
  end
end

function Base.show(io::IO, ::MIME"text/plain", C::SolveCtx)
  io = pretty(io)
  println(io, "Linear solving context of matrix")
  print(io, Indent())
  show(io, MIME"text/plain"(), matrix(C))
  print(io, Dedent())
end

function Base.show(io::IO, C::SolveCtx)
  @show_name(io, C)

  print(io, "Linear solving context")
end

@doc raw"""
    solve_init(A::MatElem)

Return a context object `C` that allows to efficiently solve linear systems
$Ax = b$ or $xA = b$ for different $b$.

# Example

```jldoctest
julia> A = QQ[1 2 3; 0 3 0; 5 0 0];

julia> C = solve_init(A)
Linear solving context of matrix
  [1//1   2//1   3//1]
  [0//1   3//1   0//1]
  [5//1   0//1   0//1]

julia> solve(C, [QQ(1), QQ(1), QQ(1)]; side = :left)
3-element Vector{Rational{BigInt}}:
 1//3
 1//9
 2//15

julia> solve(C, [QQ(1), QQ(1), QQ(1)]; side = :right)
3-element Vector{Rational{BigInt}}:
 1//5
 1//3
 2//45
```
"""
function solve_init(A::MatElem)
  return solve_context_type(base_ring(A))(A)
end

function solve_init(NF::MatrixNormalFormTrait, A::MatElem)
  return solve_context_type(NF, base_ring(A))(A)
end

# For a ring R, the following signatures of `solve_context_type` need to be
# implemented:
# 1) solve_context_type(R)
# 2) solve_context_type(::MatrixNormalFormTrait, elem_type(R))
# Version 1 should pick a matrix_normal_form_type and call 2

function solve_context_type(R::NCRing)
  return solve_context_type(matrix_normal_form_type(R), elem_type(R))
end

function solve_context_type(K::Field)
  # matrix_normal_form_type(K) would be RREFTrait, but we want to use
  # LU in solve contexts
  return solve_context_type(LUTrait(), elem_type(K))
end

function solve_context_type(K::Union{AbstractAlgebra.Rationals{BigInt}, FracField})
  # In this case, we use FFLU
  return solve_context_type(FFLUTrait(), elem_type(K))
end

function solve_context_type(A::MatElem)
  return solve_context_type(base_ring(A))
end

function solve_context_type(NF::MatrixNormalFormTrait, ::Type{T}) where {T <: NCRingElement}
  MatType = dense_matrix_type(T)
  return SolveCtx{T, typeof(NF), MatType, MatType, LazyTransposeMatElem{T, MatType}}
end

function solve_context_type(::FFLUTrait, ::Type{T}) where {T <: NCRingElement}
  # We assume that the ring in question is a fraction field and have to get the
  # type of "integral" matrices, that is, matrices over the base ring of this
  # fraction field
  # Whether this works/makes sense is up to the user
  IntMatT = dense_matrix_type(base_ring_type(T))
  return SolveCtx{T, FFLUTrait, dense_matrix_type(T), IntMatT, IntMatT}
end

solve_context_type(NF::MatrixNormalFormTrait, ::T) where {T <: NCRingElement} = solve_context_type(NF, T)
solve_context_type(NF::MatrixNormalFormTrait, ::Type{T}) where {T <: NCRing} = solve_context_type(NF, elem_type(T))
solve_context_type(NF::MatrixNormalFormTrait, ::T) where {T <: NCRing} = solve_context_type(NF, elem_type(T))
solve_context_type(NF::MatrixNormalFormTrait, ::Type{<: MatElem{T}}) where T = solve_context_type(NF, T)
solve_context_type(NF::MatrixNormalFormTrait, ::MatElem{T}) where T = solve_context_type(NF, T)

matrix_normal_form_type(C::SolveCtx{T, NF}) where {T, NF} = NF()

matrix(C::SolveCtx) = C.A

function _init_reduce(C::SolveCtx{T, RREFTrait}) where T
  if isdefined(C, :red) && isdefined(C, :trafo)
    return nothing
  end

  B, U = echelon_form_with_transformation(matrix(C))
  r = nrows(B)
  while r > 0 && is_zero_row(B, r)
    r -= 1
  end
  set_rank!(C, r)
  C.red = B
  C.trafo = U
  return nothing
end

function _init_reduce(C::SolveCtx{T, LUTrait}) where T
  if isdefined(C, :red) && isdefined(C, :lu_perm)
    return nothing
  end

  B = deepcopy(matrix(C))
  p = Perm(1:nrows(B))
  r = lu!(p, B)
  set_rank!(C, r)
  C.red = B
  C.lu_perm = p
  if r < nrows(C)
    pA = p*matrix(C)
    C.permuted_matrix = pA[r + 1:nrows(C), :]
  else
    C.permuted_matrix = zero(matrix(C), 0, ncols(C))
  end
  return nothing
end

function _init_reduce(C::SolveCtx{T, FFLUTrait}) where T
  if isdefined(C, :red)
    return nothing
  end

  A = matrix(C)
  dA = _common_denominator(A)
  Aint = matrix(parent(dA), nrows(A), ncols(A), [numerator(A[i, j]*dA) for i in 1:nrows(A) for j in 1:ncols(A)])
  p = one(SymmetricGroup(nrows(A)))
  r, dLU = fflu!(p, Aint)

  set_rank!(C, r)
  C.lu_perm = p
  d = divexact(dA, base_ring(C)(dLU))
  C.red = Aint
  C.scaling_factor = d
  if r < nrows(A)
    A2 = p*A
    A3 = A2[r + 1:nrows(A), :]
    C.permuted_matrix = A3
  else
    C.permuted_matrix = zero(A, 0, ncols(A))
  end
  return nothing
end

function _init_reduce(C::SolveCtx{T, HermiteFormTrait}) where T
  if isdefined(C, :red) && isdefined(C, :trafo)
    return nothing
  end

  R, U = hermite_form_with_transformation(matrix(C))
  C.red = R
  C.trafo = U
  return nothing
end

# For this type, we store the Howell form of
#   ( A | I_n)
#   ( 0 |  0 )
# in C.red. From this one can recover the Howell form of A with transformation,
# but also additional information for the kernel.
function _init_reduce(C::SolveCtx{T, HowellFormTrait}) where T
  if isdefined(C, :red)
    return nothing
  end
  A = matrix(C)

  B = hcat(A, identity_matrix(A, nrows(A)))
  if nrows(B) < ncols(B)
    B = vcat(B, zero(A, ncols(B) - nrows(B), ncols(B)))
  end

  howell_form!(B)
  C.red = B
  return nothing
end

function _init_reduce_transpose(C::SolveCtx{T, RREFTrait}) where T
  if isdefined(C, :red_transp) && isdefined(C, :trafo_transp)
    return nothing
  end

  B, U = echelon_form_with_transformation(lazy_transpose(matrix(C)))
  r = nrows(B)
  while r > 0 && is_zero_row(B, r)
    r -= 1
  end
  set_rank!(C, r)
  C.red_transp = B
  C.trafo_transp = U
  return nothing
end

function _init_reduce_transpose(C::SolveCtx{T, LUTrait}) where T
  if isdefined(C, :red_transp) && isdefined(C, :lu_perm_transp)
    return nothing
  end

  B = lazy_transpose(deepcopy(matrix(C)))
  p = Perm(1:nrows(B))
  r = lu!(p, B)
  set_rank!(C, r)
  C.red_transp = B
  C.lu_perm_transp = p
  if r < ncols(C)
    Ap = matrix(C)*p
    C.permuted_matrix_transp = Ap[:, r + 1:ncols(C)]
  else
    C.permuted_matrix_transp = zero(matrix(C), nrows(C), 0)
  end
  return nothing
end

function _init_reduce_transpose(C::SolveCtx{T, FFLUTrait}) where T
  if isdefined(C, :red_transp)
    return nothing
  end

  A = matrix(C)
  dA = _common_denominator(A)
  # We transpose A at this point!
  Aint = matrix(parent(dA), ncols(A), nrows(A), [numerator(A[i, j]*dA) for j in 1:ncols(A) for i in 1:nrows(A)])
  p = one(SymmetricGroup(nrows(Aint)))
  r, dLU = fflu!(p, Aint)

  set_rank!(C, r)
  C.lu_perm_transp = p
  d = divexact(dA, base_ring(C)(dLU))
  C.red_transp = Aint
  C.scaling_factor_transp = d
  if r < ncols(A)
    A2 = A*p
    A3 = A2[:, r + 1:ncols(A)]
    C.permuted_matrix_transp = A3
  else
    C.permuted_matrix_transp = zero(A, nrows(A), 0)
  end
  return nothing
end

function _init_reduce_transpose(C::SolveCtx{T, HermiteFormTrait}) where T
  if isdefined(C, :red_transp) && isdefined(C, :trafo_transp)
    return nothing
  end

  R, U = hermite_form_with_transformation(lazy_transpose(matrix(C)))
  C.red_transp = R
  C.trafo_transp = U
  return nothing
end

function _init_reduce_transpose(C::SolveCtx{T, HowellFormTrait}) where T
  if isdefined(C, :red_transp)
    return nothing
  end
  A = matrix(C)

  AT = lazy_transpose(A)
  B = hcat(AT, identity_matrix(AT, nrows(AT)))
  if nrows(B) < ncols(B)
    B = vcat(B, zero(AT, ncols(B) - nrows(B), ncols(B)))
  end

  howell_form!(B)
  C.red_transp = B
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

# Only for FFLU.
# Factor by which any solution needs to be multiplied.
# This is the chosen denominator of matrix(C) divided by the denominator returned
# by fflu!(matrix(C)).
function scaling_factor(C::SolveCtx)
  _init_reduce(C)
  return C.scaling_factor
end

function scaling_factor_of_transpose(C::SolveCtx)
  _init_reduce_transpose(C)
  return C.scaling_factor_transp
end

# Only for LU and FFLU
# Let A = matrix(C).
# Return the matrix (p*A)[rank(A) + 1:nrows(A), :] where p is lu_permutation(C).
function permuted_matrix(C::SolveCtx)
  _init_reduce(C)
  return C.permuted_matrix
end

# Only for LU and FFLU
# Let A = matrix(C).
# Return the matrix (A*p)[:, rank(A) + 1:ncols(A)] where p is lu_permutation_of_transpose(C).
function permuted_matrix_of_transpose(C::SolveCtx)
  _init_reduce_transpose(C)
  return C.permuted_matrix_transp
end

function lu_permutation(C::SolveCtx)
  _init_reduce(C)
  return C.lu_perm
end

function lu_permutation_of_transpose(C::SolveCtx)
  _init_reduce_transpose(C)
  return C.lu_perm_transp
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

function AbstractAlgebra.rank(C::SolveCtx)
  if C.rank < 0
    _init_reduce(C)
  end
  return C.rank
end

AbstractAlgebra.nrows(C::SolveCtx) = nrows(matrix(C))
AbstractAlgebra.ncols(C::SolveCtx) = ncols(matrix(C))
AbstractAlgebra.base_ring(C::SolveCtx) = base_ring(matrix(C))

function pivot_and_non_pivot_cols(C::SolveCtx)
  if !isdefined(C, :pivots)
    R = reduced_matrix(C)
    r = rank(C)
    C.pivots = pivot_and_non_pivot_cols(R, r)
  end
  return C.pivots
end

function pivot_and_non_pivot_cols_of_transpose(C::SolveCtx)
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
    solve(A::MatElem{T}, b::Vector{T}; side::Symbol = :left) where T
    solve(A::MatElem{T}, b::MatElem{T}; side::Symbol = :left) where T
    solve(C::SolveCtx{T}, b::Vector{T}; side::Symbol = :left) where T
    solve(C::SolveCtx{T}, b::MatElem{T}; side::Symbol = :left) where T

Return $x$ of same type as $b$ solving the linear system $xA = b$, if `side == :left`
(default), or $Ax = b$, if `side == :right`.

If no solution exists, an error is raised.

If a context object `C` is supplied, then the above applies for `A = matrix(C)`.

See also [`can_solve_with_solution`](@ref can_solve_with_solution(::Union{MatElem{T}, SolveCtx{T}}, ::Union{Vector{T}, MatElem{T}}) where T).

# Example

```jldoctest
julia> A = QQ[2 0 0;0 3 0;0 0 5]
[2//1   0//1   0//1]
[0//1   3//1   0//1]
[0//1   0//1   5//1]

julia> solve(A, one(A))
[1//2   0//1   0//1]
[0//1   1//3   0//1]
[0//1   0//1   1//5]
```
"""
function solve(A::Union{MatElem{T}, SolveCtx{T}}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :left) where T
  return solve(matrix_normal_form_type(A), A, b; side)
end

function solve(NF::MatrixNormalFormTrait, A::Union{MatElem{T}, SolveCtx{T}}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :left) where T
  fl, x = can_solve_with_solution(NF, A, b; side = side)
  @req fl "Unable to solve linear system"
  return x
end

@doc raw"""
    can_solve(A::MatElem{T}, b::Vector{T}; side::Symbol = :left) where T
    can_solve(A::MatElem{T}, b::MatElem{T}; side::Symbol = :left) where T
    can_solve(C::SolveCtx{T}, b::Vector{T}; side::Symbol = :left) where T
    can_solve(C::SolveCtx{T}, b::MatElem{T}; side::Symbol = :left) where T

Return `true` if the linear system $xA = b$ or $Ax = b$ with `side == :left`
(default) or `side == :right`, respectively, has a solution and `false` otherwise.

If a context object `C` is supplied, then the above applies for `A = matrix(C)`.

See also [`can_solve_with_solution`](@ref can_solve_with_solution(::Union{MatElem{T}, SolveCtx{T}}, ::Union{Vector{T}, MatElem{T}}) where T).

# Example

```jldoctest
julia> A = QQ[2 0 0;0 3 0;0 0 5]
[2//1   0//1   0//1]
[0//1   3//1   0//1]
[0//1   0//1   5//1]

julia> can_solve(A,one(A))
true

julia> A = ZZ[2 0 0;0 3 0;0 0 5]
[2   0   0]
[0   3   0]
[0   0   5]

julia> can_solve(A,one(A))
false
```
"""
function can_solve(A::Union{MatElem{T}, SolveCtx{T}}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :left) where T
  return can_solve(matrix_normal_form_type(A), A, b; side)
end

function can_solve(NF::MatrixNormalFormTrait, A::Union{MatElem{T}, SolveCtx{T}}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :left) where T
  return _can_solve_internal(NF, A, b, :only_check; side = side)[1]
end

@doc raw"""
    can_solve_with_solution(A::MatElem{T}, b::Vector{T}; side::Symbol = :left) where T
    can_solve_with_solution(A::MatElem{T}, b::MatElem{T}; side::Symbol = :left) where T
    can_solve_with_solution(C::SolveCtx{T}, b::Vector{T}; side::Symbol = :left) where T
    can_solve_with_solution(C::SolveCtx{T}, b::MatElem{T}; side::Symbol = :left) where T

Return `true` and $x$ of same type as $b$ solving the linear system $xA = b$, if
such a solution exists. Return `false` and an empty vector or matrix, if the
system has no solution.

If `side == :right`, the system $Ax = b$ is solved.

If a context object `C` is supplied, then the above applies for `A = matrix(C)`.

See also [`solve`](@ref solve(::Union{MatElem{T}, SolveCtx{T}}, ::Union{Vector{T}, MatElem{T}}) where T).
"""
function can_solve_with_solution(A::Union{MatElem{T}, SolveCtx{T}}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :left) where T
  return can_solve_with_solution(matrix_normal_form_type(A), A, b; side)
end

function can_solve_with_solution(NF::MatrixNormalFormTrait, A::Union{MatElem{T}, SolveCtx{T}}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :left) where T
  return _can_solve_internal(NF, A, b, :with_solution; side = side)[1:2]
end

@doc raw"""
    can_solve_with_solution_and_kernel(A::MatElem{T}, b::Vector{T}; side::Symbol = :left) where T
    can_solve_with_solution_and_kernel(A::MatElem{T}, b::MatElem{T}; side::Symbol = :left) where T
    can_solve_with_solution_and_kernel(C::SolveCtx{T}, b::Vector{T}; side::Symbol = :left) where T
    can_solve_with_solution_and_kernel(C::SolveCtx{T}, b::MatElem{T}; side::Symbol = :left) where T

Return `true`, $x$ of same type as $b$ solving the linear system $xA = b$,
together with a matrix $K$ giving the kernel of $A$ (i.e. $KA = 0$), if such
a solution exists. Return `false`, an empty vector or matrix and an empty matrix,
if the system has no solution.

If `side == :right`, the system $Ax = b$ is solved.

If a context object `C` is supplied, then the above applies for `A = matrix(C)`.

See also [`solve`](@ref solve(::Union{MatElem{T}, SolveCtx{T}}, ::Union{Vector{T}, MatElem{T}}) where T) and [`kernel`](@ref kernel(::Union{MatElem, SolveCtx})).
"""
function can_solve_with_solution_and_kernel(A::Union{MatElem{T}, SolveCtx{T}}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :left) where T
  return can_solve_with_solution_and_kernel(matrix_normal_form_type(A), A, b; side)
end

function can_solve_with_solution_and_kernel(NF::MatrixNormalFormTrait, A::Union{MatElem{T}, SolveCtx{T}}, b::Union{Vector{T}, MatElem{T}}; side::Symbol = :left) where T
  return _can_solve_internal(NF, A, b, :with_kernel; side = side)
end

@doc raw"""
    kernel(A::MatElem; side::Symbol = :left)
    kernel(C::SolveCtx; side::Symbol = :left)

Return a matrix $K$ whose rows generate the left kernel of $A$, that
is, $KA$ is the zero matrix.

If `side == :right`, the columns of $K$ generate the right kernel of $A$, that
is, $AK$ is the zero matrix.

If the base ring is a principal ideal domain, the rows or columns respectively of $K$
are a basis of the respective kernel.

If a context object `C` is supplied, then the above applies for `A = matrix(C)`.

#Example

```jldoctest
julia> A = QQ[2 6 0 0;1 3 0 0;0 0 5 0];

julia> kernel(A, side=:right)
[-3//1   0//1]
[ 1//1   0//1]
[ 0//1   0//1]
[ 0//1   1//1]

julia> kernel(A)
[-1//2   1//1   0//1]
```
"""
function kernel(A::Union{MatElem, SolveCtx}; side::Symbol = :left)
  return kernel(matrix_normal_form_type(A), A; side)
end

################################################################################
#
#  Kernel methods
#
################################################################################

# Catch non-matching normal forms in solve context
function kernel(::NF, ::SolveCtx{T, NF2}; side::Symbol = :left) where {T, NF <: MatrixNormalFormTrait, NF2 <: MatrixNormalFormTrait}
  error("Cannot use normal form of type $NF with solve context of type $NF2")
end

### HowellFormTrait

function kernel(::HowellFormTrait, A::MatElem; side::Symbol = :left)
  check_option(side, [:right, :left], "side")

  if side === :right
    AT = lazy_transpose(A)
    B = hcat(AT, identity_matrix(AT, nrows(AT)))
    if nrows(B) < ncols(B)
      B = vcat(B, zero(AT, ncols(B) - nrows(B), ncols(B)))
    end

    howell_form!(B)
    return lazy_transpose(_kernel_of_howell_form(B, nrows(A), ncols(A)))
  else
    B = hcat(A, identity_matrix(A, nrows(A)))
    if nrows(B) < ncols(B)
      B = vcat(B, zero(A, ncols(B) - nrows(B), ncols(B)))
    end

    howell_form!(B)
    return _kernel_of_howell_form(B, ncols(A), nrows(A))
  end
end

function kernel(::HowellFormTrait, C::SolveCtx{T, HowellFormTrait}; side::Symbol = :left) where T
  check_option(side, [:right, :left], "side")

  if side === :right
    if !isdefined(C, :kernel_right)
      B = reduced_matrix_of_transpose(C)
      C.kernel_right = lazy_transpose(_kernel_of_howell_form(B, nrows(C), ncols(C)))
    end
    return C.kernel_right
  else
    if !isdefined(C, :kernel_left)
      B = reduced_matrix(C)
      C.kernel_left = _kernel_of_howell_form(B, ncols(C), nrows(C))
    end
    return C.kernel_left
  end
end

### HermiteFormTrait

function kernel(::HermiteFormTrait, A::MatElem; side::Symbol = :left)
  check_option(side, [:right, :left], "side")

  if side === :right
    A = hnf(A)
    HH, UU = hermite_form_with_transformation(lazy_transpose(A))
    return lazy_transpose(_kernel_of_hnf(HH, UU))
  else
    A = lazy_transpose(hnf(lazy_transpose(A)))
    H, U = hermite_form_with_transformation(A)
    return _kernel_of_hnf(H, U)
  end
end

function kernel(::HermiteFormTrait, C::SolveCtx{T, HermiteFormTrait}; side::Symbol = :left) where T
  check_option(side, [:right, :left], "side")

  if side === :right
    if !isdefined(C, :kernel_right)
      C.kernel_right = lazy_transpose(_kernel_of_hnf(reduced_matrix_of_transpose(C), transformation_matrix_of_transpose(C)))
    end
    return C.kernel_right
  else
    if !isdefined(C, :kernel_left)
      C.kernel_left = _kernel_of_hnf(reduced_matrix(C), transformation_matrix(C))
    end
    return C.kernel_left
  end
end

### RREFTrait

function kernel(::RREFTrait, A::MatElem; side::Symbol = :left)
  check_option(side, [:right, :left], "side")

  if side === :left
    KK = kernel(RREFTrait(), lazy_transpose(A); side = :right)
    return lazy_transpose(KK)::typeof(A)
  end

  n, K = AbstractAlgebra.nullspace(A)
  if ncols(K) > n
    # For compatibility with `nullspace` methods in Nemo which add zero columns
    K = sub(K, 1:nrows(K), 1:n)
  end
  return K
end

# This does not really use the solve context, so the normal form type of C is
# not important
function kernel(::RREFTrait, C::SolveCtx{T, <: MatrixNormalFormTrait}; side::Symbol = :left) where T
  check_option(side, [:right, :left], "side")

  if side === :right
    if !isdefined(C, :kernel_right)
      C.kernel_right = kernel(RREFTrait(), matrix(C); side = :right)
    end
    return C.kernel_right
  else
    if !isdefined(C, :kernel_left)
      C.kernel_left = kernel(RREFTrait(), matrix(C); side = :left)
    end
    return C.kernel_left
  end
end

### LUTrait / FFLUTrait / MatrixNormalFormTrait

# We can't compute a kernel using a LU/FFLU factoring, so we have to fall back to RREF

function kernel(::LUTrait, A::Union{MatElem, SolveCtx{<: NCRingElement, LUTrait}}; side::Symbol = :left)
  return kernel(RREFTrait(), A; side)
end
function kernel(::FFLUTrait, A::Union{MatElem, SolveCtx{<: NCRingElement, FFLUTrait}}; side::Symbol = :left)
  return kernel(RREFTrait(), A; side)
end
function kernel(::MatrixInterpolateTrait, A::Union{MatElem, SolveCtx{<: NCRingElement, MatrixInterpolateTrait}}; side::Symbol = :left)
  return kernel(RREFTrait(), A; side)
end

################################################################################
#
#  Solving methods (internal)
#
################################################################################

###
# General concept:
# `_can_solve_internal` checks the sanity of the input and then calls
# `_can_solve_internal_no_check` . Only the latter function needs to be
# implemented for a given MatrixNormalFormTrait(). Specifically one needs to implement
# the signature(s)
#   _can_solve_internal_no_check(::NormalFormTrait, A::MatrixType, b::MatrixType, task::Symbol; side::Symbol)
#   _can_solve_internal_no_check(::NormalFormTrait, C::SolveCtx, b::MatrixType, task::Symbol; side::Symbol)
# Inside these functions one can assume that A (resp. C) and b have compatible
# dimensions and that `task` and `side` are set to a "legal" option.
# These functions should then (try to) solve Ax = b (side == :right) or xA = b
# (side == :left) possibly with kernel.
# They must always return a tuple (Bool, MatrixType, MatrixType).
# task may be:
# * :only_check -> It is only tested whether there is a solution, the second
#   and third return value are only for type stability
# * :with_solution -> A solution is computed, the last return value is only
#   for type stability
# * :with_kernel -> A solution and the kernel is computed
###

# Transform a right hand side of type Vector into a MatElem and do sanity checks
function _can_solve_internal(::NFTrait, A::Union{MatElem{T}, SolveCtx{T}}, b::Vector{T}, task::Symbol; side::Symbol = :left) where {T, NFTrait <: MatrixNormalFormTrait}
  check_option(task, [:only_check, :with_solution, :with_kernel], "task")
  check_option(side, [:right, :left], "side")
  @assert all(x -> parent(x) === base_ring(A), b) "Base rings do not match"

  isright = side === :right

  if isright
    check_linear_system_dim_right(A, b)
    B = matrix(base_ring(A), nrows(A), 1, b)
  else # side == :left
    check_linear_system_dim_left(A, b)
    B = matrix(base_ring(A), 1, ncols(A), b)
  end
  fl, sol, K = _can_solve_internal_no_check(NFTrait(), A, B, task; side = side)
  if isright
    x = eltype(b)[ sol[i, 1] for i in 1:nrows(sol) ]
  else # side == :left
    x = eltype(b)[ sol[1, i] for i in 1:ncols(sol) ]
  end
  return fl, x, K
end

# Do sanity checks and call _can_solve_internal_no_check
function _can_solve_internal(::NFTrait, A::Union{MatElem{T}, SolveCtx{T}}, b::MatElem{T}, task::Symbol; side::Symbol = :left) where {T, NFTrait <: MatrixNormalFormTrait}
  check_option(task, [:only_check, :with_solution, :with_kernel], "task")
  check_option(side, [:right, :left], "side")
  @assert base_ring(A) === base_ring(b) "Base rings do not match"

  if side === :right
    check_linear_system_dim_right(A, b)
  else
    check_linear_system_dim_left(A, b)
  end
  return _can_solve_internal_no_check(NFTrait(), A, b, task; side = side)
end

# Catch non-matching normal forms in solve context
function _can_solve_internal_no_check(::NFTrait, C::SolveCtx{T, NFTrait2}, b::MatElem{T}, task::Symbol; side::Symbol = :left) where {T, NFTrait <: MatrixNormalFormTrait, NFTrait2 <: MatrixNormalFormTrait}
  error("Cannot use normal form of type $NFTrait with solve context of type $NFTrait2")
end

### HowellFormTrait

function _can_solve_internal_no_check(::HowellFormTrait, A::MatElem{T}, b::MatElem{T}, task::Symbol; side::Symbol = :left) where T

  if side === :left
    # For side == :left, we pretend that A and b are transposed
    fl, _sol, _K = _can_solve_internal_no_check(HowellFormTrait(), lazy_transpose(A), lazy_transpose(b), task; side = :right)
    return fl, lazy_transpose(_sol), lazy_transpose(_K)
  end

  AT = lazy_transpose(A)
  B = hcat(AT, identity_matrix(AT, nrows(AT)))
  if nrows(B) < ncols(B)
    B = vcat(B, zero(AT, ncols(B) - nrows(B), ncols(B)))
  end

  howell_form!(B)

  m = max(nrows(AT), ncols(AT))
  H = view(B, 1:m, 1:ncols(AT))
  U = view(B, 1:m, ncols(AT) + 1:ncols(B))

  fl, sol = _can_solve_with_hnf(b, H, U, task)
  if !fl || task !== :with_kernel
    return fl, sol, zero(A, 0, 0)
  end

  N = _kernel_of_howell_form(B, nrows(A), ncols(A))
  return true, sol, lazy_transpose(N)
end

function _can_solve_internal_no_check(::HowellFormTrait, C::SolveCtx{T, HowellFormTrait}, b::MatElem{T}, task::Symbol; side::Symbol = :left) where T
  if side === :right
    A = matrix(C)
    B = reduced_matrix_of_transpose(C)
    m = max(ncols(A), nrows(A))
    H = view(B, 1:m, 1:nrows(A))
    U = view(B, 1:m, nrows(A) + 1:ncols(B))

    fl, sol = _can_solve_with_hnf(b, H, U, task)
    if !fl || task !== :with_kernel
      return fl, sol, zero(b, 0, 0)
    end

    return fl, sol, kernel(C; side = :right)
  else# side === :left
    A = matrix(C)
    B = reduced_matrix(C)
    m = max(ncols(A), nrows(A))
    H = view(B, 1:m, 1:ncols(A))
    U = view(B, 1:m, ncols(A) + 1:ncols(B))

    fl, sol_transp = _can_solve_with_hnf(lazy_transpose(b), H, U, task)
    sol = lazy_transpose(sol_transp)
    if !fl || task !== :with_kernel
      return fl, sol, zero(b, 0, 0)
    end

    return fl, sol, kernel(C; side = :left)
  end
end

### HermiteFormTrait

function _can_solve_internal_no_check(::HermiteFormTrait, A::MatElem{T}, b::MatElem{T}, task::Symbol; side::Symbol = :left) where T
  if side === :left
    # For side == :left, we pretend that A and b are transposed
    fl, _sol, _K = _can_solve_internal_no_check(HermiteFormTrait(), lazy_transpose(A), lazy_transpose(b), task; side = :right)
    return fl, lazy_transpose(_sol), lazy_transpose(_K)
  end

  H, S = hermite_form_with_transformation(lazy_transpose(A))
  fl, sol = _can_solve_with_hnf(b, H, S, task)
  if !fl || task !== :with_kernel
    return fl, sol, zero(A, 0, 0)
  end

  N = _kernel_of_hnf(H, S)
  return true, sol, lazy_transpose(N)
end

function _can_solve_internal_no_check(::HermiteFormTrait, C::SolveCtx{T, HermiteFormTrait}, b::MatElem{T}, task::Symbol; side::Symbol = :left) where T
  if side === :right
    fl, sol = _can_solve_with_hnf(b, reduced_matrix_of_transpose(C), transformation_matrix_of_transpose(C), task)
  else
    fl, _sol = _can_solve_with_hnf(lazy_transpose(b), reduced_matrix(C), transformation_matrix(C), task)
    sol = lazy_transpose(_sol)
  end
  if !fl || task !== :with_kernel
    return fl, sol, zero(b, 0, 0)
  end

  return true, sol, kernel(C; side = side)
end

### RREFTrait

function _can_solve_internal_no_check(::RREFTrait, A::MatElem{T}, b::MatElem{T}, task::Symbol; side::Symbol = :left) where T
  if side === :left
    # For side == :left, we pretend that A and b are transposed
    fl, _sol, _K = _can_solve_internal_no_check(RREFTrait(), lazy_transpose(A), lazy_transpose(b), task; side = :right)
    return fl, lazy_transpose(_sol), lazy_transpose(_K)
  end

  mu = hcat(deepcopy(A), deepcopy(b))

  rk = echelon_form!(mu)
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
    X[p[rk + i], i] = one(base_ring(A))
  end

  return true, sol, X
end

function _can_solve_internal_no_check(::RREFTrait, C::SolveCtx{T, RREFTrait}, b::MatElem{T}, task::Symbol; side::Symbol = :left) where T
  if side === :right
    fl, sol = _can_solve_with_rref(b, transformation_matrix(C), rank(C), pivot_and_non_pivot_cols(C), task)
    if !fl || task !== :with_kernel
      return fl, sol, zero(b, 0, 0)
    end

    _, K = _kernel_of_rref(reduced_matrix(C), rank(C), pivot_and_non_pivot_cols(C))
    return fl, sol, K
  else# side === :left
    fl, sol_transp = _can_solve_with_rref(lazy_transpose(b), transformation_matrix_of_transpose(C), rank(C), pivot_and_non_pivot_cols_of_transpose(C), task)
    sol = lazy_transpose(sol_transp)
    if !fl || task !== :with_kernel
      return fl, sol, zero(b, 0, 0)
    end

    _, K = _kernel_of_rref(reduced_matrix_of_transpose(C), rank(C), pivot_and_non_pivot_cols_of_transpose(C))
    return fl, sol, lazy_transpose(K)
  end
end

### LUTrait

# LUTrait with MatElem is not implemented

function _can_solve_internal_no_check(::LUTrait, C::SolveCtx{T, LUTrait}, b::MatElem{T}, task::Symbol; side::Symbol = :left) where T
  if side === :right
    p = lu_permutation(C)
    LU = reduced_matrix(C)
    r = rank(C)
    sol = AbstractAlgebra._solve_lu_precomp(p, LU, b)

    # _solve_lu_precomp only takes care of the first r rows
    # We now need to check whether the remaining rows are fine as well
    n = nrows(C)
    fl = true
    if r < n
      b2 = p*b
      fl = permuted_matrix(C)*sol == view(b2, r + 1:n, :)
    end
  else
    p = lu_permutation_of_transpose(C)
    LU = reduced_matrix_of_transpose(C)
    r = rank(C)
    sol = lazy_transpose(AbstractAlgebra._solve_lu_precomp(p, LU, lazy_transpose(b)))

    # _solve_lu_precomp only takes care of the first r rows
    # We now need to check whether the remaining rows are fine as well
    n = ncols(C)
    fl = true
    if r < n
      b2 = b*p
      fl = sol*permuted_matrix_of_transpose(C) == view(b2, :, r + 1:n)
    end
  end
  if !fl || task !== :with_kernel
    return fl, sol, zero(b, 0, 0)
  end

  return true, sol, kernel(C; side = side)
end

### FFLUTrait / MatrixInterpolateTrait

function _common_denominator(A::MatElem{T}) where T <: Union{FracElem, Rational{BigInt}}
  d = numerator(one(base_ring(A)))
  @inbounds for j in 1:ncols(A)
    for i in 1:nrows(A)
      d = lcm(d, denominator(A[i, j]))
    end
  end
  return d
end

function _can_solve_internal_no_check(NF::Union{FFLUTrait, MatrixInterpolateTrait}, A::MatElem{T}, b::MatElem{T}, task::Symbol; side::Symbol = :left) where T

  if side === :left
    fl, _sol, _K = _can_solve_internal_no_check(NF, lazy_transpose(A), lazy_transpose(b), task; side = :right)
    # This does not return LazyTransposedMat for sol because the matrices are made integral
    return fl, transpose(_sol), lazy_transpose(_K)
  end

  d = lcm(_common_denominator(A), _common_denominator(b))

  Aint = matrix(parent(d), nrows(A), ncols(A), [numerator(A[i, j]*d) for i in 1:nrows(A) for j in 1:ncols(A)])
  bint = matrix(parent(d), nrows(b), ncols(b), [numerator(b[i, j]*d) for i in 1:nrows(b) for j in 1:ncols(b)])

  flag = false
  sol_int = similar(Aint, ncols(Aint), nrows(bint))
  den = one(base_ring(A))
  if NF isa MatrixInterpolateTrait
    try
      flag, sol_int, den = _can_solve_with_solution_interpolation(Aint, bint)
    catch
      flag, sol_int, den = _can_solve_with_solution_fflu(Aint, bint)
    end
  else
    flag, sol_int, den = _can_solve_with_solution_fflu(Aint, bint)
  end

  if !flag
    return flag, zero_matrix(base_ring(A), 0, 0), zero(b, 0, 0)
  end

  sol = change_base_ring(base_ring(A), sol_int)
  sol = divexact(sol, base_ring(A)(den))

  if task === :with_kernel
    # I don't know how to compute the kernel using an (ff)lu factoring
    return flag, sol, kernel(A; side = :right)
  else
    return flag, sol, zero(A, 0, 0)
  end
end

function _can_solve_internal_no_check(::FFLUTrait, C::SolveCtx{T, FFLUTrait}, b::MatElem{T}, task::Symbol; side::Symbol = :left) where T
  # Split up in separate functions to make the compiler happy
  if side === :right
    return _can_solve_internal_no_check_right(FFLUTrait(), C, b, task)
  else
    return _can_solve_internal_no_check_left(FFLUTrait(), C, b, task)
  end
end

function _can_solve_internal_no_check_right(::FFLUTrait, C::SolveCtx{T, FFLUTrait}, b::MatElem{T}, task::Symbol) where T
  K = base_ring(C)
  db = _common_denominator(b)
  bint = matrix(parent(db), nrows(b), ncols(b), [numerator(b[i, j]*db) for i in 1:nrows(b) for j in 1:ncols(b)])
  fl, y_int = _solve_fflu_precomp(lu_permutation(C), reduced_matrix(C), bint)
  if !fl
    return fl, zero(b, 0, 0), zero(b, 0, 0)
  end
  # We have fl == true, but we still have to check whether this really is a solution
  d = scaling_factor(C)
  y = change_base_ring(K, y_int)
  y = y*divexact(d, K(db))
  if rank(C) < nrows(C)
    # We have to check whether y is also a solution for the "lower part"
    # of the system
    pb = lu_permutation(C)*b
    pA = permuted_matrix(C)
    fl = pA*y == pb[rank(C) + 1:nrows(C), :]
  end
  if task === :with_kernel
    # I don't know how to compute the kernel using an (ff)lu factoring
    return fl, y, kernel(RREFTrait(), C; side = :right)
  else
    return fl, y, zero(b, 0, 0)
  end
end

function _can_solve_internal_no_check_left(::FFLUTrait, C::SolveCtx{T, FFLUTrait}, b::MatElem{T}, task::Symbol) where T
  K = base_ring(C)
  db = _common_denominator(b)
  # bint == transpose(b)*db
  bint = matrix(parent(db), ncols(b), nrows(b), [numerator(b[i, j]*db) for j in 1:ncols(b) for i in 1:nrows(b)])
  fl, y_int = _solve_fflu_precomp(lu_permutation_of_transpose(C), reduced_matrix_of_transpose(C), bint)
  if !fl
    return fl, zero(b, 0, 0), zero(b, 0, 0)
  end
  # We have fl == true, but we still have to check whether this really is a solution
  d = scaling_factor_of_transpose(C)
  # transpose y_int
  y = matrix(K, ncols(y_int), nrows(y_int), [ K(y_int[i, j]) for j in 1:ncols(y_int) for i in 1:nrows(y_int)])
  y = y*divexact(d, K(db))
  if rank(C) < ncols(C)
    # We have to check whether y is also a solution for the "right hand part"
    # of the system
    pb = b*lu_permutation_of_transpose(C)
    pA = permuted_matrix_of_transpose(C)
    fl = y*pA == pb[:, rank(C) + 1:ncols(C)]
  end
  if task === :with_kernel
    # I don't know how to compute the kernel using an (ff)lu factoring
    return fl, y, kernel(RREFTrait(), C; side = :left)
  else
    return fl, y, zero(b, 0, 0)
  end
end

# MatrixInterpolateTrait with SolveCtx is not implemented

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

# Compute a matrix N with RN == 0 where the columns of N give a basis for the kernel.
# R must be in rref of rank r and pivots must be of length ncols(R) with the pivot
# columns in the first r entries and the non-pivot columns in the remaining entries.
# Only used by Nemo
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

# Solve Ax = b where H = U*transpose(A) is in Hermite form/Howell form.
# Takes same options for `task` as _can_solve_internal but only returns (flag, solution)
# and no kernel.
function _can_solve_with_hnf(b::MatElem{T}, H::MatElem{T}, U::MatElem{T}, task::Symbol) where T <: RingElement
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
      fl, q = divides(b[k, i], H[j, k])
      if !fl
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

# Compute a matrix N with A*transpose(N) == 0 where the rows of N give a basis for the kernel
# and H = U*transpose(A) is in HNF.
function _kernel_of_hnf(H::MatElem{T}, U::MatElem{T}) where T <: RingElement
  r = nrows(H)
  while r > 0 && is_zero_row(H, r)
    r -= 1
  end
  nullity = nrows(H) - r
  return U[r + 1:r + nullity, 1:nrows(H)]
end

# Compute a matrix N with A*transpose(N) == 0 where the rows of N generate the kernel
# and H is the Howell form of
# (A^t | I_n)
# ( 0  |  0 ).
# nrowsA and ncolsA are supposed to be the number of rows and columns of A.
function _kernel_of_howell_form(H::MatElem{T}, nrowsA::Int, ncolsA::Int) where T <: RingElement
  r = 1
  while r <= nrows(H) && !is_zero_row(H, r)
    r += 1
  end
  r -= 1
  h = view(H, 1:r, 1:nrowsA)
  s = findfirst(i -> is_zero_row(h, i), 1:nrows(h))
  if isnothing(s)
    s = nrows(h)
  else
    s -= 1
  end
  return H[s + 1:nrows(h), nrowsA + 1:nrowsA + ncolsA]
end

################################################################################
#
#  Checks
#
################################################################################

function check_option(x::Symbol, options::Vector{Symbol}, option_name::String, msg::String = "", throw_error::Bool = true)
  if msg == ""
    msg = "Unsupported argument $x for $option_name"
  end
  fl = (x in options)
  if !fl && throw_error
    throw(ArgumentError(msg))
  end
  return fl
end

# Checks whether A and b have the same number of rows
function check_linear_system_dim_right(A::Union{MatElem, SolveCtx}, b::MatElem, throw_error::Bool = true)
  fl = nrows(A) == nrows(b)
  if !fl && throw_error
    error("Incompatible number of rows in linear system (use `side = :left` for a system with the indeterminate on the left)")
  end
  return fl
end

function check_linear_system_dim_right(A::Union{MatElem, SolveCtx}, b::Vector, throw_error::Bool = true)
  fl = nrows(A) == length(b)
  if !fl && throw_error
    error("Incompatible number of rows in linear system (use `side = :left` for a system with the indeterminate on the left)")
  end
  return fl
end

# Checks whether A and b have the same number of columns
function check_linear_system_dim_left(A::Union{MatElem, SolveCtx}, b::MatElem, throw_error::Bool = true)
  fl = ncols(A) == ncols(b)
  if !fl && throw_error
    error("Incompatible number of columns in linear system (use `side = :right` for a system with the indeterminate on the right)")
  end
  return fl
end

function check_linear_system_dim_left(A::Union{MatElem, SolveCtx}, b::Vector, throw_error::Bool = true)
  fl = ncols(A) == length(b)
  if !fl && throw_error
    error("Incompatible number of columns in linear system (use `side = :right` for a system with the indeterminate on the right)")
  end
  return fl
end

end # module Solve
