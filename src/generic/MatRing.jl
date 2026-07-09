###############################################################################
#
#   MatRing.jl : Generic nxn matrices over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{MatRingElem{T}}) where T <: NCRingElement = MatRing{T}

elem_type(::Type{MatRing{T}}) where {T <: NCRingElement} = MatRingElem{T}

base_ring(a::MatRing{T}) where {T <: NCRingElement} = a.base_ring::parent_type(T)

base_ring(a::MatRingElem{T}) where {T <: NCRingElement} = base_ring(matrix(a))

@doc raw"""
    parent(a::MatRingElem{T}) where T <: NCRingElement

Return the matrix algebra containing `a`.

This is the matrix algebra over the base ring of `a` whose degree is the
number of rows (equivalently columns) of `a`.

**Examples**

```jldoctest
julia> S = matrix_ring(ZZ, 2)
Matrix ring of degree 2
  over integers

julia> A = S([1 2; 3 4])
[1   2]
[3   4]

julia> parent(A) == S
true
```
"""
parent(a::MatRingElem{T}) where T <: NCRingElement = MatRing{T}(base_ring(a), nrows(matrix(a)))

is_exact_type(::Type{MatRingElem{T}}) where T <: NCRingElement = is_exact_type(T)

is_domain_type(::Type{MatRingElem{T}}) where T <: NCRingElement = false

###############################################################################
#
#   Basic manipulation
#
###############################################################################

number_of_rows(a::MatRing) = a.n

number_of_columns(a::MatRing) = number_of_rows(a)

number_of_rows(a::MatRingElem) = nrows(matrix(a))

number_of_columns(a::MatRingElem) = ncols(matrix(a))

Base.@propagate_inbounds getindex(a::MatRingElem, r::Int, c::Int) = matrix(a)[r, c]

Base.@propagate_inbounds function setindex!(a::MatRingElem, d::NCRingElement,
                                            r::Int, c::Int)
    matrix(a)[r, c] = d
end

Base.isassigned(a::MatRingElem, i, j) = isassigned(matrix(a), i, j)

###############################################################################
#
#   Transpose
#
###############################################################################

transpose(x::MatRingElem) = MatRingElem(transpose(matrix(x)))
transpose!(x::MatRingElem) = MatRingElem(transpose!(matrix(x)))
transpose!(z::T, x::T) where T <: MatRingElem = MatRingElem(transpose!(matrix(z), matrix(x)))

###############################################################################
#
#   Solve
#
###############################################################################

function _can_solve_with_solution_lu(M::MatRingElem{T}, B::MatRingElem{T}) where {T <: RingElement}
   check_parent(M, B)
   R = base_ring(M)
   MS = matrix(M)
   BS = matrix(B)
   flag, S = _can_solve_with_solution_lu(MS, BS)
   SA = MatRingElem(S)
   return flag, SA
end

function AbstractAlgebra.can_solve_with_solution(M::MatRingElem{T}, B::MatRingElem{T}) where {T <: RingElement}
   check_parent(M, B)
   R = base_ring(M)
   MS = matrix(M)
   BS = matrix(B)
   flag, S = can_solve_with_solution(MS, BS)
   SA = MatRingElem(S)
   return flag, SA
end

function _can_solve_with_solution_fflu(M::MatRingElem{T}, B::MatRingElem{T}) where {T <: RingElement}
   check_parent(M, B)
   R = base_ring(M)
   MS = matrix(M)
   BS = matrix(B)
   flag, S, d = _can_solve_with_solution_fflu(MS, BS)
   SA = MatRingElem(S)
   return flag, SA, d
end

###############################################################################
#
#   Minimal polynomial
#
###############################################################################

@doc raw"""
    minpoly(S::Ring, M::MatRingElem{T}) where {T <: RingElement}

Return the minimal polynomial $p$ of the matrix $M$. The polynomial ring $S$
of the resulting polynomial must be supplied and the matrix must be square.
"""
function minpoly(S::Ring, M::MatRingElem{T}, charpoly_only::Bool = false) where {T <: RingElement}
   return minpoly(S, matrix(M), charpoly_only)
end

function minpoly(M::MatRingElem{T}, charpoly_only::Bool = false) where {T <: RingElement}
   R = base_ring(M)
   Rx, x = polynomial_ring(R; cached=false)
   return minpoly(Rx, M, charpoly_only)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function add!(A::T, B::T) where T <: MatRingElem
   return MatRingElem(add!(matrix(A), matrix(B)))
end

function add!(c::T, a::T, b::T) where T <: MatRingElem
  return MatRingElem(add!(matrix(c), matrix(a), matrix(b)))
end

function sub!(A::T, B::T) where T <: MatRingElem
   return MatRingElem(sub!(matrix(A), matrix(B)))
end

function sub!(c::T, a::T, b::T) where T <: MatRingElem
  return MatRingElem(sub!(matrix(c), matrix(a), matrix(b)))
end

function mul!(A::T, B::T) where T <: MatRingElem
   return MatRingElem(mul!(matrix(A), matrix(B)))
end

function mul!(c::T, a::T, b::T) where T <: MatRingElem
  return MatRingElem(mul!(matrix(c), matrix(a), matrix(b)))
end

function mul!(c::MatRingElem{T}, a::MatRingElem{T}, b::T) where T <: NCRingElement
  return MatRingElem(mul!(matrix(c), matrix(a), b))
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{S}, ::Type{S}) where {T <: NCRingElement, S <: MatRingElem{T}} = MatRingElem{T}

function promote_rule(::Type{S}, ::Type{U}) where {T <: NCRingElement, S <: MatRingElem{T}, U <: NCRingElement}
   promote_rule(T, U) == T ? MatRingElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

@doc raw"""
    (a::MatRing{T})() where {T <: NCRingElement}
    (a::MatRing{T})(b::S) where {S <: NCRingElement, T <: NCRingElement}
    (a::MatRing{T})(b::T) where {S <: NCRingElement, T <: MatRingElem{S}}
    (a::MatRing{T})(b::MatRingElem{T}) where {T <: NCRingElement}
    (a::MatRing{T})(b::MatrixElem{S}) where {S <: NCRingElement, T <: NCRingElement}
    (a::MatRing{T})(b::Matrix{S}) where {S <: NCRingElement, T <: NCRingElement}
    (a::MatRing{T})(b::Vector{S}) where {S <: NCRingElement, T <: NCRingElement}

Construct an element of the matrix algebra `a`.

The call `a()` returns the zero matrix in `a`.

If `b` is a ring element coercible into the base ring of `a`, then
`a(b)` returns the scalar matrix in `a` whose diagonal entries are `b`.

If `b` is a matrix algebra element whose parent is `a`, then `a(b)`
returns `b` unchanged.

If `b` is an algebraic matrix whose dimensions and base ring agree with
those of `a`, then `a(b)` returns the corresponding element of `a`. If
necessary, a new matrix with the appropriate implementation type is constructed.

If `b` is a Julia vector or matrix, then `a(b)` constructs an element of
`a` whose entries are obtained by coercing the entries of `b` into the base
ring of `a`. The dimensions of the Julia object must be compatible with
the degree of the matrix algebra.

# Examples

```jldoctest
julia> R, = residue_ring(ZZ, 7);

julia> A = matrix_ring(R, 3);

julia> A()
[0   0   0]
[0   0   0]
[0   0   0]

julia> A(R(3))
[3   0   0]
[0   3   0]
[0   0   3]

julia> M = A([R(1) R(2) R(3); R(4) R(5) R(6); R(0) R(1) R(2)])
[1   2   3]
[4   5   6]
[0   1   2]

julia> A(M) === M
true

julia> A([R(1), R(2), R(3), R(4), R(5), R(6), R(0), R(1), R(2)])
[1   2   3]
[4   5   6]
[0   1   2]

julia> S = matrix_space(R, 3, 3);

julia> N = S([R(1) R(0) R(0); R(0) R(2) R(0); R(0) R(0) R(3)]);

julia> A(N)
[1   0   0]
[0   2   0]
[0   0   3]
```
""" MatRing

function (a::MatRing{T})() where {T <: NCRingElement}
   R = base_ring(a)
   z = MatRingElem(zero_matrix(R, a.n, a.n))
   return z
end

function (a::MatRing{T})(b::S) where {S <: NCRingElement, T <: NCRingElement}
   R = base_ring(a)
   z = MatRingElem(scalar_matrix(R, a.n, R(b)))
   return z
end

# to resolve ambiguity for MatRing{MatRing{...}}
function (a::MatRing{T})(b::T) where {S <: NCRingElement, T <: MatRingElem{S}}
   R = base_ring(a)
   z = MatRingElem(scalar_matrix(R, a.n, R(b)))
   return z
end

function (a::MatRing{T})(b::MatRingElem{T}) where {T <: NCRingElement}
   parent(b) != a && error("Unable to coerce matrix")
   return b
end

function (a::MatRing{T})(b::MatrixElem{S}) where {S <: NCRingElement, T <: NCRingElement}
   R = base_ring(a)
   _check_dim(nrows(a), ncols(a), b)
   z = MatRingElem(matrix(R, b))
   return z
end

function (a::MatRing{T})(b::Matrix{S}) where {S <: NCRingElement, T <: NCRingElement}
   _check_dim(a.n, a.n, b)
   R = base_ring(a)
   z = MatRingElem(matrix(R, b))
   return z
end

function (a::MatRing{T})(b::Vector{S}) where {S <: NCRingElement, T <: NCRingElement}
  _check_dim(a.n, a.n, b)
   R = base_ring(a)
   z = MatRingElem(R, a.n, b)
   return z
end

###############################################################################
#
#   MatRing constructor
#
###############################################################################

function matrix_ring(R::AbstractAlgebra.NCRing, n::Int; cached::Bool = true)
   # TODO: the 'cached' argument is ignored and mainly here for backwards compatibility
   # (and perhaps future compatibility, in case we need it again)
   @req n >= 0 "n must be a non-negative integer"
   T = elem_type(R)
   return MatRing{T}(R, n)
end
