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

Return the parent object of the given matrix.
"""
parent(a::MatRingElem{T}) where T <: NCRingElement = MatRing{T}(base_ring(a), nrows(a.data))

is_exact_type(::Type{MatRingElem{T}}) where T <: NCRingElement = is_exact_type(T)

is_domain_type(::Type{MatRingElem{T}}) where T <: NCRingElement = false

###############################################################################
#
#   Basic manipulation
#
###############################################################################

number_of_rows(a::MatRingElem) = nrows(matrix(a))

number_of_columns(a::MatRingElem) = ncols(matrix(a))

Base.@propagate_inbounds getindex(a::MatRingElem, r::Int, c::Int) = a.data[r, c]

Base.@propagate_inbounds function setindex!(a::MatRingElem, d::NCRingElement,
                                            r::Int, c::Int)
    a.data[r, c] = base_ring(a)(d)
end

Base.isassigned(a::MatRingElem, i, j) = isassigned(a.data, i, j)

###############################################################################
#
#   Transpose
#
###############################################################################

function transpose(x::MatRingElem)
   return MatRingElem(transpose(matrix(x)))
end

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

function add!(A::MatRingElem{T}, B::MatRingElem{T}) where T <: NCRingElement
   #=A.data ==# add!(A.data, B.data)  ### !!struct is NOT mutable!!
   return A
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
   z = MatRingElem(matrix(R, a.n, a.n, b))
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
