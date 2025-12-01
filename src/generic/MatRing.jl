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

@doc raw"""
    parent(a::MatRingElem{T}) where T <: NCRingElement

Return the parent object of the given matrix.
"""
parent(a::MatRingElem{T}) where T <: NCRingElement = MatRing{T}(base_ring(a), size(a.entries)[1])

is_exact_type(::Type{MatRingElem{T}}) where T <: NCRingElement = is_exact_type(T)

is_domain_type(::Type{MatRingElem{T}}) where T <: NCRingElement = false

###############################################################################
#
#   Transpose
#
###############################################################################

function transpose(x::MatRingElem{T}) where T <: NCRingElement
   arr = permutedims(x.entries)
   z = MatRingElem{T}(base_ring(x), arr)
   return z
end

###############################################################################
#
#   Solve
#
###############################################################################

function _can_solve_with_solution_lu(M::MatRingElem{T}, B::MatRingElem{T}) where {T <: RingElement}
   check_parent(M, B)
   R = base_ring(M)
   MS = MatSpaceElem{T}(R, M.entries) # convert to ordinary matrix
   BS = MatSpaceElem{T}(R, B.entries)
   flag, S = _can_solve_with_solution_lu(MS, BS)
   SA = MatRingElem{T}(R, S.entries)
   return flag, SA
end

function AbstractAlgebra.can_solve_with_solution(M::MatRingElem{T}, B::MatRingElem{T}) where {T <: RingElement}
   check_parent(M, B)
   R = base_ring(M)
   # TODO: Once #1955 is resolved, the conversion to matrix and back to MatRingElem
   # should be done better
   MS = matrix(R, M.entries) # convert to ordinary matrix
   BS = matrix(R, B.entries)
   flag, S = can_solve_with_solution(MS, BS)
   SA = MatRingElem{T}(R, Array(S))
   return flag, SA
end

function _can_solve_with_solution_fflu(M::MatRingElem{T}, B::MatRingElem{T}) where {T <: RingElement}
   check_parent(M, B)
   R = base_ring(M)
   MS = MatSpaceElem{T}(R, M.entries) # convert to ordinary matrix
   BS = MatSpaceElem{T}(R, B.entries)
   flag, S, d = _can_solve_with_solution_fflu(MS, BS)
   SA = MatRingElem{T}(R, S.entries)
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
   MS = MatSpaceElem{T}(base_ring(M), M.entries) # convert to ordinary matrix
   return minpoly(S, MS, charpoly_only)
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
   A.entries .+= B.entries
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
   entries = Matrix{T}(undef, a.n, a.n)
   for i = 1:a.n
      for j = 1:a.n
         entries[i, j] = zero(R)
      end
   end
   z = MatRingElem{T}(R, entries)
   return z
end

function (a::MatRing{T})(b::S) where {S <: NCRingElement, T <: NCRingElement}
   R = base_ring(a)
   entries = Matrix{T}(undef, a.n, a.n)
   rb = R(b)
   for i = 1:a.n
      for j = 1:a.n
         if i != j
            entries[i, j] = zero(R)
         else
            entries[i, j] = rb
         end
      end
   end
   z = MatRingElem{T}(R, entries)
   return z
end

# to resolve ambiguity for MatRing{MatRing{...}}
function (a::MatRing{T})(b::T) where {S <: NCRingElement, T <: MatRingElem{S}}
   R = base_ring(a)
   entries = Matrix{T}(undef, a.n, a.n)
   rb = R(b)
   for i = 1:a.n
      for j = 1:a.n
         if i != j
            entries[i, j] = zero(R)
         else
            entries[i, j] = rb
         end
      end
   end
   z = MatRingElem{T}(R, entries)
   return z
end

function (a::MatRing{T})(b::MatRingElem{T}) where {T <: NCRingElement}
   parent(b) != a && error("Unable to coerce matrix")
   return b
end

function (a::MatRing{T})(b::MatrixElem{S}) where {S <: NCRingElement, T <: NCRingElement}
   R = base_ring(a)
   _check_dim(nrows(a), ncols(a), b)
   entries = Matrix{T}(undef, nrows(a), ncols(a))
   for i = 1:nrows(a)
      for j = 1:ncols(a)
         entries[i, j] = R(b[i, j])
      end
   end
   z = MatRingElem{T}(R, entries)
   return z
end

function (a::MatRing{T})(b::Matrix{S}) where {S <: NCRingElement, T <: NCRingElement}
   R = base_ring(a)
   _check_dim(a.n, a.n, b)
   entries = Matrix{T}(undef, a.n, a.n)
   for i = 1:a.n
      for j = 1:a.n
         entries[i, j] = R(b[i, j])
      end
   end
   z = MatRingElem{T}(R, entries)
   return z
end

function (a::MatRing{T})(b::Vector{S}) where {S <: NCRingElement, T <: NCRingElement}
   _check_dim(a.n, a.n, b)
   b = Matrix{S}(transpose(reshape(b, a.n, a.n)))
   z = a(b)
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
