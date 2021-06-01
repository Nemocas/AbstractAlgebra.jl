###############################################################################
#
#   MatrixAlgebra.jl : Generic nxn matrices over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{MatAlgElem{T}}) where T <: RingElement = MatAlgebra{T}

elem_type(::Type{MatAlgebra{T}}) where {T <: RingElement} = MatAlgElem{T}

@doc Markdown.doc"""
    parent(a::MatAlgElem{T}, cached::Bool = true) where T <: RingElement

Return the parent object of the given matrix.
"""
parent(a::MatAlgElem{T}, cached::Bool = true) where T <: RingElement =
    MatAlgebra{T}(a.base_ring, size(a.entries)[1], cached)

isexact_type(::Type{MatAlgElem{T}}) where T <: RingElement = isexact_type(T)

###############################################################################
#
#   Transpose
#
###############################################################################

@doc Markdown.doc"""
    transpose(x::MatAlgElem{T}) where T <: RingElement

Return the transpose of the given matrix.
"""
function transpose(x::MatAlgElem{T}) where T <: RingElement
   arr = permutedims(x.entries, [2, 1])
   z = MatAlgElem{T}(arr)
   z.base_ring = base_ring(x)
   return z
end

###############################################################################
#
#   Solve
#
###############################################################################

function can_solve_with_solution_lu(M::MatAlgElem{T}, B::MatAlgElem{T}) where {T <: RingElement}
   check_parent(M, B)
   R = base_ring(M)
   MS = MatSpaceElem{T}(M.entries) # convert to ordinary matrix
   MS.base_ring = R
   BS = MatSpaceElem{T}(B.entries)
   BS.base_ring = R
   flag, S = can_solve_with_solution_lu(MS, BS)
   SA = MatAlgElem{T}(S.entries)
   SA.base_ring = R
   return flag, SA
end

function can_solve_with_solution_fflu(M::MatAlgElem{T}, B::MatAlgElem{T}) where {T <: RingElement}
   check_parent(M, B)
   R = base_ring(M)
   MS = MatSpaceElem{T}(M.entries) # convert to ordinary matrix
   MS.base_ring = R
   BS = MatSpaceElem{T}(B.entries)
   BS.base_ring = R
   flag, S, d = can_solve_with_solution_fflu(MS, BS)
   SA = MatAlgElem{T}(S.entries)
   SA.base_ring = R
   return flag, SA, d
end

###############################################################################
#
#   Minimal polynomial
#
###############################################################################

@doc Markdown.doc"""
    minpoly(S::Ring, M::MatAlgElem{T}, charpoly_only::Bool = false) where {T <: RingElement}

Return the minimal polynomial $p$ of the matrix $M$. The polynomial ring $S$
of the resulting polynomial must be supplied and the matrix must be square.
"""
function minpoly(S::Ring, M::MatAlgElem{T}, charpoly_only::Bool = false) where {T <: RingElement}
   MS = MatSpaceElem{T}(M.entries) # convert to ordinary matrix
   MS.base_ring = base_ring(M)
   return minpoly(S, MS, charpoly_only)
end

###############################################################################
#
#   Unsafe operators
#
###############################################################################

function zero!(M::MatAlgElem{T}) where T <: RingElement
   n = degree(M)
   R = base_ring(M)
   for i = 1:n
      for j = 1:n
         M[i, j] = zero(R)
      end
   end
   return M
end

function mul!(A::MatAlgElem{T}, B::MatAlgElem{T},
                                C::MatAlgElem{T}) where T <: RingElement
   return B*C
end

function add!(A::MatAlgElem{T}, B::MatAlgElem{T},
                                C::MatAlgElem{T}) where T <: RingElement
   n = degree(A)
   for i = 1:n
      for j = 1:n
         A.entries[i, j] = B.entries[i, j] + C.entries[i, j]
      end
   end
   return A
end

function addeq!(A::MatAlgElem{T}, B::MatAlgElem{T}) where T <: RingElement
   n = degree(A)
   for i = 1:n
      for j = 1:n
         A.entries[i, j] += B.entries[i, j]
      end
   end
   return A
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::MatAlgebra{T})() where {T <: RingElement}
   R = base_ring(a)
   entries = Array{T}(undef, a.n, a.n)
   for i = 1:a.n
      for j = 1:a.n
         entries[i, j] = zero(R)
      end
   end
   z = MatAlgElem{T}(entries)
   z.base_ring = R
   return z
end

function (a::MatAlgebra{T})(b::S) where {S <: RingElement, T <: RingElement}
   R = base_ring(a)
   entries = Array{T}(undef, a.n, a.n)
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
   z = MatAlgElem{T}(entries)
   z.base_ring = R
   return z
end

function (a::MatAlgebra{T})(b::MatAlgElem{T}) where {T <: RingElement}
   parent(b) != a && error("Unable to coerce matrix")
   return b
end

function (a::MatAlgebra{T})(b::Array{S, 2}) where {S <: RingElement, T <: RingElement}
   R = base_ring(a)
   _check_dim(a.n, a.n, b)
   entries = Array{T}(undef, a.n, a.n)
   for i = 1:a.n
      for j = 1:a.n
         entries[i, j] = R(b[i, j])
      end
   end
   z = MatAlgElem{T}(entries)
   z.base_ring = R
   return z
end

function (a::MatAlgebra{T})(b::Array{S, 1}) where {S <: RingElement, T <: RingElement}
   _check_dim(a.n, a.n, b)
   b = Array{S, 2}(transpose(reshape(b, a.n, a.n)))
   z = a(b)
   return z
end

###############################################################################
#
#   MatrixAlgebra constructor
#
###############################################################################

function MatrixAlgebra(R::AbstractAlgebra.Ring, n::Int; cached::Bool = true)
   T = elem_type(R)
   return MatAlgebra{T}(R, n, cached)
end





