###############################################################################
#
#   MatrixAlgebra.jl : Generic nxn matrices over rings
#
###############################################################################

export MatrixAlgebra, divexact_left, divexact_right

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{MatAlgElem{T}}) where T <: RingElement = MatAlgebra{T}

elem_type(::Type{MatAlgebra{T}}) where {T <: RingElement} = MatAlgElem{T}

@doc Markdown.doc"""
    base_ring(a::AbstractAlgebra.MatAlgebra{T}) where {T <: RingElement}
> Return the base ring $R$ of the given matrix algebra.
"""
function base_ring(a::AbstractAlgebra.MatAlgebra{T}) where {T <: RingElement}
   a.base_ring::parent_type(T)
end

@doc Markdown.doc"""
    parent(a::AbstractAlgebra.MatAlgElem{T}, cached::Bool = true) where T <: RingElement
> Return the parent object of the given matrix.
"""
parent(a::AbstractAlgebra.MatAlgElem{T}, cached::Bool = true) where T <: RingElement =
    MatAlgebra{T}(a.base_ring, size(a.entries)[1], cached)

function check_parent(a::AbstractAlgebra.MatAlgElem{T}, b::AbstractAlgebra.MatAlgElem{T}, throw::Bool = true) where T <: RingElement
  fl = (base_ring(a) != base_ring(b) || degree(a) != degree(b))
  fl && throw && error("Incompatible matrix spaces in matrix operation")
  return !fl
end

isexact_type(::Type{MatAlgElem{T}}) where T <: RingElement = isexact_type(T)

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::MatAlgElem, h::UInt)
   b = 0x6413942b83a26c65%UInt
   for i in 1:nrows(a)
      for j in 1:ncols(a)
         b = xor(b, xor(hash(a[i, j], h), h))
         b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
      end
   end
   return b
end

nrows(a::AbstractAlgebra.MatAlgebra) = a.n
ncols(a::AbstractAlgebra.MatAlgebra) = nrows(a)

@doc Markdown.doc"""
    degree(a::Generic.MatAlgElem)
> Return the degree $n$ of the $n\times n$ matrix $a$..
"""
degree(a::MatAlgElem) = nrows(a)

@doc Markdown.doc"""
    degree(a::AbstractAlgebra.MatAlgebra)
> Return the degree $n$ of the given matrix algebra.
"""
degree(a::AbstractAlgebra.MatAlgebra) = nrows(a)

@doc Markdown.doc"""
    zero(a::AbstractAlgebra.MatAlgebra)
> Construct the zero matrix in the given matrix algebra.
"""
zero(a::AbstractAlgebra.MatAlgebra) = a()

@doc Markdown.doc"""
    one(a::AbstractAlgebra.MatAlgebra)
> Construct the matrix in the given matrix algebra with ones down the diagonal
> and zeroes elsewhere.
"""
one(a::AbstractAlgebra.MatAlgebra) = a(1)

isunit(a::AbstractAlgebra.MatAlgElem{T}) where T <: RingElement = isunit(det(a))

isunit(a::AbstractAlgebra.MatAlgElem{T}) where T <: FieldElement = rank(a) == degree(a)

###############################################################################
#
#   Similar and eye
#
###############################################################################


@doc Markdown.doc"""
    similar(x::Generic.MatrixElem, R::Ring=base_ring(x))
    similar(x::Generic.MatrixElem, R::Ring, r::Int, c::Int)
    similar(x::Generic.MatrixElem, r::Int, c::Int)
    similar(x::MatAlgElem, R::Ring, n::Int)
    similar(x::MatAlgElem, n::Int)

> Create a matrix over the given ring and dimensions,
> with defaults based upon the given source matrix `x`.
"""
similar(x::AbstractAlgebra.MatAlgElem, R::Ring, n::Int) = _similar(x, R, n, n)

similar(x::AbstractAlgebra.MatAlgElem, R::Ring=base_ring(x)) = similar(x, R, degree(x))

similar(x::AbstractAlgebra.MatAlgElem, n::Int) = similar(x, base_ring(x), n)

function similar(x::AbstractAlgebra.MatAlgElem{T}, R::Ring, m::Int, n::Int) where T <: RingElement
   m != n && error("Dimensions don't match in similar")
   return similar(x, R, n)
end

similar(x::AbstractAlgebra.MatAlgElem, m::Int, n::Int) = similar(x, base_ring(x), m, n)

################################################################################
#
#   Size
#
################################################################################

size(x::MatAlgElem) = tuple(nrows(x), ncols(x))

size(t::MatAlgElem, d) = d <= 2 ? size(t)[d] : 1

issquare(a::MatAlgElem) = true

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, a::AbstractAlgebra.MatAlgebra)
   print(io, "Matrix Algebra of degree ")
   print(io, a.n, " over ")
   print(IOContext(io, :compact => true), base_ring(a))
end

show_minus_one(::Type{AbstractAlgebra.MatAlgElem{T}}) where T <: RingElement = false

needs_parentheses(a::AbstractAlgebra.MatAlgElem{T}) where T <: RingElement = true

displayed_with_minus_in_front(a::AbstractAlgebra.MatAlgElem{T}) where T <: RingElement = false

###############################################################################
#
#   Binary operations
#
###############################################################################

function *(x::AbstractAlgebra.MatAlgElem{T}, y::AbstractAlgebra.MatAlgElem{T}) where {T <: RingElement}
   degree(x) != degree(y) && error("Incompatible matrix degrees")
   A = similar(x)
   C = base_ring(x)()
   for i = 1:nrows(x)
      for j = 1:ncols(y)
         A[i, j] = base_ring(x)()
         for k = 1:ncols(x)
            C = mul!(C, x[i, k], y[k, j])
            A[i, j] = addeq!(A[i, j], C)
         end
      end
   end
   return A
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

function ==(x::AbstractAlgebra.MatAlgElem, y::Union{Integer, Rational, AbstractFloat})
   n = degree(x)
   for i = 1:n
      if x[i, i] != y
         return false
      end
   end
   for i = 1:n
      for j = 1:n
         if i != j && !iszero(x[i, j])
            return false
         end
      end
   end
   return true
end

==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.MatAlgElem) = y == x

function ==(x::AbstractAlgebra.MatAlgElem{T}, y::T) where {T <: RingElem}
   n = degree(x)
   for i = 1:n
      if x[i, i] != y
         return false
      end
   end
   for i = 1:n
      for j = 1:n
         if i != j && !iszero(x[i, j])
            return false
         end
      end
   end
   return true
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact_left(f::AbstractAlgebra.MatAlgElem{T},
                       g::AbstractAlgebra.MatAlgElem{T}) where T <: RingElement
   ginv, d = pseudo_inv(g)
   return divexact(ginv*f, d)
end

function divexact_right(f::AbstractAlgebra.MatAlgElem{T},
                       g::AbstractAlgebra.MatAlgElem{T}) where T <: RingElement
   ginv, d = pseudo_inv(g)
   return divexact(f*ginv, d)
end

function divexact_left(f::AbstractAlgebra.MatAlgElem{T},
                       g::AbstractAlgebra.MatAlgElem{T}) where T <: FieldElement
   return inv(g)*f
end

function divexact_right(f::AbstractAlgebra.MatAlgElem{T},
                       g::AbstractAlgebra.MatAlgElem{T}) where T <: FieldElement
   return f*inv(g)
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact_left(x::AbstractAlgebra.MatAlgElem{T}, y::T) where {T <: RingElem}
   return divexact(x, y)
end

function divexact_right(x::AbstractAlgebra.MatAlgElem{T}, y::T) where {T <: RingElem}
   return divexact(x, y)
end

function divexact_left(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat})
   return divexact(x, y)
end

function divexact_right(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat})
   return divexact(x, y)
end

###############################################################################
#
#   Transpose
#
###############################################################################

@doc Markdown.doc"""
    transpose(x::MatAlgElem{T}) where T <: RingElement
> Return the transpose of the given matrix.
"""
function transpose(x::MatAlgElem{T}) where T <: RingElement
   arr = permutedims(x.entries, [2, 1])
   z = MatAlgElem{T}(arr)
   z.base_ring = base_ring(x)
   return z
end

@doc Markdown.doc"""
    gram(x::AbstractAlgebra.MatAlgElem)
> Return the Gram matrix of $x$, i.e. if $x$ is an $r\times c$ matrix return
> the $r\times r$ matrix whose entries $i, j$ are the dot products of the
> $i$-th and $j$-th rows, respectively.
"""
function gram(x::AbstractAlgebra.MatAlgElem)
   n = degree(x)
   z = similar(x)
   for i = 1:n
      for j = 1:n
         z[i, j] = zero(base_ring(x))
         for k = 1:n
            z[i, j] += x[i, k] * x[j, k]
         end
      end
   end
   return z
end

###############################################################################
#
#   Random generation
#
###############################################################################

function rand(S::AbstractAlgebra.MatAlgebra, v...)
   M = S()
   n = degree(M)
   R = base_ring(S)
   for i = 1:n
      for j = 1:n
         M[i, j] = rand(R, v...)
      end
   end
   return M
end

function randmat_triu(S::AbstractAlgebra.MatAlgebra, v...)
   M = S()
   n = degree(M)
   R = base_ring(S)
   for i = 1:n
      for j = 1:i - 1
         M[i, j] = R()
      end
      for j = i:n
         M[i, j] = rand(R, v...)
      end
      while iszero(M[i, i])
         M[i, i] = rand(R, v...)
      end
   end
   return M
end

function randmat_with_rank(S::Generic.MatAlgebra{T}, rank::Int, v...) where {T <: AbstractAlgebra.RingElement}
   M = S()
   n = degree(M)
   R = base_ring(S)
   for i = 1:rank
      for j = 1:i - 1
         M[i, j] = R()
      end
      M[i, i] = rand(R, v...)
      while iszero(M[i, i])
         M[i, i] = rand(R, v...)
      end
      for j = i + 1:n
         M[i, j] = rand(R, v...)
      end
   end
   for i = rank + 1:n
      for j = 1:n
         M[i, j] = R()
      end
   end
   if n > 1
      for i = 1:4*n
         r1 = rand(1:n)
         r2 = rand(1:n - 1)
         r2 = r2 >= r1 ? r2 + 1 : r2
         d = rand(-5:5)
         for j = 1:n
            M[r1, j] = M[r1, j] + d*M[r2, j]
         end
      end
   end
   return M
end

###############################################################################
#
#   Solve
#
###############################################################################

function solve_lu(M::MatAlgElem{T}, B::MatAlgElem{T}) where {T <: RingElement}
   check_parent(M, B)
   R = base_ring(M)
   MS = MatSpaceElem{T}(M.entries) # convert to ordinary matrix
   MS.base_ring = R
   BS = MatSpaceElem{T}(B.entries)
   BS.base_ring = R
   S = solve_lu(MS, BS)
   SA = MatAlgElem{T}(S.entries)
   SA.base_ring = R
   return SA
end

function solve_fflu(M::MatAlgElem{T}, B::MatAlgElem{T}) where {T <: RingElement}
   check_parent(M, B)
   R = base_ring(M)
   MS = MatSpaceElem{T}(M.entries) # convert to ordinary matrix
   MS.base_ring = R
   BS = MatSpaceElem{T}(B.entries)
   BS.base_ring = R
   S, d = solve_fflu(MS, BS)
   SA = MatAlgElem{T}(S.entries)
   SA.base_ring = R
   return SA, d
end

###############################################################################
#
#   Minimal polynomial
#
###############################################################################


@doc Markdown.doc"""
    minpoly(S::Ring, M::MatAlgElem{T}, charpoly_only::Bool = false) where {T <: RingElement}
> Return the minimal polynomial $p$ of the matrix $M$. The polynomial ring $S$
> of the resulting polynomial must be supplied and the matrix must be square.
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
         M.entries[i, j] = R()
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
         addeq!(A.entries[i, j], B.entries[i, j])
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
#   MatrixSpace constructor
#
###############################################################################

@doc Markdown.doc"""
    MatrixAlgebra(R::AbstractAlgebra.Ring, n::Int, cached::Bool = true)
> Return parent object corresponding to the ring of $n\times n$ matrices over
> the ring $R$. If `cached == true` (the default), the returned parent object
> is cached so that it can returned by future calls to the constructor with the
> same degree and base ring.
"""
function MatrixAlgebra(R::AbstractAlgebra.Ring, n::Int, cached::Bool = true)
   T = elem_type(R)
   return MatAlgebra{T}(R, n, cached)
end
