###############################################################################
#
#   MatrixAlgebra.jl : Generic nxn matrices over rings
#
###############################################################################

export MatrixAlgebra, dimension

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{MatAlgElem{T}}) where T <: RingElement = MatAlgebra{T}

elem_type(::Type{MatAlgebra{T}}) where {T <: RingElement} = MatAlgElem{T}

Markdown.doc"""
    base_ring{T <: RingElement}(S::AbstractAlgebra.MatAlgebra{T})
> Return the base ring $R$ of the given matrix algebra.
"""
function base_ring(a::AbstractAlgebra.MatAlgebra{T}) where {T <: RingElement}
   a.base_ring::parent_type(T)
end

function check_parent(a::AbstractAlgebra.MatAlgElem{T}, b::AbstractAlgebra.MatAlgElem{T}) where T <: RingElement
  (base_ring(a) != base_ring(b) || dimension(a) != dimension(b)) &&
                error("Incompatible matrix spaces in matrix operation")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::MatAlgElem, h::UInt)
   b = 0x6413942b83a26c65%UInt
   for i in 1:rows(a)
      for j in 1:cols(a)
         b = xor(b, xor(hash(a[i, j], h), h))
         b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
      end
   end
   return b
end

dimension(a::MatAlgElem) = size(a.entries, 1)

Markdown.doc"""
    zero(a::AbstractAlgebra.MatAlgebra)
> Construct the zero matrix in the given matrix space.
"""
zero(a::AbstractAlgebra.MatAlgebra) = a()

Markdown.doc"""
    one(a::AbstractAlgebra.MatAlgebra)
> Construct the matrix in the given matrix space with ones down the diagonal
> and zeroes elsewhere.
"""
one(a::AbstractAlgebra.MatAlgebra) = a(1)

###############################################################################
#
#   Similar and eye
#
###############################################################################

function similar(x::MatAlgElem{T}) where T <: RingElement
   R = base_ring(x)
   M = similar(x.entries)
   for i in 1:size(M, 1)
      for j in 1:size(M, 2)
         M[i, j] = zero(R)
      end
   end
   z = MatAlgElem{T}(M)
   z.base_ring = R
   return z
end

function similar(x::MatAlgElem{T}, n::Int) where T <: RingElement
   R = base_ring(x)
   M = similar(x.entries, n, n)
   for i in 1:size(M, 1)
      for j in 1:size(M, 2)
         M[i, j] = zero(R)
      end
   end
   z = MatAlgElem{T}(M)
   z.base_ring = R
   return z
end

################################################################################
#
#   Size
#
################################################################################

issquare(a::MatAlgElem) = true

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, a::AbstractAlgebra.MatAlgebra)
   print(io, "Matrix Algebra of dimension ")
   print(io, a.n, " over ")
   print(io, base_ring(a))
end

show_minus_one(::Type{AbstractAlgebra.MatAlgElem{T}}) where {T <: RingElement} = false

###############################################################################
#
#   Binary operations
#
###############################################################################

function *(x::AbstractAlgebra.MatAlgElem{T}, y::AbstractAlgebra.MatAlgElem{T}) where {T <: RingElement}
   dimension(x) != dimension(y) && error("Incompatible matrix dimensions")
   A = similar(x)
   C = base_ring(x)()
   for i = 1:rows(x)
      for j = 1:cols(y)
         A[i, j] = base_ring(x)()
         for k = 1:cols(x)
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
   n = dimension(x)
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
   n = dimension(x)
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

Markdown.doc"""
    MatrixAlgebra(R::AbstractAlgebra.Ring, n::Int, cached::Bool = true)
> Return parent object corresponding to the ring of $n\times n$ matrices over
> the ring $R$. If `cached == true` (the default), the returned parent object
> is cached so that it can returned by future calls to the constructor with the
> same dimensions and base ring.
"""
function MatrixAlgebra(R::AbstractAlgebra.Ring, n::Int, cached::Bool = true)
   T = elem_type(R)
   return MatAlgebra{T}(R, n, cached)
end

