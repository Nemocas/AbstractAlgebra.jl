###############################################################################
#
#   MatrixAlgebra.jl : nxn matrices over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring_type(::Type{MatAlgebra{T}}) where T <: NCRingElement = parent_type(T)

function base_ring(a::MatAlgebra{T}) where {T <: NCRingElement}
   a.base_ring::parent_type(T)
end

function check_parent(a::MatAlgElem{T}, b::MatAlgElem{T}, throw::Bool = true) where T <: NCRingElement
  fl = (base_ring(a) != base_ring(b) || degree(a) != degree(b))
  fl && throw && error("Incompatible matrix spaces in matrix operation")
  return !fl
end

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

number_of_rows(a::MatAlgebra) = a.n
number_of_columns(a::MatAlgebra) = number_of_rows(a)

@doc raw"""
    degree(a::MatAlgebra)

Return the degree $n$ of the given matrix algebra.
"""
degree(a::MatAlgebra) = nrows(a)

@doc raw"""
    degree(a::MatAlgElem{T}) where T <: RingElement

Return the degree $n$ of the given matrix algebra.
"""
degree(a::MatAlgElem{T}) where T <: NCRingElement = degree(parent(a))

zero(a::MatAlgebra) = a()

one(a::MatAlgebra) = a(1)

is_unit(a::MatAlgElem{T}) where T <: RingElement = is_unit(det(a))

is_unit(a::MatAlgElem{T}) where T <: FieldElement = rank(a) == degree(a)

# proof over a commutative ring: use adj(A)*A = det(A)*I = A*adj(A)
is_zero_divisor(a::MatAlgElem{T}) where T <: RingElement = is_zero_divisor(det(a))

is_zero_divisor(a::MatAlgElem{T}) where T <: FieldElement = rank(a) != degree(a)

function is_zero_divisor_with_annihilator(a::MatAlgElem{T}) where T <: RingElement
   f, b = is_zero_divisor_with_annihilator(det(a))
   throw(NotImplementedError(:adj, a)) #return f, b*adj(A)
end


function characteristic(a::MatAlgebra)
   return characteristic(base_ring(a))
end

###############################################################################
#
#   Similar and zero
#
###############################################################################

@doc raw"""
    similar(x::Generic.MatrixElem, R::NCRing=base_ring(x))
    similar(x::Generic.MatrixElem, R::NCRing, r::Int, c::Int)
    similar(x::Generic.MatrixElem, r::Int, c::Int)
    similar(x::MatAlgElem, R::NCRing, n::Int)
    similar(x::MatAlgElem, n::Int)

Create an uninitialized matrix over the given ring and dimensions,
with defaults based upon the given source matrix `x`.
"""
similar(x::MatAlgElem, R::NCRing, n::Int) = _similar(x, R, n, n)

similar(x::MatAlgElem, R::NCRing=base_ring(x)) = similar(x, R, degree(x))

similar(x::MatAlgElem, n::Int) = similar(x, base_ring(x), n)

function similar(x::MatAlgElem{T}, R::NCRing, m::Int, n::Int) where T <: NCRingElement
   m != n && error("Dimensions don't match in similar")
   return similar(x, R, n)
end

similar(x::MatAlgElem, m::Int, n::Int) = similar(x, base_ring(x), m, n)

@doc raw"""
    zero(x::MatrixElem, R::NCRing=base_ring(x))
    zero(x::MatrixElem, R::NCRing, r::Int, c::Int)
    zero(x::MatrixElem, r::Int, c::Int)
    zero(x::MatAlgElem, R::NCRing, n::Int)
    zero(x::MatAlgElem, n::Int)

Create a zero matrix over the given ring and dimensions,
with defaults based upon the given source matrix `x`.
"""
zero(x::MatAlgElem, R::NCRing, n::Int) = zero!(similar(x, R, n))
zero(x::MatAlgElem, n::Int) = zero!(similar(x, n))

################################################################################
#
#  Copy and deepcopy
#
################################################################################

function copy(d::MatAlgElem{T}) where T <: NCRingElement
   z = similar(d)
   for i = 1:nrows(d)
      for j = 1:ncols(d)
         z[i, j] = d[i, j]
      end
   end
   return z
end

function deepcopy_internal(d::MatAlgElem{T}, dict::IdDict) where T <: NCRingElement
   z = similar(d)
   for i = 1:nrows(d)
      for j = 1:ncols(d)
         z[i, j] = deepcopy_internal(d[i, j], dict)
      end
   end
   return z
end

################################################################################
#
#   is_square
#
################################################################################

is_square(a::MatAlgElem) = true

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, ::MIME"text/plain", a::MatAlgebra)
  print(io, "Matrix algebra of")
  print(io, " degree ", a.n)
  println(io)
  io = pretty(io)
  print(io, Indent(), "over ")
  print(io, Lowercase(), base_ring(a))
end

function show(io::IO, a::MatAlgebra)
   if get(io, :supercompact, false)
      print(io, "Matrix algebra")
   else
      io = pretty(io)
      print(io, "Matrix algebra of ")
      print(io, "degree ", a.n, " over ")
      print(IOContext(io, :supercompact => true), Lowercase(), base_ring(a))
   end
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function *(x::MatAlgElem{T}, y::MatAlgElem{T}) where {T <: NCRingElement}
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

function ==(x::MatAlgElem, y::Union{Integer, Rational, AbstractFloat})
   n = degree(x)
   for i = 1:n
      if x[i, i] != y
         return false
      end
   end
   for i = 1:n
      for j = 1:n
         if i != j && !is_zero_entry(x, i, j)
            return false
         end
      end
   end
   return true
end

==(x::Union{Integer, Rational, AbstractFloat}, y::MatAlgElem) = y == x

function ==(x::MatAlgElem{T}, y::T) where T <: NCRingElem
   n = degree(x)
   for i = 1:n
      if x[i, i] != y
         return false
      end
   end
   for i = 1:n
      for j = 1:n
         if i != j && !is_zero_entry(x, i, j)
            return false
         end
      end
   end
   return true
end

==(x::T, y::MatAlgElem{T}) where T <: NCRingElem = y == x

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact_left(f::MatAlgElem{T},
                       g::MatAlgElem{T}; check::Bool=true) where T <: RingElement
   ginv, d = pseudo_inv(g)
   return divexact(ginv*f, d; check=check)
end

function divexact_right(f::MatAlgElem{T},
                       g::MatAlgElem{T}; check::Bool=true) where T <: RingElement
   ginv, d = pseudo_inv(g)
   return divexact(f*ginv, d; check=check)
end

function divexact_left(f::MatAlgElem{T},
                       g::MatAlgElem{T}; check::Bool=true) where T <: FieldElement
   return inv(g)*f
end

function divexact_right(f::MatAlgElem{T},
                       g::MatAlgElem{T}; check::Bool=true) where T <: FieldElement
   return f*inv(g)
end

###############################################################################
#
#   Gram
#
###############################################################################

@doc raw"""
    gram(x::MatAlgElem)

Return the Gram matrix of $x$, i.e. if $x$ is an $r\times c$ matrix return
the $r\times r$ matrix whose entries $i, j$ are the dot products of the
$i$-th and $j$-th rows, respectively.
"""
function gram(x::MatAlgElem)
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


RandomExtensions.maketype(S::MatAlgebra, _) = elem_type(S)

function RandomExtensions.make(S::MatAlgebra, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1]) # forward to default Make constructor
   else
      Make(S, make(R, vs...))
   end
end

Random.Sampler(::Type{RNG}, S::MatAlgebra, n::Random.Repetition
               ) where {RNG<:AbstractRNG} =
   Random.Sampler(RNG, make(S), n)

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{<:MatAlgElem,
                                         <:MatAlgebra}})
   S, v = sp[][1:end]
   M = S()
   n = degree(M)
   R = base_ring(S)
   for i = 1:n
      for j = 1:n
         M[i, j] = rand(rng, v)
      end
   end
   return M
end

rand(rng::AbstractRNG, S::MatAlgebra, v...) = rand(rng, make(S, v...))

rand(S::MatAlgebra, v...) = rand(Random.GLOBAL_RNG, S, v...)

# resolve ambiguities
rand(rng::AbstractRNG, S::MatAlgebra, dims::Integer...) =
   rand(rng, make(S), dims...)

rand(S::MatAlgebra, dims::Integer...) = rand(Random.GLOBAL_RNG, S, dims...)

function randmat_triu(rng::AbstractRNG, S::MatAlgebra, v...)
   M = S()
   n = degree(M)
   R = base_ring(S)
   for i = 1:n
      for j = 1:i - 1
         M[i, j] = R()
      end
      for j = i:n
         M[i, j] = rand(rng, R, v...)
      end
      while is_zero_entry(M, i, i)
         M[i, i] = rand(rng, R, v...)
      end
   end
   return M
end

randmat_triu(S::MatAlgebra, v...) = randmat_triu(Random.GLOBAL_RNG, S, v...)

function randmat_with_rank(rng::AbstractRNG, S::MatAlgebra{T}, rank::Int, v...) where {T <: RingElement}
   M = S()
   n = degree(M)
   R = base_ring(S)
   for i = 1:rank
      for j = 1:i - 1
         M[i, j] = R()
      end
      M[i, i] = rand(rng, R, v...)
      while is_zero_entry(M, i, i)
         M[i, i] = rand(rng, R, v...)
      end
      for j = i + 1:n
         M[i, j] = rand(rng, R, v...)
      end
   end
   for i = rank + 1:n
      for j = 1:n
         M[i, j] = R()
      end
   end
   if n > 1
      for i = 1:4*n
         r1 = rand(rng, 1:n)
         r2 = rand(rng, 1:n - 1)
         r2 = r2 >= r1 ? r2 + 1 : r2
         d = rand(rng, -5:5)
         for j = 1:n
            M[r1, j] = M[r1, j] + d*M[r2, j]
         end
      end
   end
   return M
end

randmat_with_rank(S::MatAlgebra{T}, rank::Int, v...) where {T <: RingElement} =
   randmat_with_rank(Random.GLOBAL_RNG, S, rank, v...)

###############################################################################
#
#   Identity matrix
#
###############################################################################

function identity_matrix(M::MatAlgElem{T}, n::Int) where T <: NCRingElement
   R = base_ring(M)
   arr = Matrix{T}(undef, n, n)
   for i in 1:n
      for j in 1:n
         arr[i, j] = i == j ? one(R) : zero(R)
      end
   end
   z = Generic.MatAlgElem{T}(R, arr)
   return z
end

@doc raw"""
    identity_matrix(M::MatAlgElem{T}) where T <: RingElement

Return the identity matrix over the same base ring as $M$ and with the
same dimensions.
"""
function identity_matrix(M::MatAlgElem{T}) where T <: NCRingElement
   return identity_matrix(M, nrows(M))
end

###############################################################################
#
#   MatrixAlgebra constructor
#
###############################################################################

@doc raw"""
    MatrixAlgebra(R::Ring, n::Int)

Return parent object corresponding to the ring of $n\times n$ matrices over
the ring $R$.
"""
function MatrixAlgebra(R::NCRing, n::Int)
   Generic.MatrixAlgebra(R, n)
end
