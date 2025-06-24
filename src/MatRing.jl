###############################################################################
#
#   MatRing.jl : nxn matrices over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring_type(::Type{<:MatRing{T}}) where T <: NCRingElement = parent_type(T)

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::MatRingElem, h::UInt)
   b = 0x6413942b83a26c65%UInt
   for i in 1:nrows(a)
      for j in 1:ncols(a)
         b = xor(b, xor(hash(a[i, j], h), h))
         b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
      end
   end
   return b
end

number_of_rows(a::MatRing) = a.n
number_of_columns(a::MatRing) = number_of_rows(a)

@doc raw"""
    degree(a::MatRing)

Return the degree $n$ of the given matrix algebra.
"""
degree(a::MatRing) = nrows(a)

@doc raw"""
    degree(a::MatRingElem{T}) where T <: RingElement

Return the degree $n$ of the given matrix algebra.
"""
degree(a::MatRingElem{T}) where T <: NCRingElement = degree(parent(a))

zero(a::MatRing) = a()

one(a::MatRing) = a(1)

is_unit(a::MatRingElem{T}) where T <: RingElement = is_unit(det(a))

is_unit(a::MatRingElem{T}) where T <: FieldElement = rank(a) == degree(a)

ConformanceTests._implements(::Type{MatRingElem{T}}, ::typeof(is_unit)) where {T <: RingElement} = _implements(T, is_unit)

# proof over a commutative ring: use adj(A)*A = det(A)*I = A*adj(A)
is_zero_divisor(a::MatRingElem{T}) where T <: RingElement = is_zero_divisor(det(a))

is_zero_divisor(a::MatRingElem{T}) where T <: FieldElement = rank(a) != degree(a)

ConformanceTests._implements(::Type{MatRingElem{T}}, ::typeof(is_zero_divisor)) where {T <: RingElement} = _implements(T, is_zero_divisor)

function is_zero_divisor_with_annihilator(a::MatRingElem{T}) where T <: RingElement
   f, b = is_zero_divisor_with_annihilator(det(a))
   throw(NotImplementedError(:adj, a)) #return f, b*adj(A)
end


function characteristic(a::MatRing)
   iszero(nrows(a)) && return 1
   return characteristic(base_ring(a))
end

is_finite(R::MatRing) = iszero(nrows(a)) || is_finite(base_ring(R))

###############################################################################
#
#   Similar and zero
#
###############################################################################

@doc raw"""
    similar(x::MatRingElem, R::NCRing, n::Int)
    similar(x::MatRingElem, R::NCRing)
    similar(x::MatRingElem, n::Int)
    similar(x::MatRingElem)

Create an uninitialized matrix ring element over the given ring and dimension,
with defaults based upon the given source matrix ring element `x`.
"""
function similar(x::MatRingElem, R::NCRing=base_ring(x), n::Int=degree(x))
   TT = elem_type(R)
   M = Matrix{TT}(undef, (n, n))
   return Generic.MatRingElem{TT}(R, M)
end

similar(x::MatRingElem, n::Int) = similar(x, base_ring(x), n)

# TODO: deprecate these:
function similar(x::MatRingElem{T}, R::NCRing, m::Int, n::Int) where T <: NCRingElement
   m != n && error("Dimensions don't match in similar")
   return similar(x, R, n)
end

similar(x::MatRingElem, m::Int, n::Int) = similar(x, base_ring(x), m, n)

@doc raw"""
    zero(x::MatRingElem, R::NCRing, n::Int)
    zero(x::MatRingElem, R::NCRing)
    zero(x::MatRingElem, n::Int)
    zero(x::MatRingElem)

Create a zero matrix ring element over the given ring and dimension,
with defaults based upon the given source matrix ring element `x`.
"""
zero(x::MatRingElem, R::NCRing=base_ring(x), n::Int=degree(x)) = zero!(similar(x, R, n))
zero(x::MatRingElem, n::Int) = zero!(similar(x, n))

# TODO: deprecate these
zero(x::MatRingElem, R::NCRing, r::Int, c::Int) = zero!(similar(x, R, r, c))
zero(x::MatRingElem, r::Int, c::Int) = zero!(similar(x, r, c))

################################################################################
#
#  Copy and deepcopy
#
################################################################################

function copy(d::MatRingElem{T}) where T <: NCRingElement
   z = similar(d)
   for i = 1:nrows(d)
      for j = 1:ncols(d)
         z[i, j] = d[i, j]
      end
   end
   return z
end

function deepcopy_internal(d::MatRingElem{T}, dict::IdDict) where T <: NCRingElement
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

is_square(a::MatRingElem) = true

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, mime::MIME"text/plain", a::MatRing)
  print(io, "Matrix ring of")
  print(io, " degree ", a.n)
  println(io)
  io = pretty(io)
  print(io, Indent(), "over ")
  print(io, Lowercase(), base_ring(a))
end

function show(io::IO, a::MatRing)
   if is_terse(io)
      print(io, "Matrix ring")
   else
      io = pretty(io)
      print(io, "Matrix ring of ")
      print(io, "degree ", a.n, " over ")
      print(terse(io), Lowercase(), base_ring(a))
   end
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function *(x::MatRingElem{T}, y::MatRingElem{T}) where {T <: NCRingElement}
   degree(x) != degree(y) && error("Incompatible matrix degrees")
   A = similar(x)
   C = base_ring(x)()
   for i = 1:nrows(x)
      for j = 1:ncols(y)
         A[i, j] = base_ring(x)()
         for k = 1:ncols(x)
            C = mul!(C, x[i, k], y[k, j])
            A[i, j] = add!(A[i, j], C)
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

function ==(x::MatRingElem, y::Union{Integer, Rational, AbstractFloat})
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

==(x::Union{Integer, Rational, AbstractFloat}, y::MatRingElem) = y == x

function ==(x::MatRingElem{T}, y::T) where T <: NCRingElem
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

==(x::T, y::MatRingElem{T}) where T <: NCRingElem = y == x

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact_left(f::MatRingElem{T},
                       g::MatRingElem{T}; check::Bool=true) where T <: RingElement
   ginv, d = pseudo_inv(g)
   return divexact(ginv*f, d; check=check)
end

function divexact_right(f::MatRingElem{T},
                       g::MatRingElem{T}; check::Bool=true) where T <: RingElement
   ginv, d = pseudo_inv(g)
   return divexact(f*ginv, d; check=check)
end

function divexact_left(f::MatRingElem{T},
                       g::MatRingElem{T}; check::Bool=true) where T <: FieldElement
   return inv(g)*f
end

function divexact_right(f::MatRingElem{T},
                       g::MatRingElem{T}; check::Bool=true) where T <: FieldElement
   return f*inv(g)
end

###############################################################################
#
#   Gram
#
###############################################################################

@doc raw"""
    gram(x::MatRingElem)

Return the Gram matrix of $x$, i.e. if $x$ is an $r\times c$ matrix return
the $r\times r$ matrix whose entries $i, j$ are the dot products of the
$i$-th and $j$-th rows, respectively.
"""
function gram(x::MatRingElem)
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


RandomExtensions.maketype(S::MatRing, _) = elem_type(S)

function RandomExtensions.make(S::MatRing, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1]) # forward to default Make constructor
   else
      Make(S, make(R, vs...))
   end
end

Random.Sampler(::Type{RNG}, S::MatRing, n::Random.Repetition
               ) where {RNG<:AbstractRNG} =
   Random.Sampler(RNG, make(S), n)

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{<:MatRingElem,
                                         <:MatRing}})
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

rand(rng::AbstractRNG, S::MatRing, v...) = rand(rng, make(S, v...))

rand(S::MatRing, v...) = rand(Random.default_rng(), S, v...)

# resolve ambiguities
rand(rng::AbstractRNG, S::MatRing, dims::Integer...) =
   rand(rng, make(S), dims...)

rand(S::MatRing, dims::Integer...) = rand(Random.default_rng(), S, dims...)

function randmat_triu(rng::AbstractRNG, S::MatRing, v...)
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

randmat_triu(S::MatRing, v...) = randmat_triu(Random.default_rng(), S, v...)

function randmat_with_rank(rng::AbstractRNG, S::MatRing{T}, rank::Int, v...) where {T <: RingElement}
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

randmat_with_rank(S::MatRing{T}, rank::Int, v...) where {T <: RingElement} =
   randmat_with_rank(Random.default_rng(), S, rank, v...)

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(S::MatRing)
  R = base_ring(S)
  return S(elem_type(R)[ConformanceTests.generate_element(R) for i in 1:nrows(S), j in 1:ncols(S)])
end

###############################################################################
#
#   Identity matrix
#
###############################################################################

function identity_matrix(M::MatRingElem{T}, n::Int) where T <: NCRingElement
   R = base_ring(M)
   arr = Matrix{T}(undef, n, n)
   for i in 1:n
      for j in 1:n
         arr[i, j] = i == j ? one(R) : zero(R)
      end
   end
   z = Generic.MatRingElem{T}(R, arr)
   return z
end

@doc raw"""
    identity_matrix(M::MatRingElem{T}) where T <: RingElement

Return the identity matrix over the same base ring as $M$ and with the
same dimensions.
"""
function identity_matrix(M::MatRingElem{T}) where T <: NCRingElement
   return identity_matrix(M, nrows(M))
end

###############################################################################
#
#   MatRing constructor
#
###############################################################################

@doc raw"""
    matrix_ring(R::Ring, n::Int)

Return parent object corresponding to the ring of $n\times n$ matrices over
the ring $R$.
"""
function matrix_ring(R::NCRing, n::Int)
   Generic.matrix_ring(R, n)
end
