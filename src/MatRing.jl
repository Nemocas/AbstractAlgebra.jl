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

base_ring_type(::Type{<:MatRingElem{T}}) where T <: NCRingElement = parent_type(T)

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::MatRingElem, h::UInt)
  b = 0x6413942b83a26c65 % UInt
  return xor(b, hash(matrix(a), h))
end

vector_space_dim(R::MatRing{T}) where {T <: Union{FieldElem, Rational{BigInt}}} = nrows(R)^2

@doc raw"""
    degree(R::MatRing)

Return the degree $n$ of the given matrix algebra.

The degree is the number of rows of the square matrices belonging to `R`.

# Examples

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S = matrix_ring(R, 3)
Matrix ring of degree 3
  over univariate polynomial ring in t over rationals

julia> degree(S)
3
```
"""
degree(R::MatRing) = nrows(R)

@doc raw"""
    degree(M::MatRingElem{T}) where T <: NCRingElement

Return the degree $n$ of the parent matrix algebra of `M`.

# Examples

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S = matrix_ring(R, 3)
Matrix ring of degree 3
  over univariate polynomial ring in t over rationals

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t + 1       t             1]
[  t^2       t             t]
[   -2   t + 2   t^2 + t + 1]

julia> degree(A)
3
```
"""
degree(M::MatRingElem{T}) where T <: NCRingElement = degree(parent(M))

zero(R::MatRing) = R()

one(R::MatRing) = R(1)

is_unit(M::MatRingElem{T}) where T <: RingElement = is_unit(det(M))

is_unit(M::MatRingElem{T}) where T <: FieldElement = rank(M) == degree(M)

# proof over a commutative ring: use adj(A)*A = det(A)*I = A*adj(A)
is_zero_divisor(M::MatRingElem{T}) where T <: RingElement = is_zero_divisor(det(M))

is_zero_divisor(M::MatRingElem{T}) where T <: FieldElement = rank(M) != degree(M)

function is_zero_divisor_with_annihilator(M::MatRingElem{T}) where T <: RingElement
   f, b = is_zero_divisor_with_annihilator(det(M))
   throw(NotImplementedError(:adj, M)) #return f, b*adj(A)
end

characteristic(R::MatRing) = iszero(nrows(R)) ? 1 : characteristic(base_ring(R))
is_known(::typeof(characteristic), R::MatRing) = iszero(nrows(R)) || is_known(characteristic, base_ring(R))

is_finite(R::MatRing) = iszero(nrows(R)) || is_finite(base_ring(R))
is_known(::typeof(is_finite), R::MatRing) = iszero(nrows(R)) || is_known(is_finite, base_ring(R))

###############################################################################
#
#   Similar, zero and iszero, one and isone
#
###############################################################################

@doc raw"""
    similar(M::MatRingElem, R::NCRing, n::Int)
    similar(M::MatRingElem, R::NCRing)
    similar(M::MatRingElem, n::Int)
    similar(M::MatRingElem)

Create an uninitialized matrix ring element over the given ring and dimension,
with defaults based upon the given source matrix ring element `M`.
"""
function similar(M::MatRingElem, R::NCRing=base_ring(M), n::Int=degree(M))
   @req (n >= 0)  "Matrix dimension must be non-negative"
   @req (n < 2^30)  "Matrix dimension is excessively large"
   return Generic.MatRingElem(R, n, fill(0,n^2)) # n^2 cannot overflow given check in line above
end

similar(M::MatRingElem, n::Int) = similar(M, base_ring(M), n)

# TODO: deprecate these:
function similar(M::MatRingElem{T}, R::NCRing, m::Int, n::Int) where T <: NCRingElement
   m != n && error("Dimensions don't match in similar")
   return similar(M, R, n)
end

similar(M::MatRingElem, m::Int, n::Int) = similar(M, base_ring(M), m, n)

@doc raw"""
    zero(M::MatRingElem, R::NCRing, n::Int)
    zero(M::MatRingElem, R::NCRing)
    zero(M::MatRingElem, n::Int)
    zero(M::MatRingElem)

Create a zero matrix ring element over the given ring and dimension,
with defaults based upon the given source matrix ring element `M`.
"""
zero(M::MatRingElem, R::NCRing=base_ring(M), n::Int=degree(M)) = zero!(similar(M, R, n))
zero(M::MatRingElem, n::Int) = zero!(similar(M, n))

# TODO: deprecate these
zero(M::MatRingElem, R::NCRing, r::Int, c::Int) = zero!(similar(M, R, r, c))
zero(M::MatRingElem, r::Int, c::Int) = zero!(similar(M, r, c))

iszero(M::MatRingElem{T}) where T <: NCRingElement = iszero(matrix(M))

one(M::MatRingElem{T}) where T <: NCRingElement = one(parent(M))

isone(M::MatRingElem{T}) where T <: NCRingElement = isone(matrix(M))

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(M::MatRingElem{T}) where T <: NCRingElement = canonical_unit(matrix(M))

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

is_square(::MatRingElem) = true   # FIXME: remove this once we untangled MatRingElem and MatrixElement etc.

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, mime::MIME"text/plain", R::MatRing)
  print(io, "Matrix ring of")
  print(io, " degree ", R.n)
  println(io)
  io = pretty(io)
  print(io, Indent(), "over ")
  print(io, Lowercase(), base_ring(R))
end

function show(io::IO, R::MatRing)
   if is_terse(io)
      print(io, "Matrix ring")
   else
      io = pretty(io)
      print(io, "Matrix ring of ")
      print(io, "degree ", R.n, " over ")
      print(terse(io), Lowercase(), base_ring(R))
   end
end

###############################################################################
#
#   Basic arithmetic -- delegate to MatElem
#
###############################################################################

*(x::T, y::T) where {T <: MatRingElem} = Generic.MatRingElem(matrix(x) * matrix(y))

+(x::T, y::T) where {T <: MatRingElem} = Generic.MatRingElem(matrix(x) + matrix(y))

-(x::T, y::T) where {T <: MatRingElem} = Generic.MatRingElem(matrix(x) - matrix(y))

^(a::MatRingElem{T}, b::Int) where T <: NCRingElement = Generic.MatRingElem(matrix(a)^b)

==(x::T, y::T) where {T <: MatRingElem} = matrix(x) == matrix(y)

isequal(x::T, y::T) where {T <: MatRingElem} = isequal(matrix(x), matrix(y))

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(x::JuliaRingElement, y::MatRingElem{T}) where T <: NCRingElement
  return Generic.MatRingElem(x * matrix(y))
end

function *(x::T, y::MatRingElem{T}) where {T <: NCRingElem}
  return Generic.MatRingElem(x * matrix(y))
end

function *(x::MatRingElem{T}, y::JuliaRingElement) where T <: NCRingElement
  return Generic.MatRingElem(matrix(x) * y)
end

function *(x::MatRingElem{T}, y::T) where {T <: NCRingElem}
  return Generic.MatRingElem(matrix(x) * y)
end

function +(x::JuliaRingElement, y::MatRingElem{T}) where T <: NCRingElement
  return Generic.MatRingElem(x + matrix(y))
end

function +(x::T, y::MatRingElem{T}) where {T <: NCRingElem}
  return Generic.MatRingElem(x + matrix(y))
end

function +(x::MatRingElem{T}, y::JuliaRingElement) where T <: NCRingElement
  return Generic.MatRingElem(matrix(x) + y)
end

function +(x::MatRingElem{T}, y::T) where {T <: NCRingElem}
  return Generic.MatRingElem(matrix(x) + y)
end

function -(x::JuliaRingElement, y::MatRingElem{T}) where T <: NCRingElement
  return Generic.MatRingElem(x - matrix(y))
end

function -(x::T, y::MatRingElem{T}) where {T <: NCRingElem}
  return Generic.MatRingElem(x - matrix(y))
end

function -(x::MatRingElem{T}, y::JuliaRingElement) where T <: NCRingElement
  return Generic.MatRingElem(matrix(x) - y)
end

function -(x::MatRingElem{T}, y::T) where {T <: NCRingElem}
  return Generic.MatRingElem(matrix(x) - y)
end

function *(x::MatRingElem{T}, y::Vector{T}) where T <: NCRingElement
  return matrix(x) * y
end

function *(x::Vector{T}, y::MatRingElem{T}) where T <: NCRingElement
  return x * matrix(y)
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

function ==(x::MatRingElem{T}, y::JuliaRingElement) where T <: NCRingElement
   return matrix(x) == y
end

==(x::JuliaRingElement, y::MatRingElem{T}) where T <: NCRingElement = y == x

function ==(x::MatRingElem{T}, y::T) where {T <: NCRingElem}
   return matrix(x) == y
end

==(x::T, y::MatRingElem{T}) where {T <: NCRingElem} = y == x

function ==(x::MatElem{T}, y::MatRingElem) where {T <: NCRingElement}
  error("Equality comparison of MatElem with MatRingElem unsupported. Call `x == matrix(y)` instead.")
end

function ==(x::MatRingElem, y::MatElem{T}) where {T <: NCRingElement}
  error("Equality comparison of MatRingElem with MatElem unsupported. Call `matrix(x) == y` instead.")
end

function isequal(x::MatElem{T}, y::MatRingElem{T}) where {T <: NCRingElement}
  error("Equality comparison of MatElem with MatRingElem unsupported. Call `isequal(x, matrix(y))` instead.")
end

function isequal(x::MatRingElem{T}, y::MatElem{T}) where {T <: NCRingElement}
  error("Equality comparison of MatRingElem with MatElem unsupported. Call `isequal(matrix(x), y)` instead.")
end

# to resolve ambiguities:
function ==(x::MatElem{T}, y::T) where {T <: MatRingElem}
   for i = 1:min(nrows(x), ncols(x))
      if x[i, i] != y
         return false
      end
   end
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if i != j && !is_zero_entry(x, i, j)
            return false
         end
      end
   end
   return true
end

==(x::T, y::MatElem{T}) where {T <: MatRingElem} = y == x

###############################################################################
#
#   Exact division
#
###############################################################################

# TO DO: using pseudo_inv is not ideal
#   consider case M = matrix(ZZmod4, 2,2, [2,1,0,1])
#   pseudo_inv(M) gives error, but copuld give matrix(ZZ4, 2,2, [1,-1,0,2]) with denom=2
# HINT: Consider using solve(f,g;side=:right)  or side=:left
# The unused kwargs in the field cases are necessary for generic code to compile.

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
    gram(M::MatRingElem)

Return the Gram matrix of $M$, i.e. if $M$ is an $r\times c$ matrix return
the $r\times r$ matrix whose entries $i, j$ are the dot products of the
$i$-th and $j$-th rows, respectively.
"""
function gram(M::MatRingElem)
   n = degree(M)
   z = similar(M)
   for i = 1:n
      for j = 1:n
         z[i, j] = zero(base_ring(M))
         for k = 1:n
            z[i, j] += M[i, k] * M[j, k]
         end
      end
   end
   return z
end

###############################################################################
#
#   Delegations to underlying matrix
#
###############################################################################

det(M::MatRingElem) = det(matrix(M))
rank(M::MatRingElem) = rank(matrix(M))

is_symmetric(M::MatRingElem) = is_symmetric(matrix(M))
is_alternating(M::MatRingElem) = is_alternating(matrix(M))
is_skew_symmetric(M::MatRingElem) = is_skew_symmetric(matrix(M))

is_upper_triangular(M::MatRingElem) = is_upper_triangular(matrix(M))
is_lower_triangular(M::MatRingElem) = is_lower_triangular(matrix(M))
is_diagonal(M::MatRingElem) = is_diagonal(matrix(M))

function kronecker_product(x::MatRingElem{T}, y::MatRingElem{T}) where {T <: RingElement}
  return Generic.MatRingElem(kronecker_product(matrix(x), matrix(y)))
end

function tr(x::MatRingElem{T}) where T <: NCRingElement
  return tr(matrix(x))
end

function map_entries!(f::S, dst::MatRingElem{T}, src::MatRingElem{U}) where {S, T <: NCRingElement, U <: NCRingElement}
  map_entries!(f, matrix(dst), matrix(src))
  return dst
end

function map_entries(f::S, M::MatrixElem{T}) where {S, T <: NCRingElement}
  return Generic.MatRingElem(map_entries(f, matrix(M)))
end

function pseudo_inv(M::MatRingElem{T}) where {T <: RingElement}
  X,d = pseudo_inv(matrix(M))
  return Generic.MatRingElem(X), d
end

function Base.inv(M::MatRingElem{T}) where {T <: RingElement}
  return Generic.MatRingElem(inv(matrix(M)))
end

function is_invertible_with_inverse(M::MatRingElem{T}; side::Symbol = :left) where {T <: RingElement}
  flag, inv = is_invertible_with_inverse(matrix(M); side = side)
  return flag, Generic.MatRingElem(inv)
end

is_invertible(M::MatRingElem{T}) where {T <: RingElement} = is_unit(det(M))

is_invertible(M::MatRingElem{T}) where {T <: FieldElement} = ncols(M) == rank(M)

function is_nilpotent(M::MatRingElem{T}) where {T <: RingElement}
  return is_nilpotent(matrix(M))
end

function hessenberg!(M::MatRingElem{T}) where {T <: RingElement}
  return Generic.MatRingElem(hessenberg!(matrix(M)))
end

function hessenberg(M::MatRingElem{T}) where {T <: RingElement}
  return Generic.MatRingElem(hessenberg(matrix(M)))
end

function is_hessenberg(M::MatRingElem{T}) where {T <: RingElement}
  return is_hessenberg(matrix(M))
end

function charpoly_hessenberg!(S::Ring, M::MatRingElem{T}) where {T <: RingElement}
  a = matrix(M)
  return charpoly_hessenberg!(S, a)  ## !!! WARNING !!!  may not be correct
end

function charpoly(S::PolyRing{T}, M::MatRingElem{T}) where {T <: RingElement}
  return charpoly(S, matrix(M))
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

# Sampler for a MatRing not needing arguments (passed via make)
# this allows to obtain the Sampler in simple cases without having to know about make
# (when one can do `rand(M)`, one can expect to be able to do `rand(Sampler(rng, M))`)
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

################################################################################
#
#  Promotion
#
################################################################################

function Base.promote(x::MatRingElem{S},
                      y::MatRingElem{T}) where {S <: NCRingElement,
                                                T <: NCRingElement}
   U = promote_rule_sym(S, T)
   if U === S
      return x, change_base_ring(base_ring(x), y)
   elseif U === T
      return change_base_ring(base_ring(y), x), y
   else
      error("Cannot promote to common type")
   end
end

# matrix * vec and vec * matrx
function Base.promote(x::MatRingElem{S},
                      y::Vector{T}) where {S <: NCRingElement,
                                           T <: NCRingElement}
   U = promote_rule_sym(S, T)
   if U === S
      return x, map(base_ring(x), y)::Vector{S}  # Julia needs help here
   elseif U === T && length(y) != 0
      return change_base_ring(parent(y[1]), x), y
   else
      error("Cannot promote to common type")
   end
end

function Base.promote(x::Vector{S},
                      y::MatRingElem{T}) where {S <: NCRingElement,
                                                T <: NCRingElement}
   yy, xx = promote(y, x)
   return xx, yy
end

*(x::MatRingElem, y::Vector) = *(promote(x, y)...)

*(x::Vector, y::MatRingElem) = *(promote(x, y)...)

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
  @req (n >= 0)  "Matrix dimension must be non-negative"
  @req (n < 2^30)  "Matrix dimension is excessively large"
  R = base_ring(M)
  return Generic.MatRingElem(identity_matrix(R,n))
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
    matrix_ring(R::NCRing, n::Int)

Return the matrix algebra (or matrix ring) of degree $n$ over the base ring $R$.

The returned parent object represents the ring of all $n \times n$ matrices
over $R$.

# Examples

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S = matrix_ring(R, 3)
Matrix ring of degree 3
  over univariate polynomial ring in t over rationals
```
"""
function matrix_ring(R::NCRing, n::Int)
   Generic.matrix_ring(R, n)
end
