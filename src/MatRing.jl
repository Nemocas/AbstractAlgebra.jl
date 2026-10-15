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

base_ring_type(::Type{MatRing{T}}) where T <: NCRingElement = parent_type(T)

base_ring_type(::Type{MatRingElem{T}}) where T <: NCRingElement = parent_type(T)

parent_type(::Type{MatRingElem{T}}) where T <: NCRingElement = MatRing{T}

elem_type(::Type{MatRing{T}}) where {T <: NCRingElement} = MatRingElem{T}

base_ring(a::MatRing{T}) where {T <: NCRingElement} = a.base_ring::parent_type(T)

base_ring(a::MatRingElem{T}) where {T <: NCRingElement} = base_ring(matrix(a))

matrix(A::MatRingElem{T}) where {T <: NCRingElement} = A.data::dense_matrix_type(T)

@doc raw"""
    parent(A::MatRingElem{T}) where T <: NCRingElement

Return the matrix algebra containing `A`.

This is the matrix algebra over the base ring of `A` whose degree is the
number of rows (equivalently columns) of `A`.

# Examples

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
parent(A::MatRingElem{T}) where T <: NCRingElement = MatRing{T}(base_ring(A), nrows(matrix(A)))

is_exact_type(::Type{MatRingElem{T}}) where T <: NCRingElement = is_exact_type(T)

is_domain_type(::Type{MatRingElem{T}}) where T <: NCRingElement = false

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::MatRingElem, h::UInt)
  b = 0x6413942b83a26c65 % UInt
  return xor(b, hash(matrix(a), h))
end

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
    degree(A::MatRingElem{T}) where T <: NCRingElement

Return the degree $n$ of the parent matrix algebra of `A`.

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
degree(A::MatRingElem{T}) where T <: NCRingElement = degree(parent(A))

zero(R::MatRing) = R()

one(R::MatRing) = R(1)

is_unit(A::MatRingElem{T}) where T <: RingElement = is_unit(det(A))

is_unit(A::MatRingElem{T}) where T <: FieldElement = rank(A) == degree(A)

# proof over a commutative ring: use adj(A)*A = det(A)*I = A*adj(A)
is_zero_divisor(A::MatRingElem{T}) where T <: RingElement = is_zero_divisor(det(A))

is_zero_divisor(A::MatRingElem{T}) where T <: FieldElement = rank(A) != degree(A)

function is_zero_divisor_with_annihilator(A::MatRingElem{T}) where T <: RingElement
   f, b = is_zero_divisor_with_annihilator(det(A))
   throw(NotImplementedError(:adj, A)) #return f, b*adj(A)
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
    similar(A::MatRingElem, R::NCRing, n::Int)
    similar(A::MatRingElem, R::NCRing)
    similar(A::MatRingElem, n::Int)
    similar(A::MatRingElem)

Create an uninitialized matrix ring element over the given ring and dimension,
with defaults based upon the given source matrix ring element `A`.
"""
function similar(A::MatRingElem, R::NCRing=base_ring(A), n::Int=degree(A))
   @req (n >= 0)  "Matrix dimension must be non-negative"
   @req (n < 2^30)  "Matrix dimension is excessively large"
   return MatRingElem(R, n, fill(0,n^2)) # n^2 cannot overflow given check in line above
end

similar(A::MatRingElem, n::Int) = similar(A, base_ring(A), n)

@doc raw"""
    zero(A::MatRingElem, R::NCRing, n::Int)
    zero(A::MatRingElem, R::NCRing)
    zero(A::MatRingElem, n::Int)
    zero(A::MatRingElem)

Create a zero matrix ring element over the given ring and dimension,
with defaults based upon the given source matrix ring element `A`.
"""
zero(A::MatRingElem, R::NCRing=base_ring(A), n::Int=degree(A)) = zero!(similar(A, R, n))
zero(A::MatRingElem, n::Int) = zero!(similar(A, n))

iszero(A::MatRingElem{T}) where T <: NCRingElement = iszero(matrix(A))

one(A::MatRingElem{T}) where T <: NCRingElement = one(parent(A))

isone(A::MatRingElem{T}) where T <: NCRingElement = isone(matrix(A))

@inline is_zero_entry(A::MatRingElem, i::Int, j::Int) = is_zero_entry(matrix(A), i, j)

@inline is_positive_entry(A::MatRingElem, i::Int, j::Int) = is_positive_entry(matrix(A), i, j)

@inline is_negative_entry(A::MatRingElem, i::Int, j::Int) = is_negative_entry(matrix(A), i, j)

Base.@propagate_inbounds is_zero_row(A::MatRingElem, i::Int) = is_zero_row(matrix(A), i)

Base.@propagate_inbounds is_zero_column(A::MatRingElem, j::Int) = is_zero_column(matrix(A), j)

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(A::MatRingElem{T}) where T <: NCRingElement = canonical_unit(matrix(A))

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

function mutable_copy(d::MatRingElem{T}) where T <: NCRingElement
   z = similar(d)
   for i = 1:nrows(d)
      for j = 1:ncols(d)
         z[i, j] = mutable_copy(d[i, j])
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

-(x::MatRingElem) = MatRingElem(-matrix(x))

*(x::T, y::T) where {T <: MatRingElem} = MatRingElem(matrix(x) * matrix(y))

+(x::T, y::T) where {T <: MatRingElem} = MatRingElem(matrix(x) + matrix(y))

-(x::T, y::T) where {T <: MatRingElem} = MatRingElem(matrix(x) - matrix(y))

^(a::MatRingElem{T}, b::Int) where T <: NCRingElement = MatRingElem(matrix(a)^b)

==(x::T, y::T) where {T <: MatRingElem} = matrix(x) == matrix(y)

isequal(x::T, y::T) where {T <: MatRingElem} = isequal(matrix(x), matrix(y))

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(x::JuliaRingElement, y::MatRingElem{T}) where T <: NCRingElement
  return MatRingElem(x * matrix(y))
end

function *(x::T, y::MatRingElem{T}) where {T <: NCRingElem}
  return MatRingElem(x * matrix(y))
end

function *(x::MatRingElem{T}, y::JuliaRingElement) where T <: NCRingElement
  return MatRingElem(matrix(x) * y)
end

function *(x::MatRingElem{T}, y::T) where {T <: NCRingElem}
  return MatRingElem(matrix(x) * y)
end

function +(x::JuliaRingElement, y::MatRingElem{T}) where T <: NCRingElement
  return MatRingElem(x + matrix(y))
end

function +(x::T, y::MatRingElem{T}) where {T <: NCRingElem}
  return MatRingElem(x + matrix(y))
end

function +(x::MatRingElem{T}, y::JuliaRingElement) where T <: NCRingElement
  return MatRingElem(matrix(x) + y)
end

function +(x::MatRingElem{T}, y::T) where {T <: NCRingElem}
  return MatRingElem(matrix(x) + y)
end

function -(x::JuliaRingElement, y::MatRingElem{T}) where T <: NCRingElement
  return MatRingElem(x - matrix(y))
end

function -(x::T, y::MatRingElem{T}) where {T <: NCRingElem}
  return MatRingElem(x - matrix(y))
end

function -(x::MatRingElem{T}, y::JuliaRingElement) where T <: NCRingElement
  return MatRingElem(matrix(x) - y)
end

function -(x::MatRingElem{T}, y::T) where {T <: NCRingElem}
  return MatRingElem(matrix(x) - y)
end

*(P::Perm, A::MatRingElem) = MatRingElem(P * matrix(A))

*(A::MatRingElem, P::Perm) = MatRingElem(matrix(A) * P)

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

function divexact(x::MatRingElem{T}, y::JuliaRingElement; check::Bool=true) where T <: NCRingElement
   return MatRingElem(divexact(matrix(x), y; check=check))
end

function divexact(x::MatRingElem{T}, y::T; check::Bool=true) where {T <: RingElem}
   return MatRingElem(divexact(matrix(x), y; check=check))
end

function divexact_left(x::MatRingElem{T}, y::T; check::Bool=true) where {T <: NCRingElem}
   return MatRingElem(divexact_left(matrix(x), y; check=check))
end

function divexact_right(x::MatRingElem{T}, y::T; check::Bool=true) where {T <: NCRingElem}
   return MatRingElem(divexact_right(matrix(x), y; check=check))
end

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
    gram(A::MatRingElem)

Return the Gram matrix of $A$, i.e. if $A$ is an $r\times c$ matrix return
the $r\times r$ matrix whose entries $i, j$ are the dot products of the
$i$-th and $j$-th rows, respectively.
"""
function gram(A::MatRingElem)
   n = degree(A)
   z = similar(A)
   for i = 1:n
      for j = 1:n
         z[i, j] = zero(base_ring(A))
         for k = 1:n
            z[i, j] += A[i, k] * A[j, k]
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

det(A::MatRingElem) = det(matrix(A))
rank(A::MatRingElem) = rank(matrix(A))

is_symmetric(A::MatRingElem) = is_symmetric(matrix(A))
is_alternating(A::MatRingElem) = is_alternating(matrix(A))
is_skew_symmetric(A::MatRingElem) = is_skew_symmetric(matrix(A))

is_upper_triangular(A::MatRingElem) = is_upper_triangular(matrix(A))
is_lower_triangular(A::MatRingElem) = is_lower_triangular(matrix(A))
is_diagonal(A::MatRingElem) = is_diagonal(matrix(A))

function kronecker_product(x::MatRingElem{T}, y::MatRingElem{T}) where {T <: RingElement}
  return MatRingElem(kronecker_product(matrix(x), matrix(y)))
end

function tr(x::MatRingElem{T}) where T <: NCRingElement
  return tr(matrix(x))
end

(S::MatSpace{T})(A::MatRingElem{T}) where {T <: NCRingElement} = S(matrix(A))

content(A::MatRingElem{T}) where T <: RingElement = content(matrix(A))

function change_base_ring(R::NCRing, A::MatRingElem{T}) where {T <: NCRingElement}
  return MatRingElem(change_base_ring(R, matrix(A)))
end

function map_entries!(f::S, dst::MatRingElem{T}, src::MatRingElem{U}) where {S, T <: NCRingElement, U <: NCRingElement}
  map_entries!(f, matrix(dst), matrix(src))
  return dst
end

Base.map!(f::S, dst::MatRingElem{T}, src::MatRingElem{U}) where {S, T <: NCRingElement, U <: NCRingElement} = map_entries!(f, dst, src)

function map_entries(f::S, A::MatRingElem{T}) where {S, T <: NCRingElement}
  return MatRingElem(map_entries(f, matrix(A)))
end

Base.map(f::S, A::MatRingElem{T}) where {S, T <: NCRingElement} = map_entries(f, A)

function pseudo_inv(A::MatRingElem{T}) where {T <: RingElement}
  X,d = pseudo_inv(matrix(A))
  return MatRingElem(X), d
end

function Base.inv(A::MatRingElem{T}) where {T <: RingElement}
  return MatRingElem(inv(matrix(A)))
end

function is_invertible_with_inverse(A::MatRingElem{T}; side::Symbol = :left) where {T <: RingElement}
  flag, inv = is_invertible_with_inverse(matrix(A); side = side)
  return flag, MatRingElem(inv)
end

is_invertible(A::MatRingElem{T}) where {T <: RingElement} = is_unit(det(A))

is_invertible(A::MatRingElem{T}) where {T <: FieldElement} = ncols(A) == rank(A)

function is_nilpotent(A::MatRingElem{T}) where {T <: RingElement}
  return is_nilpotent(matrix(A))
end

function hessenberg!(A::MatRingElem{T}) where {T <: RingElement}
  hessenberg!(matrix(A))
  return A
end

function hessenberg(A::MatRingElem{T}) where {T <: RingElement}
  return MatRingElem(hessenberg(matrix(A)))
end

function is_hessenberg(A::MatRingElem{T}) where {T <: RingElement}
  return is_hessenberg(matrix(A))
end

function charpoly_hessenberg!(S::Ring, A::MatRingElem{T}) where {T <: RingElement}
  a = matrix(A)
  return charpoly_hessenberg!(S, a)  ## !!! WARNING !!!  may not be correct
end

function charpoly(S::PolyRing{T}, A::MatRingElem{T}) where {T <: RingElement}
  return charpoly(S, matrix(A))
end

function lu(A::MatRingElem{T}, P = SymmetricGroup(nrows(A))) where {T <: FieldElement}
  S = parent(A)
  r, p, L, U = lu(matrix(A), P)
  return r, p, S(L), S(U)
end

function fflu(A::MatRingElem{T}, P = SymmetricGroup(nrows(A))) where {T <: RingElement}
  S = parent(A)
  r, d, p, L, U = fflu(matrix(A), P)
  return r, d, p, S(L), S(U)
end

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

function _can_solve_with_solution_lu(A::MatRingElem{T}, B::MatRingElem{T}) where {T <: RingElement}
   check_parent(A, B)
   R = base_ring(A)
   MS = matrix(A)
   BS = matrix(B)
   flag, S = _can_solve_with_solution_lu(MS, BS)
   SA = MatRingElem(S)
   return flag, SA
end

function can_solve_with_solution(A::MatRingElem{T}, B::MatRingElem{T}) where {T <: RingElement}
   check_parent(A, B)
   R = base_ring(A)
   MS = matrix(A)
   BS = matrix(B)
   flag, S = can_solve_with_solution(MS, BS)
   SA = MatRingElem(S)
   return flag, SA
end

function _can_solve_with_solution_fflu(A::MatRingElem{T}, B::MatRingElem{T}) where {T <: RingElement}
   check_parent(A, B)
   R = base_ring(A)
   MS = matrix(A)
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
    minpoly(S::Ring, A::MatRingElem{T}) where {T <: RingElement}

Return the minimal polynomial $p$ of the matrix $A$. The polynomial ring $S$
of the resulting polynomial must be supplied and the matrix must be square.
"""
function minpoly(S::Ring, A::MatRingElem{T}, charpoly_only::Bool = false) where {T <: RingElement}
   return minpoly(S, matrix(A), charpoly_only)
end

function minpoly(A::MatRingElem{T}, charpoly_only::Bool = false) where {T <: RingElement}
   R = base_ring(A)
   Rx, x = polynomial_ring(R; cached=false)
   return minpoly(Rx, A, charpoly_only)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(A::MatRingElem)
   return MatRingElem(zero!(matrix(A)))
end

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

function mul!(z::Vector{T}, x::MatRingElem{T}, y::Vector{T}) where T <: NCRingElement
  return mul!(z, matrix(x), y)
end

function mul!(z::Vector{T}, x::Vector{T}, y::MatRingElem{T}) where T <: NCRingElement
  return mul!(z, x, matrix(y))
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

function identity_matrix(A::MatRingElem{T}, n::Int) where T <: NCRingElement
  @req (n >= 0)  "Matrix dimension must be non-negative"
  @req (n < 2^30)  "Matrix dimension is excessively large"
  R = base_ring(A)
  return MatRingElem(identity_matrix(R,n))
end

@doc raw"""
    identity_matrix(A::MatRingElem{T}) where T <: RingElement

Return the identity matrix over the same base ring as $A$ and with the
same dimensions.
"""
function identity_matrix(A::MatRingElem{T}) where T <: NCRingElement
   return identity_matrix(A, nrows(A))
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
    (a::MatRing{T})(b::MatElem{S}) where {S <: NCRingElement, T <: NCRingElement}
    (a::MatRing{T})(b::Matrix{S}) where {S <: NCRingElement, T <: NCRingElement}
    (a::MatRing{T})(b::Vector{S}) where {S <: NCRingElement, T <: NCRingElement}

Construct an element of the matrix algebra `a`.

The call `a()` returns the zero matrix in `a`.

If `b` is a ring element coercible into the base ring of `a`, then
`a(b)` returns the scalar matrix in `a` whose diagonal entries are `b`.

If `b` is a matrix algebra element whose parent is `a`, then `a(b)`
returns `b` unchanged.

If `b` is a matrix whose dimensions and base ring agree with
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

function (a::MatRing{T})(b::MatElem{S}) where {S <: NCRingElement, T <: NCRingElement}
   R = base_ring(a)
   _check_dim(nrows(a), ncols(a), b)
   z = MatRingElem(matrix(R, b))
   return z
end

(a::MatRing{T})(b::MatRingElem{S}) where {S <: NCRingElement, T <: NCRingElement} = a(matrix(b))

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
function matrix_ring(R::NCRing, n::Int; cached::Bool = true)
   # TODO: the 'cached' argument is ignored and mainly here for backwards compatibility
   # (and perhaps future compatibility, in case we need it again)
   @req n >= 0 "n must be a non-negative integer"
   T = elem_type(R)
   return MatRing{T}(R, n)
end
