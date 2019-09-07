###############################################################################
#
#   Matrix.jl : Generic mxn matrices over rings
#
###############################################################################

export MatrixSpace, fflu!, fflu, solve_triu, isrref, charpoly_danilevsky!,
       charpoly_danilevsky_ff!, hessenberg!, hessenberg, ishessenberg,
       identity_matrix, charpoly_hessenberg!, invert_cols, invert_cols!,
       invert_rows, invert_rows!, matrix, minpoly, typed_hvcat, typed_hcat,
       powers, randmat_triu, randmat_with_rank, similarity!, solve,
       solve_rational, hnf, hnf_kb, hnf_kb_with_transform, hnf_with_transform,
       issquare, snf, snf_with_transform, weak_popov,
       weak_popov_with_transform, can_solve_left_reduced_triu,
       extended_weak_popov, extended_weak_popov_with_transform, rank,
       rank_profile_popov, hnf_via_popov, hnf_via_popov_with_transform, popov,
       popov_with_transform, det_popov, _check_dim, nrows, ncols, gram, rref,
       rref!, swap_cols, swap_cols!, swap_rows, swap_rows!, hnf_kb,
       hnf_kb_with_transform, hnf_cohen, hnf_cohen_with_transform, snf_kb,
       snf_kb_with_transform, find_pivot_popov, inv!, zero_matrix,
       kronecker_product, minors, tr, lu, lu!, pseudo_inv, dense_matrix_type,
       kernel, iszero_row, iszero_column, left_kernel, right_kernel, ishnf,
       solve_left

###############################################################################
#
#   Similar and eye
#
###############################################################################

function _similar(x::Union{Mat{T}, MatAlgElem{T}}, R::Ring, r::Int, c::Int) where T <: RingElement
   TT = elem_type(R)
   M = similar(x.entries, TT, r, c)
   for i in 1:size(M, 1)
      for j in 1:size(M, 2)
         M[i, j] = zero(R)
      end
   end
   z = x isa Mat ? MatSpaceElem{TT}(M) : MatAlgElem{TT}(M)
   z.base_ring = R
   return z
end

similar(x::Mat, R::Ring, r::Int, c::Int) = _similar(x, R, r, c)

similar(x::Mat, R::Ring=base_ring(x))= similar(x, R, nrows(x), ncols(x))

similar(x::Mat, r::Int, c::Int) = similar(x, base_ring(x), r, c)

@doc Markdown.doc"""
    eye(x::Generic.MatrixElem)
> Return the identity matrix with the same shape as $x$.
"""
function eye(x::MatrixElem)
  z = similar(x)
  for i in 1:nrows(x)
    z[i, i] = one(base_ring(x))
  end
  return z
end

@doc Markdown.doc"""
    eye(x::Generic.MatrixElem, d::Int)
> Return the $d$-by-$d$ identity matrix with the same base ring as $x$.
"""
function eye(x::MatrixElem, d::Int)
  z = similar(x, d, d)
  for i in 1:nrows(z)
    z[i, i] = one(base_ring(x))
  end
  return z
end

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{S}) where {T <: RingElement, S <: Mat{T}} = MatSpace{T}

elem_type(::Type{MatSpace{T}}) where {T <: RingElement} = MatSpaceElem{T}

@doc Markdown.doc"""
    base_ring(a::AbstractAlgebra.MatSpace{T}) where {T <: RingElement}
> Return the base ring $R$ of the given matrix space.
"""
base_ring(a::AbstractAlgebra.MatSpace{T}) where {T <: RingElement} = a.base_ring::parent_type(T)

@doc Markdown.doc"""
    base_ring(a::Generic.MatrixElem{T}) where {T <: RingElement}
> Return the base ring $R$ of the matrix space that the supplied matrix $r$
> belongs to.
"""
base_ring(a::MatrixElem{T}) where {T <: RingElement} = a.base_ring::parent_type(T)

@doc Markdown.doc"""
    parent(a::AbstractAlgebra.MatElem{T}, cached::Bool = true) where T <: RingElement
> Return the parent object of the given matrix.
"""
parent(a::AbstractAlgebra.MatElem{T}, cached::Bool = true) where T <: RingElement =
    MatSpace{T}(a.base_ring, size(a.entries)..., cached)

dense_matrix_type(::Type{T}) where T <: RingElement = MatSpaceElem{T}

@doc Markdown.doc"""
    dense_matrix_type(R::Ring)
> Return the type of matrices over the given ring.
"""
dense_matrix_type(R::Ring) = dense_matrix_type(elem_type(R))

function check_parent(a::AbstractAlgebra.MatElem, b::AbstractAlgebra.MatElem, throw::Bool = true)
  fl = (base_ring(a) != base_ring(b) || nrows(a) != nrows(b) || ncols(a) != ncols(b))
  fl && throw && error("Incompatible matrix spaces in matrix operation")
  return !fl
end

function _check_dim(r::Int, c::Int, arr::AbstractArray{T, 2}, transpose::Bool = false) where {T}
  if !transpose
    size(arr) != (r, c) && throw(ErrorConstrDimMismatch(r, c, size(arr)...))
  else
    size(arr) != (c, r) && throw(ErrorConstrDimMismatch(r, c, (reverse(size(arr)))...))
  end
  return nothing
end

function _check_dim(r::Int, c::Int, arr::AbstractArray{T, 1}) where {T}
  length(arr) != r*c && throw(ErrorConstrDimMismatch(r, c, length(arr)))
  return nothing
end

function _checkbounds(i::Int, j::Int)
   j >= 1 && j <= i
end

function _checkbounds(A, i::Int, j::Int)
  (_checkbounds(nrows(A), i) && _checkbounds(ncols(A), j)) ||
            Base.throw_boundserror(A, (i, j))
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

@doc Markdown.doc"""
    nrows(a::Generic.MatrixSpace)
> Return the number of rows of the given matrix space.
"""
nrows(a::MatSpace) = a.nrows

@doc Markdown.doc"""
    ncols(a::Generic.MatrixSpace)
> Return the number of columns of the given matrix space.
"""
ncols(a::MatSpace) = a.ncols

function Base.hash(a::AbstractAlgebra.MatElem, h::UInt)
   b = 0x3e4ea81eb31d94f4%UInt
   for i in 1:nrows(a)
      for j in 1:ncols(a)
         b = xor(b, xor(hash(a[i, j], h), h))
         b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
      end
   end
   return b
end

@doc Markdown.doc"""
    nrows(a::Generic.MatrixElem)
> Return the number of rows of the given matrix.
"""
nrows(a::MatrixElem) = size(a.entries, 1)

@doc Markdown.doc"""
    ncols(a::Generic.MatrixElem)
> Return the number of columns of the given matrix.
"""
ncols(a::MatrixElem) = size(a.entries, 2)

@doc Markdown.doc"""
    length(a::Generic.MatrixElem)
> Return the number of entries in the given matrix.
"""
length(a::MatrixElem) = nrows(a) * ncols(a)

@doc Markdown.doc"""
    isempty(a::Generic.MatrixElem)
> Return `true` if `a` does not contain any entry (i.e. `length(a) == 0`), and `false` otherwise.
"""
isempty(a::MatrixElem) = (nrows(a) == 0) | (ncols(a) == 0)

Base.@propagate_inbounds function getindex(a::MatrixElem, r::Int, c::Int)
   return a.entries[r, c]
end

Base.@propagate_inbounds function setindex!(a::MatrixElem, d::T, r::Int,
                                            c::Int) where T <: RingElement
    a.entries[r, c] = base_ring(a)(d)
end

@doc Markdown.doc"""
    zero(a::AbstractAlgebra.MatSpace)
> Construct the zero matrix in the given matrix space.
"""
zero(a::AbstractAlgebra.MatSpace) = a()

@doc Markdown.doc"""
    one(a::AbstractAlgebra.MatSpace)
> Construct the matrix in the given matrix space with ones down the diagonal
> and zeroes elsewhere.
"""
one(a::AbstractAlgebra.MatSpace) = a(1)

@doc Markdown.doc"""
    iszero(a::Generic.MatrixElem)
> Return `true` if the supplied matrix $a$ is the zero matrix, otherwise
> return `false`.
"""
function iszero(a::MatrixElem)
   for i = 1:nrows(a)
      for j = 1:ncols(a)
         if !iszero(a[i, j])
            return false
         end
      end
  end
  return true
end

@doc Markdown.doc"""
    isone(a::Generic.MatrixElem)
> Return `true` if the supplied matrix $a$ is diagonal with ones along the
> diagonal, otherwise return `false`.
"""
function isone(a::MatrixElem)
   for i = 1:nrows(a)
      for j = 1:ncols(a)
         if i == j
            if !isone(a[i, j])
               return false
            end
         else
            if !iszero(a[i, j])
               return false
            end
         end
      end
  end
  return true
end

@doc Markdown.doc"""
    iszero_row(M::MatrixElem{T}, i::Int) where T <: RingElement
> Return `true` if the $i$-th row of the matrix $M$ is zero.
"""
function iszero_row(M::MatrixElem{T}, i::Int) where T <: RingElement
  for j in 1:ncols(M)
    if !iszero(M[i, j])
      return false
    end
  end
  return true
end

@doc Markdown.doc"""
    iszero_column(M::MatrixElem{T}, i::Int) where T <: RingElement
> Return `true` if the $i$-th column of the matrix $M$ is zero.
"""
function iszero_column(M::MatrixElem{T}, i::Int) where T <: RingElement
  for j in 1:nrows(M)
    if !iszero(M[j, i])
      return false
    end
  end
  return true
end

function copy(d::MatrixElem)
   c = similar(d)
   for i = 1:nrows(d)
      for j = 1:ncols(d)
         c[i, j] = d[i, j]
      end
   end
   return c
end

function deepcopy_internal(d::MatrixElem, dict::IdDict)
   c = similar(d)
   for i = 1:nrows(d)
      for j = 1:ncols(d)
         c[i, j] = deepcopy(d[i, j])
      end
   end
   return c
end

function deepcopy_internal(d::MatSpaceView{T}, dict::IdDict) where T <: RingElement
   return MatSpaceView(deepcopy(d.entries), d.base_ring)
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::MatrixElem) = canonical_unit(a[1, 1])

###############################################################################
#
#   Sub
#
###############################################################################

function sub(M::AbstractAlgebra.MatElem, rows::UnitRange{Int}, cols::UnitRange{Int})
  Generic._checkbounds(M, rows.start, cols.start)
  Generic._checkbounds(M, rows.stop, cols.stop)
  z = similar(M, length(rows), length(cols))
  startr = first(rows)
  startc = first(cols)
  for i in rows
    for j in cols
      z[i - startr + 1, j - startc + 1] = deepcopy(M[i, j])
    end
  end
  return z
end

@doc Markdown.doc"""
    sub(M::AbstractAlgebra.MatElem, r1::Int, c1::Int, r2::Int, c2::Int)
> Return a copy of the submatrix of $M$ from $(r1, c1)$ to $(r2, c2)$ inclusive. Note
> that is the copy is modified, the original matrix is not.
"""
function sub(M::AbstractAlgebra.MatElem, r1::Int, c1::Int, r2::Int, c2::Int)
  return sub(M, r1:r2, c1:c2)
end

@doc Markdown.doc"""
    sub(M::AbstractAlgebra.MatElem, rows::Array{Int,1}, cols::Array{Int,1})
> Return a copy of the submatrix $A$ of $M$ defined by A[i,j] = M[rows[i], cols[j]]
> for i=1,...,length(rows) and j=1,...,length(cols)
"""
function sub(M::AbstractAlgebra.MatElem, rows::Array{Int,1}, cols::Array{Int,1})
   z = similar(M, length(rows), length(cols))
   for i in 1:length(rows)
      for j in 1:length(cols)
         Generic._checkbounds(M, rows[i], cols[j])
         z[i, j] = deepcopy(M[rows[i], cols[j]])
      end
   end
   return z
end

getindex(x::AbstractAlgebra.MatElem, r::UnitRange{Int}, c::UnitRange{Int}) = sub(x, r, c)

getindex(x::AbstractAlgebra.MatElem, r::UnitRange, ::Colon) = sub(x, r, 1:ncols(x))

getindex(x::AbstractAlgebra.MatElem, ::Colon, c::UnitRange{Int}) = sub(x, 1:nrows(x), c)

getindex(x::AbstractAlgebra.MatElem, ::Colon, ::Colon) = sub(x, 1:nrows(x), 1:ncols(x))

function Base.view(M::Mat{T}, rows::UnitRange{Int}, cols::UnitRange{Int}) where T <: RingElement
   return MatSpaceView(view(M.entries, rows, cols), M.base_ring)
end

function Base.view(M::AbstractAlgebra.MatElem{T}, rows::Colon, cols::UnitRange{Int64}) where T <: RingElement
   return view(M, 1:nrows(M), cols)
end

function Base.view(M::AbstractAlgebra.MatElem{T}, rows::UnitRange{Int64}, cols::Colon) where T <: RingElement
   return view(M, rows, 1:ncols(M))
end

function Base.view(M::AbstractAlgebra.MatElem{T}, rows::Colon, cols::Colon) where T <: RingElement
   return view(M, 1:nrows(M), 1:ncols(M))
end

################################################################################
#
#   Size
#
################################################################################

size(x::MatElem) = tuple(nrows(x), ncols(x))

size(t::MatElem, d) = d <= 2 ? size(t)[d] : 1

issquare(a::MatElem) = (nrows(a) == ncols(a))

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, a::AbstractAlgebra.MatSpace)
   print(io, "Matrix Space of ")
   print(io, a.nrows, " rows and ", a.ncols, " columns over ")
   print(IOContext(io, :compact => true), base_ring(a))
end

function show(io::IO, a::MatrixElem)
   isempty(a) && return print(io, "$r by $c matrix")

   r = nrows(a)
   c = ncols(a)

   # preprint each element to know the widths so as to align the columns
   strings = String[sprint(print, a[i,j], context = :compact => true) for i=1:r, j=1:c]
   maxs = maximum(length, strings, dims=1)

   for i = 1:r
      print(io, "[")
      for j = 1:c
         s = strings[i, j]
         s = ' '^(maxs[j] - length(s)) * s
         print(io, s)
         if j != c
            print(io, "  ")
         end
      end
      print(io, "]")
      if i != r
         println(io, "")
      end
   end
end

show_minus_one(::Type{AbstractAlgebra.MatElem{T}}) where {T <: RingElement} = false

###############################################################################
#
#   Unary operations
#
###############################################################################

@doc Markdown.doc"""
    -(x::Generic.MatrixElem)
> Return $-x$.
"""
function -(x::MatrixElem)
   z = similar(x)
   for i in 1:nrows(x)
      for j in 1:ncols(x)
         z[i, j] = -x[i, j]
      end
   end
   return z
end

###############################################################################
#
#   Binary operations
#
###############################################################################

@doc Markdown.doc"""
    +(x::Generic.MatrixElem{T}, y::Generic.MatrixElem{T}) where {T <: RingElement}
> Return $x + y$.
"""
function +(x::MatrixElem{T}, y::MatrixElem{T}) where {T <: RingElement}
   check_parent(x, y)
   r = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         r[i, j] = x[i, j] + y[i, j]
      end
   end
   return r
end

@doc Markdown.doc"""
    -(x::Generic.MatrixElem{T}, y::Generic.MatrixElem{T}) where {T <: RingElement}
> Return $x - y$.
"""
function -(x::MatrixElem{T}, y::MatrixElem{T}) where {T <: RingElement}
   check_parent(x, y)
   r = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         r[i, j] = x[i, j] - y[i, j]
      end
   end
   return r
end

@doc Markdown.doc"""
    *(x::AbstractAlgebra.MatElem{T}, y::AbstractAlgebra.MatElem{T}) where {T <: RingElement}
> Return $x\times y$.
"""
function *(x::AbstractAlgebra.MatElem{T}, y::AbstractAlgebra.MatElem{T}) where {T <: RingElement}
   ncols(x) != nrows(y) && error("Incompatible matrix dimensions")
   A = similar(x, nrows(x), ncols(y))
   C = base_ring(x)()
   for i = 1:nrows(x)
      for j = 1:ncols(y)
         A[i, j] = base_ring(x)()
         for k = 1:ncols(x)
            A[i, j] = addmul_delayed_reduction!(A[i, j], x[i, k], y[k, j], C)
         end
         A[i, j] = reduce!(A[i, j])
      end
   end
   return A
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

@doc Markdown.doc"""
    *(x::Union{Integer, Rational, AbstractFloat}, y::Generic.MatrixElem)
> Return $x\times y$.
"""
function *(x::Union{Integer, Rational, AbstractFloat}, y::MatrixElem)
   z = similar(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         z[i, j] = x*y[i, j]
      end
   end
   return z
end

@doc Markdown.doc"""
    *(x::T, y::Generic.MatrixElem{T}) where {T <: RingElem}
> Return $x\times y$.
"""
function *(x::T, y::MatrixElem{T}) where {T <: RingElem}
   z = similar(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         z[i, j] = x*y[i, j]
      end
   end
   return z
end

@doc Markdown.doc"""
    *(x::Generic.MatrixElem, y::Union{Integer, Rational, AbstractFloat})
> Return $x\times y$.
"""
*(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat}) = y*x

@doc Markdown.doc"""
    *(x::Generic.MatrixElem{T}, y::T) where {T <: RingElem}
> Return $x\times y$.
"""
*(x::MatrixElem{T}, y::T) where {T <: RingElem} = y*x

@doc Markdown.doc"""
    +(x::Union{Integer, Rational, AbstractFloat}, y::Generic.MatrixElem)
> Return $S(x) + y$ where $S$ is the parent of $y$.
"""
function +(x::Union{Integer, Rational, AbstractFloat}, y::MatrixElem)
   z = similar(y)
   R = base_ring(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         if i != j
            z[i, j] = deepcopy(y[i, j])
         else
            z[i, j] = y[i, j] + R(x)
         end
      end
   end
   return z
end

@doc Markdown.doc"""
    +(x::Generic.MatrixElem, y::Union{Integer, Rational, AbstractFloat})
> Return $x + S(y)$ where $S$ is the parent of $x$.
"""
+(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat}) = y + x

@doc Markdown.doc"""
    +(x::T, y::Generic.MatrixElem{T}) where {T <: RingElem}
> Return $S(x) + y$ where $S$ is the parent of $y$.
"""
function +(x::T, y::MatrixElem{T}) where {T <: RingElem}
   z = similar(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         if i != j
            z[i, j] = deepcopy(y[i, j])
         else
            z[i, j] = y[i, j] + x
         end
      end
   end
   return z
end

@doc Markdown.doc"""
    +(x::Generic.MatrixElem{T}, y::T) where {T <: RingElem}
> Return $x + S(y)$ where $S$ is the parent of $x$.
"""
+(x::MatrixElem{T}, y::T) where {T <: RingElem} = y + x

@doc Markdown.doc"""
    -(x::Union{Integer, Rational, AbstractFloat}, y::Generic.MatrixElem)
> Return $S(x) - y$ where $S$ is the parent of $y$.
"""
function -(x::Union{Integer, Rational, AbstractFloat}, y::MatrixElem)
   z = similar(y)
   R = base_ring(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         if i != j
            z[i, j] = -y[i, j]
         else
            z[i, j] = R(x) - y[i, j]
         end
      end
   end
   return z
end

@doc Markdown.doc"""
    -(x::Generic.MatrixElem, y::Union{Integer, Rational, AbstractFloat})
> Return $x - S(y)$, where $S$ is the parent of $x$.
"""
function -(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat})
   z = similar(x)
   R = base_ring(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if i != j
            z[i, j] = deepcopy(x[i, j])
         else
            z[i, j] = x[i, j] - R(y)
         end
      end
   end
   return z
end

@doc Markdown.doc"""
    -(x::T, y::Generic.MatrixElem{T}) where {T <: RingElem}
> Return $S(x) - y$ where $S$ is the parent of $y$.
"""
function -(x::T, y::MatrixElem{T}) where {T <: RingElem}
   z = similar(y)
   R = base_ring(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         if i != j
            z[i, j] = -y[i, j]
         else
            z[i, j] = x - y[i, j]
         end
      end
   end
   return z
end

@doc Markdown.doc"""
    -(x::Generic.MatrixElem{T}, y::T) where {T <: RingElem}
> Return $x - S(y)$, where $S$ is the parent of $a$.
"""
function -(x::MatrixElem{T}, y::T) where {T <: RingElem}
   z = similar(x)
   R = base_ring(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if i != j
            z[i, j] = deepcopy(x[i, j])
         else
            z[i, j] = x[i, j] - y
         end
      end
   end
   return z
end

###############################################################################
#
#   Powering
#
###############################################################################

Base.literal_pow(::typeof(^), x::T, ::Val{p}) where {p, T <: MatElem} = x^p

@doc Markdown.doc"""
    ^(a::Generic.MatrixElem, b::Int)
> Return $a^b$. We require $b \geq 0$ and that the matrix $a$ is square.
"""
function ^(a::MatrixElem, b::Int)
   b < 0 && throw(DomainError(b, "Negative exponent in power"))
   !issquare(a) && error("Incompatible matrix dimensions in power")
   # special case powers of x for constructing polynomials efficiently
   if b == 0
      return eye(a)
   elseif b == 1
      return deepcopy(a)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit != 0
         z = z*z
         if (UInt(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end

@doc Markdown.doc"""
    powers(a::Generic.MatrixElem, d::Int)
> Return an array of matrices $M$ wher $M[i + 1] = a^i$ for $i = 0..d$
"""
function powers(a::MatrixElem, d::Int)
   !issquare(a) && error("Dimensions do not match in powers")
   d <= 0 && throw(DomainError(d, "Negative dimension in powers"))
   A = Array{typeof(a)}(undef, d + 1)
   A[1] = eye(a)
   if d > 1
      c = a
      A[2] = a
      for i = 2:d
         c *= a
         A[i + 1] = c
      end
   end
   return A
end

###############################################################################
#
#   Comparisons
#
###############################################################################

@doc Markdown.doc"""
    ==(x::Generic.MatrixElem{T}, y::Generic.MatrixElem{T}) where {T <: RingElement}
> Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
> that power series to different precisions may still be arithmetically
> equal to the minimum of the two precisions.
"""
function ==(x::MatrixElem{T}, y::MatrixElem{T}) where {T <: RingElement}
   b = check_parent(x, y, false)
   !b && return false
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if x[i, j] != y[i, j]
            return false
         end
      end
   end
   return true
end

@doc Markdown.doc"""
    isequal(x::Generic.MatrixElem{T}, y::Generic.MatrixElem{T}) where {T <: RingElement}
> Return `true` if $x == y$ exactly, otherwise return `false`. This function is
> useful in cases where the entries of the matrices are inexact, e.g. power
> series. Only if the power series are precisely the same, to the same precision,
> are they declared equal by this function.
"""
function isequal(x::MatrixElem{T}, y::MatrixElem{T}) where {T <: RingElement}
   b = check_parent(x, y, false)
   !b && return false
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if !isequal(x[i, j], y[i, j])
            return false
         end
      end
   end
   return true
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

@doc Markdown.doc"""
    ==(x::Generic.MatrixElem, y::Union{Integer, Rational, AbstractFloat})
> Return `true` if $x == S(y)$ arithmetically, where $S$ is the parent of $x$,
> otherwise return `false`.
"""
function ==(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat})
   for i = 1:min(nrows(x), ncols(x))
      if x[i, i] != y
         return false
      end
   end
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if i != j && !iszero(x[i, j])
            return false
         end
      end
   end
   return true
end

@doc Markdown.doc"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::Generic.MatrixElem)
> Return `true` if $S(x) == y$ arithmetically, where $S$ is the parent of $y$,
> otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::MatrixElem) = y == x

@doc Markdown.doc"""
    ==(x::Generic.MatrixElem{T}, y::T) where {T <: RingElem}
> Return `true` if $x == S(y)$ arithmetically, where $S$ is the parent of $x$,
> otherwise return `false`.
"""
function ==(x::MatrixElem{T}, y::T) where {T <: RingElem}
   for i = 1:min(nrows(x), ncols(x))
      if x[i, i] != y
         return false
      end
   end
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         if i != j && !iszero(x[i, j])
            return false
         end
      end
   end
   return true
end

@doc Markdown.doc"""
    ==(x::T, y::Generic.MatrixElem{T}) where {T <: RingElem}
> Return `true` if $S(x) == y$ arithmetically, where $S$ is the parent of $y$,
> otherwise return `false`.
"""
==(x::T, y::MatrixElem{T}) where {T <: RingElem} = y == x

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact(x::Generic.MatrixElem, y::Union{Integer, Rational, AbstractFloat})
> Return $x/y$, i.e. the matrix where each of the entries has been divided by
> $y$. Each division is expected to be exact.
"""
function divexact(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat})
   z = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         z[i, j] = divexact(x[i, j], y)
      end
   end
   return z
end

@doc Markdown.doc"""
    divexact(x::Generic.MatrixElem{T}, y::T) where {T <: RingElem}
> Return $x/y$, i.e. the matrix where each of the entries has been divided by
> $y$. Each division is expected to be exact.
"""
function divexact(x::MatrixElem{T}, y::T) where {T <: RingElem}
   z = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         z[i, j] = divexact(x[i, j], y)
      end
   end
   return z
end

###############################################################################
#
#   Kronecker product
#
###############################################################################

function kronecker_product(x::MatElem{T}, y::MatElem{T}) where {T <: RingElement}
    base_ring(parent(x)) == base_ring(parent(y)) || error("Incompatible matrix spaces in matrix operation")
    z = similar(x, nrows(x)*nrows(y), ncols(x)*ncols(y))
    for ix in 1:nrows(x)
       ixr = (ix-1)*nrows(y)
       for jx in 1:ncols(x)
          jxc = (jx-1)*ncols(y)
          for iy in 1:nrows(y)
             for jy in 1:ncols(y)
               z[ixr+iy,jxc+jy] = x[ix,jx]*y[iy,jy]
             end
          end
       end
    end
    return z
end

###############################################################################
#
#   Transpose
#
###############################################################################

@doc Markdown.doc"""
    transpose(x::Mat)
> Return the transpose of the given matrix.
"""
function transpose(x::Mat)
   return matrix(base_ring(x), permutedims(x.entries, [2, 1]))
end

###############################################################################
#
#   Gram matrix
#
###############################################################################

@doc Markdown.doc"""
    gram(x::AbstractAlgebra.MatElem)
> Return the Gram matrix of $x$, i.e. if $x$ is an $r\times c$ matrix return
> the $r\times r$ matrix whose entries $i, j$ are the dot products of the
> $i$-th and $j$-th rows, respectively.
"""
function gram(x::AbstractAlgebra.MatElem)
   z = similar(x, nrows(x), nrows(x))
   for i = 1:nrows(x)
      for j = 1:nrows(x)
         z[i, j] = zero(base_ring(x))
         for k = 1:ncols(x)
            z[i, j] += x[i, k] * x[j, k]
         end
      end
   end
   return z
end

###############################################################################
#
#   Trace
#
###############################################################################

@doc Markdown.doc"""
    tr(x::Generic.MatrixElem)
> Return the trace of the matrix $a$, i.e. the sum of the diagonal elements. We
> require the matrix to be square.
"""
function tr(x::MatrixElem)
   !issquare(x) && error("Not a square matrix in trace")
   d = zero(base_ring(x))
   for i = 1:nrows(x)
      d = addeq!(d, x[i, i])
   end
   return d
end

###############################################################################
#
#   Content
#
###############################################################################

@doc Markdown.doc"""
    content(x::Generic.MatrixElem)
> Return the content of the matrix $a$, i.e. the greatest common divisor of all
> its entries, assuming it exists.
"""
function content(x::MatrixElem)
  d = zero(base_ring(x))
  for i = 1:nrows(x)
     for j = 1:ncols(x)
        d = gcd(d, x[i, j])
        if isone(d)
           return d
        end
     end
  end
  return d
end

###############################################################################
#
#   Permutation
#
###############################################################################

@doc Markdown.doc"""
    *(P::Generic.perm, x::Generic.MatrixElem)
> Apply the pemutation $P$ to the rows of the matrix $x$ and return the result.
"""
function *(P::Generic.perm, x::MatrixElem)
   z = similar(x)
   m = nrows(x)
   n = ncols(x)
   for i = 1:m
      for j = 1:n
         z[P[i], j] = x[i, j]
      end
   end
   return z
end

###############################################################################
#
#   LU factorisation
#
###############################################################################

function lu!(P::Generic.perm, A::MatrixElem{T}) where {T <: FieldElement}
   m = nrows(A)
   n = ncols(A)
   rank = 0
   r = 1
   c = 1
   R = base_ring(A)
   t = R()
   while r <= m && c <= n
      if c != 1
         # reduction of lower right square was delayed, reduce left col now
         for i = r:m
            A[i, c] = reduce!(A[i, c])
         end
      end
      if iszero(A[r, c])
         i = r + 1
         while i <= m
            if !iszero(A[i, c])
               for j = 1:n
                  A[i, j], A[r, j] = A[r, j], A[i, j]
               end
               P[r], P[i] = P[i], P[r]
               break
            end
            i += 1
         end
         if i > m
            c += 1
            continue
         end
      end
      rank += 1
      d = -inv(A[r, c])
      # reduction of lower right square was delayed, reduce top row now
      if c != 1
         for j = c + 1:n
            A[r, j] = reduce!(A[r, j])
         end
      end
      for i = r + 1:m
         q = A[i, c]*d
         for j = c + 1:n
            t = mul_red!(t, A[r, j], q, false)
            A[i, j] = addeq!(A[i, j], t)
         end
         A[i, c] = R()
         A[i, rank] = -q
      end
      r += 1
      c += 1
   end
   inv!(P)
   return rank
end

@doc Markdown.doc"""
    lu(A::Generic.MatrixElem{T}, P = PermGroup(rows(A))) where {T <: FieldElement}
> Return a tuple $r, p, L, U$ consisting of the rank of $A$, a permutation
> $p$ of $A$ belonging to $P$, a lower triangular matrix $L$ and an upper
> triangular matrix $U$ such that $p(A) = LU$, where $p(A)$ stands for the
> matrix whose rows are the given permutation $p$ of the rows of $A$.
"""
function lu(A::MatrixElem{T}, P = PermGroup(nrows(A))) where {T <: FieldElement}
   m = nrows(A)
   n = ncols(A)
   P.n != m && error("Permutation does not match matrix")
   p = P()
   R = base_ring(A)
   U = deepcopy(A)
   L = similar(A, m, m)
   rank = lu!(p, U)
   for i = 1:m
      for j = 1:n
         if i > j
            L[i, j] = U[i, j]
            U[i, j] = R()
         elseif i == j
            L[i, j] = R(1)
         elseif j <= m
            L[i, j] = R()
         end
      end
   end
   return rank, p, L, U
end

function fflu!(P::Generic.perm, A::MatrixElem{T}) where {T <: RingElement}
   if !isdomain_type(T)
      error("Not implemented")
   end
   m = nrows(A)
   n = ncols(A)
   rank = 0
   r = 1
   c = 1
   R = base_ring(A)
   d = R(1)
   d2 = R(1)
   if m == 0 || n == 0
      return 0, d
   end
   t = R()
   while r <= m && c <= n
      if iszero(A[r, c])
         i = r + 1
         while i <= m
            if !iszero(A[i, c])
               for j = 1:n
                  A[i, j], A[r, j] = A[r, j], A[i, j]
               end
               P[r], P[i] = P[i], P[r]
               break
            end
            i += 1
         end
         if i > m
            c += 1
            continue
         end
      end
      rank += 1
      q = -A[r, c]
      for i = r + 1:m
         for j = c + 1:n
            A[i, j] = mul_red!(A[i, j], A[i, j], q, false)
            t = mul_red!(t, A[i, c], A[r, j], false)
            A[i, j] = addeq!(A[i, j], t)
            A[i, j] = reduce!(A[i, j])
            if r > 1
               A[i, j] = divexact(A[i, j], d)
            else
               A[i, j] = -A[i, j]
            end
         end
      end
      d = -A[r, c]
      d2 = A[r, c]
      r += 1
      c += 1
   end
   inv!(P)
   return rank, d2
end

function fflu!(P::Generic.perm, A::MatrixElem{T}) where {T <: Union{FieldElement, ResElem}}
   m = nrows(A)
   n = ncols(A)
   rank = 0
   r = 1
   c = 1
   R = base_ring(A)
   d = R(1)
   d2 = R(1)
   if m == 0 || n == 0
      return 0, d
   end
   t = R()
   while r <= m && c <= n
      if iszero(A[r, c])
         i = r + 1
         while i <= m
            if !iszero(A[i, c])
               for j = 1:n
                  A[i, j], A[r, j] = A[r, j], A[i, j]
               end
               P[r], P[i] = P[i], P[r]
               break
            end
            i += 1
         end
         if i > m
            c += 1
            continue
         end
      end
      rank += 1
      q = -A[r, c]
      for i = r + 1:m
         for j = c + 1:n
            A[i, j] = mul_red!(A[i, j], A[i, j], q, false)
            t = mul_red!(t, A[i, c], A[r, j], false)
            A[i, j] = addeq!(A[i, j], t)
            A[i, j] = reduce!(A[i, j])
            if r > 1
               A[i, j] = mul!(A[i, j], A[i, j], d)
            else
               A[i, j] = -A[i, j]
            end
         end
      end
      d = -inv(A[r, c])
      d2 = A[r, c]
      r += 1
      c += 1
   end
   inv!(P)
   return rank, d2
end

@doc Markdown.doc"""
    fflu(A::Generic.MatrixElem{T}, P = PermGroup(nrows(A))) where {T <: RingElement}
> Return a tuple $r, d, p, L, U$ consisting of the rank of $A$, a
> denominator $d$, a permutation $p$ of $A$ belonging to $P$, a lower
> triangular matrix $L$ and an upper triangular matrix $U$ such that
> $p(A) = LDU$, where $p(A)$ stands for the matrix whose rows are the given
> permutation $p$ of the rows of $A$ and such that $D$ is the diagonal matrix
> diag$(p_1, p_1p_2, \ldots, p_{n-2}p_{n-1}, p_{n-1})$ where the $p_i$ are the
> inverses of the diagonal entries of $U$. The denominator $d$ is set to
> $\pm \mbox{det}(S)$ where $S$ is an appropriate submatrix of $A$ ($S = A$ if
> $A$ is square) and the sign is decided by the parity of the permutation.
"""
function fflu(A::MatrixElem{T}, P = PermGroup(nrows(A))) where {T <: RingElement}
   m = nrows(A)
   n = ncols(A)
   P.n != m && error("Permutation does not match matrix")
   p = P()
   R = base_ring(A)
   U = deepcopy(A)
   L = similar(A, m, m)
   rank, d = fflu!(p, U)
   for i = 1:m
      for j = 1:n
         if i > j
            L[i, j] = U[i, j]
            U[i, j] = R()
         elseif i == j
            L[i, j] = U[i, j]
         elseif j <= m
            L[i, j] = R()
         end
      end
   end
   if m > 0
      L[m, m] = R(1)
   end
   return rank, d, p, L, U
end

###############################################################################
#
#   Reduced row-echelon form
#
###############################################################################

function rref!(A::MatrixElem{T}) where {T <: RingElement}
   m = nrows(A)
   n = ncols(A)
   R = base_ring(A)
   P = PermGroup(m)()
   rank, d = fflu!(P, A)
   for i = rank + 1:m
      for j = 1:n
         A[i, j] = R()
      end
   end
   if rank > 1
      t = R()
      q = R()
      d = -d
      pivots = zeros(Int, n)
      np = rank
      j = k = 1
      for i = 1:rank
         while iszero(A[i, j])
            pivots[np + k] = j
            j += 1
            k += 1
         end
         pivots[i] = j
         j += 1
      end
      while k <= n - rank
         pivots[np + k] = j
         j += 1
         k += 1
      end
      for k = 1:n - rank
         for i = rank - 1:-1:1
            t = mul_red!(t, A[i, pivots[np + k]], d, false)
            for j = i + 1:rank
               t = addmul_delayed_reduction!(t, A[i, pivots[j]], A[j, pivots[np + k]], q)
            end
            t = reduce!(t)
            A[i, pivots[np + k]] = divexact(-t, A[i, pivots[i]])
         end
      end
      d = -d
      for i = 1:rank
         for j = 1:rank
            if i == j
               A[j, pivots[i]] = d
            else
               A[j, pivots[i]] = R()
            end
         end
      end
   end
   return rank, d
end

@doc Markdown.doc"""
    rref(M::Generic.MatrixElem{T}) where {T <: RingElement}
> Return a tuple $(r, d, A)$ consisting of the rank $r$ of $M$ and a
> denominator $d$ in the base ring of $M$ and a matrix $A$ such that $A/d$ is
> the reduced row echelon form of $M$. Note that the denominator is not usually
> minimal.
"""
function rref(M::MatrixElem{T}) where {T <: RingElement}
   A = deepcopy(M)
   r, d = rref!(A)
   return r, d, A
end

function rref!(A::MatrixElem{T}) where {T <: FieldElement}
   m = nrows(A)
   n = ncols(A)
   R = base_ring(A)
   P = PermGroup(m)()
   rnk = lu!(P, A)
   if rnk == 0
      return 0
   end
   for i = 1:m
      for j = 1:min(rnk, i - 1)
         A[i, j] = R()
      end
   end
   U = zero_matrix(R, rnk, rnk)
   V = zero_matrix(R, rnk, n - rnk)
   pivots = zeros(Int, n)
   np = rnk
   j = k = 1
   for i = 1:rnk
      while iszero(A[i, j])
         pivots[np + k] = j
         j += 1
         k += 1
      end
      pivots[i] = j
      j += 1
   end
   while k <= n - rnk
      pivots[np + k] = j
      j += 1
      k += 1
   end
   for i = 1:rnk
      for j = 1:i
         U[j, i] = A[j, pivots[i]]
      end
      for j = i + 1:rnk
         U[j, i] = R()
      end
   end
   for i = 1:n - rnk
      for j = 1:rnk
         V[j, i] = A[j, pivots[np + i]]
      end
   end
   V = solve_triu(U, V, false)
   for i = 1:rnk
      for j = 1:i
         A[j, pivots[i]] = i == j ? R(1) : R()
      end
   end
   for i = 1:n - rnk
      for j = 1:rnk
         A[j, pivots[np + i]] = V[j, i]
      end
   end
   return rnk
end

@doc Markdown.doc"""
    rref(M::Generic.MatrixElem{T}) where {T <: FieldElement}
> Return a tuple $(r, A)$ consisting of the rank $r$ of $M$ and a reduced row
> echelon form $A$ of $M$.
"""
function rref(M::MatrixElem{T}) where {T <: FieldElement}
   A = deepcopy(M)
   r = rref!(A)
   return r, A
end

@doc Markdown.doc"""
    isrref(M::Generic.MatrixElem{T}) where {T <: RingElement}
> Return `true` if $M$ is in reduced row echelon form, otherwise return
> `false`.
"""
function isrref(M::MatrixElem{T}) where {T <: RingElement}
   m = nrows(M)
   n = ncols(M)
   c = 1
   for r = 1:m
      for i = 1:c - 1
         if !iszero(M[r, i])
            return false
         end
      end
      while c <= n && iszero(M[r, c])
         c += 1
      end
      if c <= n
         for i = 1:r - 1
            if !iszero(M[i, c])
               return false
            end
         end
      end
   end
   return true
end

@doc Markdown.doc"""
    isrref(M::Generic.MatrixElem{T}) where {T <: FieldElement}
> Return `true` if $M$ is in reduced row echelon form, otherwise return
> `false`.
"""
function isrref(M::MatrixElem{T}) where {T <: FieldElement}
   m = nrows(M)
   n = ncols(M)
   c = 1
   for r = 1:m
      for i = 1:c - 1
         if !iszero(M[r, i])
            return false
         end
      end
      while c <= n && iszero(M[r, c])
         c += 1
      end
      if c <= n
         if !isone(M[r, c])
            return false
         end
         for i = 1:r - 1
            if !iszero(M[i, c])
               return false
            end
         end
      end
   end
   return true
end

# Reduce row m of the matrix A, assuming the first m - 1 rows are in Gauss
# form. However those rows may not be in order. The i-th entry of the array
# P is the row of A which has a pivot in the i-th column. If no such row
# exists, the entry of P will be 0. The function returns the column in which
# the m-th row has a pivot after reduction. This will always be chosen to be
# the first available column for a pivot from the left. This information is
# also updated in P. The i-th entry of the array L contains the last column
# of A for which the row i is nonzero. This speeds up reduction in the case
# that A is chambered on the right. Otherwise the entries can all be set to
# the number of columns of A. The entries of L must be monotonic increasing.

function reduce_row!(A::MatrixElem{T}, P::Array{Int}, L::Array{Int}, m::Int) where {T <: FieldElement}
   R = base_ring(A)
   n = ncols(A)
   t = R()
   for i = 1:n
      # reduction of row was delayed, reduce next element now
      if i != 1
         A[m, i] = reduce!(A[m, i])
      end
      if !iszero(A[m, i])
         h = -A[m, i]
         r = P[i]
         if r != 0
            A[m, i] = R()
            for j = i + 1:L[r]
               t = mul_red!(t, A[r, j], h, false)
               A[m, j] = addeq!(A[m, j], t)
            end
         else
            # reduce remainder of row for return
            for j = i + 1:L[m]
               A[m, j] = reduce!(A[m, j])
            end
            h = inv(A[m, i])
            A[m, i] = R(1)
            for j = i + 1:L[m]
               A[m, j] = mul!(A[m, j], A[m, j], h)
            end
            P[i] = m
            return i
         end
      end
   end
   return 0
end

function reduce_row!(A::MatrixElem{T}, P::Array{Int}, L::Array{Int}, m::Int) where {T <: RingElement}
   R = base_ring(A)
   n = ncols(A)
   t = R()
   c = R(1)
   c1 = 0
   for i = 1:n
      if !iszero(A[m, i])
         h = -A[m, i]
         r = P[i]
         if r != 0
            d = A[r, i]
            A[m, i] = R()
            for j = i + 1:L[r]
               t = mul_red!(t, A[r, j], h, false)
               A[m, j] = mul_red!(A[m, j], A[m, j], d, false)
               A[m, j] = addeq!(A[m, j], t)
               A[m, j] = reduce!(A[m, j])
            end
            for j = L[r] + 1:L[m]
               A[m, j] = mul!(A[m, j], A[m, j], d)
            end
            if c1 > 0 && P[c1] < P[i]
               for j = i + 1:L[m]
                  A[m, j] = divexact(A[m, j], c)
               end
            end
            c1 = i
            c = d
         else
            P[i] = m
            return i
         end
      else
         r = P[i]
         if r != 0
            for j = i + 1:L[m]
               A[m, j] = mul!(A[m, j], A[m, j], A[r, i])
            end
         end
      end
   end
   return 0
end

###############################################################################
#
#   Determinant
#
###############################################################################

function det_clow(M::MatrixElem{T}) where {T <: RingElement}
   R = base_ring(M)
   n = nrows(M)
   if n == 0
      return one(R)
   end
   A = Array{T}(undef, n, n)
   B = Array{T}(undef, n, n)
   C = R()
   for i = 1:n
      for j = 1:n
         A[i, j] = i == j ? R(1) : R(0)
         B[i, j] = R()
      end
   end
   for k = 1:n - 1
      for i = 1:n
         for j = 1:i
            if !iszero(A[i, j])
               for m = j + 1:n
                  C = mul!(C, A[i, j], M[i, m])
                  B[m, j] = addeq!(B[m, j], C)
               end
               for m = j + 1:n
                  C = mul!(C, A[i, j], M[i, j])
                  B[m, m] = addeq!(B[m, m], -C)
               end
            end
         end
      end
      Temp = A
      A = B
      B = Temp
      if k != n - 1
         for i = 1:n
            for j = 1:i
               B[i, j] = R()
            end
         end
      end
   end
   D = R()
   for i = 1:n
      for j = 1:i
         if !iszero(A[i, j])
            D -= A[i, j]*M[i, j]
         end
      end
   end
   return isodd(n) ? -D : D
end

function det_df(M::MatrixElem{T}) where {T <: RingElement}
   R = base_ring(M)
   S, z = PolynomialRing(R, "z")
   n = nrows(M)
   p = charpoly(S, M)
   d = coeff(p, 0)
   return isodd(n) ? -d : d
end

function det_fflu(M::MatrixElem{T}) where {T <: RingElement}
   n = nrows(M)
   if n == 0
      return base_ring(M)()
   end
   A = deepcopy(M)
   P = PermGroup(n)()
   r, d = fflu!(P, A)
   return r < n ? base_ring(M)() : (parity(P) == 0 ? d : -d)
end

@doc Markdown.doc"""
    det(M::Generic.MatrixElem{T}) where {T <: FieldElement}
> Return the determinant of the matrix $M$. We assume $M$ is square.
"""
function det(M::MatrixElem{T}) where {T <: FieldElement}
   !issquare(M) && error("Not a square matrix in det")
   if nrows(M) == 0
      return one(base_ring(M))
   end
   return det_fflu(M)
end

@doc Markdown.doc"""
    det(M::Generic.MatrixElem{T}) where {T <: RingElement}
> Return the determinant of the matrix $M$. We assume $M$ is square.
"""
function det(M::MatrixElem{T}) where {T <: RingElement}
   !issquare(M) && error("Not a square matrix in det")
   if nrows(M) == 0
      return one(base_ring(M))
   end
   try
      return det_fflu(M)
   catch
      return det_df(M)
   end
end

function det_interpolation(M::MatrixElem{T}) where {T <: PolyElem}
   n = nrows(M)
   !isdomain_type(elem_type(typeof(base_ring(base_ring(M))))) &&
          error("Generic interpolation requires a domain type")
   R = base_ring(M)
   if n == 0
      return R()
   end
   maxlen = 0
   for i = 1:n
      for j = 1:n
         maxlen = max(maxlen, length(M[i, j]))
      end
   end
   if maxlen == 0
      return R()
   end
   bound = n*(maxlen - 1) + 1
   x = Array{elem_type(base_ring(R))}(undef, bound)
   d = Array{elem_type(base_ring(R))}(undef, bound)
   X = zero_matrix(base_ring(R), n, n)
   b2 = AbstractAlgebra.div(bound, 2)
   pt1 = base_ring(R)(1 - b2)
   for i = 1:bound
      x[i] = base_ring(R)(i - b2)
      (x[i] == pt1 && i != 1) && error("Not enough interpolation points in ring")
      for j = 1:n
         for k = 1:n
            X[j, k] = evaluate(M[j, k], x[i])
         end
      end
      d[i] = det(X)
   end
   return interpolate(R, x, d)
end

function det(M::MatrixElem{T}) where {S <: FinFieldElem, T <: PolyElem{S}}
   !issquare(M) && error("Not a square matrix in det")
   if nrows(M) == 0
      return one(base_ring(M))
   end
   return det_popov(M)
end

function det(M::MatrixElem{T}) where {T <: PolyElem}
   !issquare(M) && error("Not a square matrix in det")
   if nrows(M) == 0
      return one(base_ring(M))
   end
   try
      return det_interpolation(M)
   catch
      # no point trying fflu, since it probably fails
      # for same reason as det_interpolation
      return det_df(M)
   end
end

###############################################################################
#
#   Minors
#
###############################################################################

@doc Markdown.doc"""
    combinations(n::Int,k::Int)
> Return an array consisting of k-combinations of {1,...,n} as arrays.
"""
function combinations(n, k)
   if (k == 0)
      r = Array{Array{Int, 1}, 1}(undef, 1)
      r[1] = []
      return r
   elseif k > n
      return Array{Array{Int, 1 }, 1}(undef, 0)
   elseif k == n
      return [collect(1:k)]
   else
      return vcat(combinations(n - 1, k), [append!(l, [n]) for l in combinations(n - 1, k - 1)])
   end
end

@doc Markdown.doc"""
    minors(A::AbstractAlgebra.MatElem, k::Int)
> Return an array consisting of the k-minors of A
"""
function minors(A::AbstractAlgebra.MatElem, k::Int)
   row_indices = combinations(nrows(A), k)
   col_indices = combinations(ncols(A), k)
   mins = Array{elem_type(base_ring(A)), 1}(undef, 0)
   for ri in row_indices
      for ci in col_indices
         push!(mins, det(AbstractAlgebra.Generic.sub(A, ri, ci)))
      end
   end
   return(mins)
end

###############################################################################
#
#   Rank
#
###############################################################################

@doc Markdown.doc"""
    rank(M::Generic.MatrixElem{T}) where {T <: RingElement}
> Return the rank of the matrix $M$.
"""
function rank(M::MatrixElem{T}) where {T <: RingElement}
   n = nrows(M)
   if n == 0
      return 0
   end
   A = deepcopy(M)
   P = PermGroup(n)()
   r, d = fflu!(P, A)
   return r
end

@doc Markdown.doc"""
    rank(M::Generic.MatrixElem{T}) where {T <: FieldElement}
> Return the rank of the matrix $M$.
"""
function rank(M::MatrixElem{T}) where {T <: FieldElement}
   n = nrows(M)
   if n == 0
      return 0
   end
   A = deepcopy(M)
   P = PermGroup(n)()
   return lu!(P, A)
end

###############################################################################
#
#   Linear solving
#
###############################################################################

function solve_fflu(A::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   base_ring(A) != base_ring(b) && error("Base rings don't match in solve_fflu")
   !issquare(A) && error("Non-square matrix in solve_fflu")
   nrows(A) != nrows(b) && error("Dimensions don't match in solve_fflu")
   FFLU = deepcopy(A)
   p = PermGroup(nrows(A))()
   r, d = fflu!(p, FFLU)
   r < nrows(A) && error("Singular matrix in solve_fflu")
   return solve_fflu_precomp(p, FFLU, b), d
end

function solve_fflu_precomp(p::Generic.perm, FFLU::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   x = p * b
   n = nrows(x)
   m = ncols(x)
   R = base_ring(FFLU)

   t = base_ring(b)()
   s = base_ring(b)()
   minus_one = R(-1)

   for k in 1:m
      for i in 1:(n - 1)
         t = mul!(t, x[i, k], minus_one)
         for j in (i + 1):n
            if i == 1
              x[j, k] = mul_red!(R(), x[j, k], FFLU[i, i], false)
            else
              x[j, k] = mul_red!(x[j, k], x[j, k], FFLU[i, i], false)
            end
            s = mul_red!(s, FFLU[j, i], t, false)
            x[j, k] = addeq!(x[j, k], s)
            x[j, k] = reduce!(x[j, k])
            if i > 1
                x[j, k] = divexact(x[j, k], FFLU[i - 1, i - 1])
            end
         end
      end

      for i in (n - 1):-1:1
         if i > 1
            x[i, k] = mul!(x[i, k], x[i, k], FFLU[n, n])
         else
            x[i, k] = x[i, k] * FFLU[n, n]
         end
         for j in (i + 1):n
            t = mul!(t, x[j, k], FFLU[i, j])
            t = mul!(t, t, minus_one)
            x[i, k] = addeq!(x[i, k], t)
         end
         x[i, k] = divexact(x[i, k], FFLU[i, i])
      end
   end
   return x
end

function solve_lu(A::MatElem{T}, b::MatElem{T}) where {T <: FieldElement}
   base_ring(A) != base_ring(b) && error("Base rings don't match in solve_lu")
   !issquare(A) && error("Non-square matrix in solve_lu")
   nrows(A) != nrows(b) && error("Dimensions don't match in solve_lu")

   if nrows(A) == 0 || ncols(A) == 0
      return b
   end

   LU = deepcopy(A)
   p = PermGroup(nrows(A))()
   r = lu!(p, LU)
   r < nrows(A) && error("Singular matrix in solve_lu")
   return solve_lu_precomp(p, LU, b)
end

function solve_lu_precomp(p::Generic.perm, LU::MatElem{T}, b::MatrixElem{T}) where {T <: FieldElement}
   x = p * b
   n = nrows(x)
   m = ncols(x)
   R = base_ring(LU)

   t = base_ring(b)()
   s = base_ring(b)()

   for k in 1:m
      x[1, k] = deepcopy(x[1, k])
      for i in 2:n
         for j in 1:(i - 1)
            # x[i, k] = x[i, k] - LU[i, j] * x[j, k]
            t = mul_red!(t, -LU[i, j], x[j, k], false)
            if j == 1
               x[i, k] = x[i, k] + t # LU[i, j] * x[j, k]
            else
               x[i, k] = addeq!(x[i, k], t)
            end
         end
         x[i, k] = reduce!(x[i, k])
      end

      # Now every entry of x is a proper copy, so we can change the entries
      # as much as we want.

      x[n, k] = divexact(x[n, k], LU[n, n])

      for i in (n - 1):-1:1
         for j in (i + 1):n
            # x[i, k] = x[i, k] - x[j, k] * LU[i, j]
            t = mul_red!(t, x[j, k], -LU[i, j], false)
            x[i, k] = addeq!(x[i, k], t)
         end
         x[i, k] = reduce!(x[i, k])
         x[i, k] = divexact(x[i, k], LU[i, i])
      end
   end
   return x
end

function solve_ff(M::MatrixElem{T}, b::MatrixElem{T}) where {T <: FieldElement}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   !issquare(M) && error("Non-square matrix in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   m = nrows(M)
   x, d = solve_fflu(M, b)
   for i in 1:nrows(x)
      for j in 1:ncols(x)
         x[i, j] = divexact(x[i, j], d)
      end
   end
   return x
end

function solve_with_det(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: RingElement}
   # We cannot use solve_fflu directly, since it forgot about the (parity of
   # the) permutation.
   nrows(M) != ncols(M) && error("Non-square matrix")
   R = base_ring(M)
   FFLU = deepcopy(M)
   p = PermGroup(nrows(M))()
   r, d = fflu!(p, FFLU)
   if r < nrows(M)
      error("Singular matrix in solve_with_det")
   end
   x = solve_fflu_precomp(p, FFLU, b)
   # Now M*x = d*b, but d is only sign(P) * det(M)
   if parity(p) != 0
      minus_one = R(-1)
      for k in 1:ncols(x)
         for i in 1:nrows(x)
            # We are allowed to modify x in-place.
            x[i, k] = mul!(x[i, k], x[i, k], minus_one)
         end
      end
      d = mul!(d, d, minus_one)
   end
   return x, d
end

function solve_with_det(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: PolyElem}
   x, d = solve_interpolation(M, b)
   return x, d
end

function solve_ff(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: RingElement}
   m = nrows(M)
   n = ncols(M)
   if m == 0 || n == 0
      return b, base_ring(M)()
   end
   return solve_fflu(M, b)
end

function solve_interpolation(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: PolyElem}
   m = nrows(M)
   h = ncols(b)
   if m == 0
      return b, base_ring(M)()
   end
   R = base_ring(M)
   maxlen = 0
   for i = 1:m
      for j = 1:m
         maxlen = max(maxlen, length(M[i, j]))
      end
   end
   maxlen == 0 && error("Matrix is singular in solve")
   maxlenb = 0
   for i = 1:m
      for j = 1:h
         maxlenb = max(maxlenb, length(b[i, j]))
      end
   end
   # bound from xd = (M*)b where d is the det
   bound = (maxlen - 1)*(m - 1) + max(maxlenb, maxlen)
   tmat = matrix(base_ring(R), 0, 0, elem_type(base_ring(R))[])
   V = Array{typeof(tmat)}(undef, bound)
   d = Array{elem_type(base_ring(R))}(undef, bound)
   y = Array{elem_type(base_ring(R))}(undef, bound)
   bj = Array{elem_type(base_ring(R))}(undef, bound)
   X = similar(tmat, m, m)
   Y = similar(tmat, m, h)
   x = similar(b)
   b2 = AbstractAlgebra.div(bound, 2)
   pt1 = base_ring(R)(1 - b2)
   l = 1
   i = 1
   pt = 1
   while l <= bound
      y[l] = base_ring(R)(pt - b2)
      (y[l] == pt1 && pt != 1) && error("Not enough interpolation points in ring")
      for j = 1:m
         for k = 1:m
            X[j, k] = evaluate(M[j, k], y[l])
         end
         for k = 1:h
            Y[j, k] = evaluate(b[j, k], y[l])
         end
      end
      try
         V[l], d[l] = solve_with_det(X, Y)
         l = l + 1
      catch e
         if !(e isa ErrorException)
            rethrow(e)
         end
         i = i + 1
      end

      # We tested bound evaluation points and either an impossible
      # inverse was encountered, or the matrix was singular for all
      # the values.

      if i > bound && l == 1
         error("Impossible inverse or too many singular matrices in solve_interpolation")
      end

      pt = pt + 1
   end
   for k = 1:h
      for i = 1:m
         for j = 1:bound
            bj[j] = V[j][i, k]
         end
         x[i, k] = interpolate(R, y, bj)
      end
   end
   return x, interpolate(R, y, d)
end

@doc Markdown.doc"""
    solve(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: FieldElement}
> Given a non-singular $n\times n$ matrix over a field and an $n\times m$
> matrix over the same field, return $x$ an
> $n\times m$ matrix $x$ such that $Ax = b$.
> If $A$ is singular an exception is raised.
"""
function solve(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: FieldElement}
    return solve_ringelem(M, b)
end

@doc Markdown.doc"""
    solve_rational(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where T <: RingElement
> Given a non-singular $n\times n$ matrix over a ring and an $n\times m$
> matrix over the same ring, return a tuple $x, d$ consisting of an
> $n\times m$ matrix $x$ and a denominator $d$ such that $Ax = db$. The
> denominator will be the determinant of $A$ up to sign. If $A$ is singular an
> exception is raised.
"""
function solve_rational(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where T <: RingElement
   return solve_ringelem(M, b)
end

function solve_ringelem(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: RingElement}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   nrows(M) != ncols(M) && error("Non-square matrix in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   return solve_ff(M, b)
end

function solve_rational(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: PolyElem}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   nrows(M) != ncols(M) && error("Non-square matrix in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   try
      return solve_interpolation(M, b)
   catch e
      if !isa(e, ErrorException)
         rethrow(e)
      end
      return solve_ff(M, b)
   end
end

@doc Markdown.doc"""
    solve_left(a::AbstractAlgebra.MatElem{S}, b::AbstractAlgebra.MatElem{S}) where S <: RingElement
> Given an $r\times n$ matrix $a$ over a ring and an $m\times n$ matrix $b$
> over the same ring, return an $m\times r$ matrix $x$ such that $xa = b$. If
> no such matrix exists, an exception is raised.
"""
function solve_left(a::AbstractAlgebra.MatElem{S}, b::AbstractAlgebra.MatElem{S}) where S <: RingElement
  @assert ncols(a) == ncols(b)
  H, T = hnf_with_transform(a)
  b = deepcopy(b)
  z = similar(a, nrows(b), nrows(a))
  l = min(ncols(a), nrows(a))
  t = base_ring(a)()
  for i = 1:nrows(b)
    for j = 1:l
      k = 1
      while k <= ncols(H) && iszero(H[j, k])
        k += 1
      end
      if k > ncols(H)
        continue
      end
      q, r = AbstractAlgebra.divrem(b[i, k], H[j, k])
      r != 0 && error("Unable to solve linear system")
      z[i, j] = q
      q = -q
      for h = k:ncols(H)
        t = mul!(t, q, H[j, h])
        b[i, h] = addeq!(b[i, h], t)
      end
    end
  end
  b != 0 && error("Unable to solve linear system")
  return z*T
end

# Find the pivot columns of an rref matrix
function find_pivot(A::AbstractAlgebra.MatElem{T}) where T <: RingElement
  p = Int[]
  j = 0
  for i = 1:nrows(A)
    j += 1
    if j > ncols(A)
      return p
    end
    while iszero(A[i, j])
      j += 1
      if j > ncols(A)
        return p
      end
    end
    push!(p, j)
  end
  return p
end

function solve_left(A::AbstractAlgebra.MatElem{T}, B::AbstractAlgebra.MatElem{T}) where T <: FieldElement
  R = base_ring(A)
  ncols(A) != ncols(B) && error("Incompatible matrices")
  mu = zero_matrix(R, ncols(A), nrows(A) + nrows(B))
  for i = 1:ncols(A)
     for j = 1:nrows(A)
        mu[i, j] = A[j, i]
     end
     for j = 1:nrows(B)
        mu[i, nrows(A) + j] = B[j, i]
     end
  end
  rk, mu = rref(mu)
  p = find_pivot(mu)
  if any(i -> i > nrows(A), p)
    error("Unable to solve linear system")
  end
  sol = zero_matrix(R, nrows(B), nrows(A))
  for i = 1:length(p)
    for j = 1:nrows(B)
      sol[j, p[i]] = mu[i, nrows(A) + j]
    end
  end
  return sol
end

###############################################################################
#
#   Upper triangular solving
#
###############################################################################

@doc Markdown.doc"""
    solve_triu(U::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}, unit::Bool = false) where {T <: FieldElement}
> Given a non-singular $n\times n$ matrix over a field which is upper
> triangular, and an $n\times m$ matrix over the same field, return an
> $n\times m$ matrix $x$ such that $Ax = b$. If $A$ is singular an exception
> is raised. If unit is true then $U$ is assumed to have ones on its
> diagonal, and the diagonal will not be read.
"""
function solve_triu(U::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}, unit::Bool = false) where {T <: FieldElement}
   n = nrows(U)
   m = ncols(b)
   R = base_ring(U)
   X = similar(b)
   Tinv = Array{elem_type(R)}(undef, n)
   tmp = Array{elem_type(R)}(undef, n)
   if unit == false
      for i = 1:n
         Tinv[i] = inv(U[i, i])
      end
   end
   t = R()
   for i = 1:m
      for j = 1:n
         tmp[j] = X[j, i]
      end
      for j = n:-1:1
         s = R()
         for k = j + 1:n
            s = addmul_delayed_reduction!(s, U[j, k], tmp[k], t)
         end
         s = reduce!(s)
         s = b[j, i] - s
         if unit == false
            s = mul!(s, s, Tinv[j])
         end
         tmp[j] = s
      end
      for j = 1:n
         X[j, i] = tmp[j]
      end
   end
   return X
end

###############################################################################
#
#   Can solve
#
###############################################################################

@doc Markdown.doc"""
    can_solve_left_reduced_triu(r::AbstractAlgebra.MatElem{T},
                          M::AbstractAlgebra.MatElem{T}) where T <: RingElement
> Return a tuple `flag, x` where `flag` is set to true if $xM = r$ has a
> solution, where $M$ is an $m\times n$ matrix in (upper triangular) Hermite
> normal form or reduced row echelon form and $r$ and $x$ are row vectors with
> $m$ columns. If there is no solution, flag is set to `false` and $x$ is set
> to the zero row.
"""
function can_solve_left_reduced_triu(r::AbstractAlgebra.MatElem{T},
                          M::AbstractAlgebra.MatElem{T}) where T <: RingElement
   ncols(r) != ncols(M) && error("Incompatible matrices")
   r = deepcopy(r) # do not destroy input
   m = ncols(r)
   n = nrows(M)
   if n == 0
      return true, r
   end
   R = base_ring(r)
   x = zero_matrix(R, 1, n)
   j = 1 # row in M
   k = 1 # column in M
   t = R()
   for i = 1:m # column in r
      if iszero(r[1, i])
         continue
      end
      while k <= i && j <= n
         if iszero(M[j, k])
            k += 1
         elseif k < i
            j += 1
         else
            break
         end
      end
      if k != i
         return false, x
      end
      x[1, j], r[1, i] = AbstractAlgebra.divrem(r[1, i], M[j, k])
      if !iszero(r[1, i])
         return false, x
      end
      q = -x[1, j]
      for l = i + 1:m
         t = mul!(t, q, M[j, l])
         r[1, l] = addeq!(r[1, l], t)
      end
   end
   return true, x
end

###############################################################################
#
#   Inverse
#
###############################################################################

@doc Markdown.doc"""
    pseudo_inv(M::Generic.MatrixElem{T}) where {T <: RingElement}
> Given a non-singular $n\times n$ matrix $M$ over a ring return a tuple $X, d$
> consisting of an $n\times n$ matrix $X$ and a denominator $d$ such that
> $MX = dI_n$, where $I_n$ is the $n\times n$ identity matrix. The denominator
> will be the determinant of $M$ up to sign. If $M$ is singular an exception
> is raised.
"""
function pseudo_inv(M::MatrixElem{T}) where {T <: RingElement}
   issquare(M) || throw(DomainError(M, "Can not invert non-square Matrix"))
   X, d = solve_fflu(M, eye(M))
   return X, d
end

@doc Markdown.doc"""
    inv(M::Generic.MatrixElem{T}) where {T <: FieldElement}
> Given a non-singular $n\times n$ matrix over a field, return an
> $n\times n$ matrix $X$ such that $MX = I_n$ where $I_n$ is the $n\times n$
> identity matrix. If $M$ is singular an exception is raised.
"""
function inv(M::MatrixElem{T}) where {T <: FieldElement}
   issquare(M) || throw(DomainError(M, "Can not invert non-square Matrix"))
   A = solve_lu(M, eye(M))
   return A
end

@doc Markdown.doc"""
    inv(M::Generic.MatrixElem{T}) where {T <: RingElement}
> Given a non-singular $n\times n$ matrix over a ring, return an
> $n\times n$ matrix $X$ such that $MX = I_n$, where $I_n$ is the $n\times n$
> identity matrix. If $M$ is not invertible over the base ring an exception is
> raised.
"""
function inv(M::MatrixElem{T}) where {T <: RingElement}
   issquare(M) || throw(DomainError(M, "Can not invert non-square Matrix"))
   X, d = pseudo_inv(M)
   isunit(d) || throw(DomainError(M, "Matrix is not invertible."))
   return divexact(X, d)
end

###############################################################################
#
#   Nullspace
#
###############################################################################

@doc Markdown.doc"""
    nullspace(M::AbstractAlgebra.MatElem{T}) where {T <: RingElement}
> Return a tuple $(\nu, N)$ consisting of the nullity $\nu$ of $M$ and
> a basis $N$ (consisting of column vectors) for the right nullspace of $M$,
> i.e. such that $MN$ is the zero matrix. If $M$ is an $m\times n$ matrix
> $N$ will be an $n\times \nu$ matrix. Note that the nullspace is taken to be
> the vector space kernel over the fraction field of the base ring if the
> latter is not a field. In AbstractAlgebra we use the name ``kernel'' for a
> function to compute an integral kernel.
"""
function nullspace(M::AbstractAlgebra.MatElem{T}) where {T <: RingElement}
   n = ncols(M)
   rank, d, A = rref(M)
   nullity = n - rank
   R = base_ring(M)
   U = similar(M, n, nullity)
   if rank == 0
      for i = 1:nullity
         U[i, i] = R(1)
      end
   elseif nullity != 0
      pivots = zeros(Int, rank)
      nonpivots = zeros(Int, nullity)
      j = k = 1
      for i = 1:rank
         while iszero(A[i, j])
            nonpivots[k] = j
            j += 1
            k += 1
         end
         pivots[i] = j
         j += 1
      end
      while k <= nullity
         nonpivots[k] = j
         j += 1
         k += 1
      end
      d = -A[1, pivots[1]]
      for i = 1:nullity
         for j = 1:rank
            U[pivots[j], i] = A[j, nonpivots[i]]
         end
         U[nonpivots[i], i] = d
      end
   end
   return nullity, U
end

@doc Markdown.doc"""
    nullspace(M::AbstractAlgebra.MatElem{T}) where {T <: FieldElement}
> Return a tuple $(\nu, N)$ consisting of the nullity $\nu$ of $M$ and
> a basis $N$ (consisting of column vectors) for the right nullspace of $M$,
> i.e. such that $MN$ is the zero matrix. If $M$ is an $m\times n$ matrix
> $N$ will be an $n\times \nu$ matrix.
"""
function nullspace(M::AbstractAlgebra.MatElem{T}) where {T <: FieldElement}
   m = nrows(M)
   n = ncols(M)
   rank, A = rref(M)
   nullity = n - rank
   R = base_ring(M)
   X = similar(M, n, nullity)
   if rank == 0
      for i = 1:nullity
         X[i, i] = R(1)
      end
   elseif nullity != 0
      pivots = zeros(Int, max(m, n))
      np = rank
      j = k = 1
      for i = 1:rank
         while iszero(A[i, j])
            pivots[np + k] = j
            j += 1
            k += 1
         end
         pivots[i] = j
         j += 1
      end
      while k <= nullity
         pivots[np + k] = j
         j += 1
         k += 1
      end
      for i = 1:nullity
         for j = 1:rank
            X[pivots[j], i] = -A[j, pivots[np + i]]
         end
         X[pivots[np + i], i] = R(1)
      end
   end
   return nullity, X
end

###############################################################################
#
#   Kernel
#
###############################################################################

@doc Markdown.doc"""
    left_kernel(a::AbstractAlgebra.MatElem{T}) where T <: RingElement
> Return a tuple `n, M` where $M$ is a matrix whose rows generate the kernel
> of $M$ and $n$ is the rank of the kernel. The transpose of the output of this
> function is guaranteed to be in flipped upper triangular format (i.e. upper
> triangular format if columns and rows are reversed).
"""
function left_kernel(x::AbstractAlgebra.MatElem{T}) where T <: RingElement
   !isdomain_type(elem_type(base_ring(x))) && error("Not implemented")
   R = base_ring(x)
   H, U = hnf_with_transform(x)
   i = nrows(H)
   zero_rows = false
   while i > 0 && iszero_row(H, i)
      zero_rows = true
      i -= 1
   end
   if zero_rows
      return nrows(U) - i, U[i + 1:nrows(U), 1:ncols(U)]
   else
      return 0, zero_matrix(R, 0, ncols(U))
   end
end

function left_kernel(M::AbstractAlgebra.MatElem{T}) where T <: FieldElement
  n, N = nullspace(transpose(M))
  return n, transpose(N)
end

@doc Markdown.doc"""
    right_kernel(a::AbstractAlgebra.MatElem{T}) where T <: RingElement
> Return a tuple `n, M` where $M$ is a matrix whose columns generate the
> kernel of $M$ and $n$ is the rank of the kernel.
"""
function right_kernel(x::AbstractAlgebra.MatElem{T}) where T <: RingElement
   n, M = left_kernel(transpose(x))
   return n, transpose(M)
end

function right_kernel(M::AbstractAlgebra.MatElem{T}) where T <: FieldElement
   return nullspace(M)
end

@doc Markdown.doc"""
    kernel(a::MatElem{T}; side::Symbol = :right) where T <: RingElement
> Return a tuple $(n, M)$, where n is the rank of the kernel and $M$ is a
> basis for it. If side is $:right$ or not specified, the right kernel is
> computed, i.e. the matrix of columns whose span gives the right kernel
> space. If side is $:left$, the left kernel is computed, i.e. the matrix
> of rows whose span is the left kernel space.
"""
function kernel(A::AbstractAlgebra.MatElem{T}; side::Symbol = :right) where T <: RingElement
   if side == :right
      return right_kernel(A)
   elseif side == :left
      return left_kernel(A)
   else
      error("Unsupported argument: :$side for side: Must be :left or :right")
   end
end

###############################################################################
#
#   Hessenberg form
#
###############################################################################

function hessenberg!(A::MatrixElem{T}) where {T <: RingElement}
   !issquare(A) && error("Dimensions don't match in hessenberg")
   R = base_ring(A)
   n = nrows(A)
   u = R()
   t = R()
   for m = 2:n - 1
      i = m + 1
      while i <= n && iszero(A[i, m - 1])
         i += 1
      end
      if i != n + 1
         if !iszero(A[m, m - 1])
            i = m
         end
         h = -inv(A[i, m - 1])
         if i > m
            for j = m - 1:n
               A[i, j], A[m, j] = A[m, j], A[i, j]
            end
            for j = 1:n
               A[j, i], A[j, m] = A[j, m], A[j, i]
            end
         end
         for i = m + 1:n
            if !iszero(A[i, m - 1])
               u = mul!(u, A[i, m - 1], h)
               for j = m:n
                  t = mul!(t, u, A[m, j])
                  A[i, j] = addeq!(A[i, j], t)
               end
               u = -u
               for j = 1:n
                  t = mul!(t, u, A[j, i])
                  A[j, m] = addeq!(A[j, m], t)
               end
               A[i, m - 1] = R()
            end
         end
      end
   end
end

@doc Markdown.doc"""
    hessenberg(A::Generic.MatrixElem{T}) where {T <: RingElement}
> Return the Hessenberg form of $M$, i.e. an upper Hessenberg matrix
> which is similar to $M$. The upper Hessenberg form has nonzero entries
> above and on the diagonal and in the diagonal line immediately below the
> diagonal.
"""
function hessenberg(A::MatrixElem{T}) where {T <: RingElement}
   !issquare(A) && error("Dimensions don't match in hessenberg")
   M = deepcopy(A)
   hessenberg!(M)
   return M
end

@doc Markdown.doc"""
    ishessenberg(A::Generic.MatrixElem{T}) where {T <: RingElement}
> Return `true` if $M$ is in Hessenberg form, otherwise returns `false`.
"""
function ishessenberg(A::MatrixElem{T}) where {T <: RingElement}
   if !issquare(A)
      return false
   end
   n = nrows(A)
   for i = 3:n
      for j = 1:i - 2
         if !iszero(A[i, j])
            return false
         end
      end
   end
   return true
end

###############################################################################
#
#   Characteristic polynomial
#
###############################################################################

function charpoly_hessenberg!(S::Ring, A::MatrixElem{T}) where {T <: RingElement}
   !issquare(A) && error("Dimensions don't match in charpoly")
   R = base_ring(A)
   base_ring(S) != base_ring(A) && error("Cannot coerce into polynomial ring")
   n = nrows(A)
   if n == 0
      return S(1)
   end
   if n == 1
      return gen(S) - A[1, 1]
   end
   hessenberg!(A)
   P = Array{elem_type(S)}(undef, n + 1)
   P[1] = S(1)
   x = gen(S)
   for m = 1:n
      P[m + 1] = (x - A[m, m])*P[m]
      t = R(1)
      for i = 1:m - 1
         t = mul!(t, t, A[m - i + 1, m - i])
         P[m + 1] -= t*A[m - i, m]*P[m - i]
      end
   end
   return P[n + 1]
end

function charpoly_danilevsky_ff!(S::Ring, A::MatrixElem{T}) where {T <: RingElement}
   !issquare(A) && error("Dimensions don't match in charpoly")
   R = base_ring(A)
   base_ring(S) != base_ring(A) && error("Cannot coerce into polynomial ring")
   n = nrows(A)
   if n == 0
      return S(1)
   end
   if n == 1
      return gen(S) - A[1, 1]
   end
   d = R(1)
   t = R()
   V = Array{T}(undef, n)
   W = Array{T}(undef, n)
   pol = S(1)
   i = 1
   while i < n
      h = A[n - i + 1, n - i]
      while iszero(h)
         k = 1
         while k < n - i && iszero(A[n - i + 1, n - i - k])
            k += 1
         end
         if k == n - i
            b = S()
            fit!(b, i + 1)
            b = setcoeff!(b, i, R(1))
            for kk = 1:i
               b = setcoeff!(b, kk - 1, -A[n - i + 1, n - kk + 1]*d)
            end
            pol *= b
            n -= i
            i = 1
            if n == 1
               pol *= (gen(S) - A[1, 1]*d)
               return pol
            end
         else
            for j = 1:n
               A[n - i - k, j], A[n - i, j] = A[n - i, j], A[n - i - k, j]
            end
            for j = 1:n - i + 1
               A[j, n - i - k], A[j, n - i] = A[j, n - i], A[j, n - i - k]
            end
         end
         h = A[n - i + 1, n - i]
      end
      for j = 1:n
         V[j] = -A[n - i + 1, j]
         W[j] = deepcopy(A[n - i + 1, j])
      end
      for j = 1:n - i
         for kk = 1:n - i - 1
            t = mul_red!(t, A[j, n - i], V[kk], false)
            A[j, kk] = mul_red!(A[j, kk], A[j, kk], h, false)
            A[j, kk] = addeq!(A[j, kk], t)
            A[j, kk] = reduce!(A[j, kk])
         end
         for kk = n - i + 1:n
            t = mul_red!(t, A[j, n - i], V[kk], false)
            A[j, kk] = mul_red!(A[j, kk], A[j, kk], h, false)
            A[j, kk] = addeq!(A[j, kk], t)
            A[j, kk] = reduce!(A[j, kk])
         end
      end
      for kk = 1:n
         A[n - i + 1, kk] = R()
      end
      for j = 1:n - i
         for kk = 1:n - i - 1
            A[j, kk] = mul!(A[j, kk], A[j, kk], d)
         end
         for kk = n - i + 1:n
            A[j, kk] = mul!(A[j, kk], A[j, kk], d)
         end
      end
      A[n - i + 1, n - i] = deepcopy(h)
      for j = 1:n - i - 1
         s = R()
         for kk = 1:n - i
            s = addmul_delayed_reduction!(s, A[kk, j], W[kk], t)
         end
         s = reduce!(s)
         A[n - i, j] = s
      end
      for j = n - i:n - 1
         s = R()
         for kk = 1:n - i
            s = addmul_delayed_reduction!(s, A[kk, j], W[kk], t)
         end
         s = addmul_delayed_reduction!(s, h, W[j + 1], t)
         s = reduce!(s)
         A[n - i, j] = s
      end
      s = R()
      for kk = 1:n - i
         s = addmul_delayed_reduction!(s, A[kk, n], W[kk], t)
      end
      s = reduce!(s)
      A[n - i, n] = s
      for kk = 1:n
         A[n - i, kk] = mul!(A[n - i, kk], A[n - i, kk], d)
      end
      d = inv(h)
      parent(d) # To work around a bug in julia
      i += 1
   end
   b = S()
   fit!(b, n + 1)
   b = setcoeff!(b, n, R(1))
   for i = 1:n
      c = -A[1, n - i + 1]*d
      b = setcoeff!(b, i - 1, c)
   end
   return pol*b
end

function charpoly_danilevsky!(S::Ring, A::MatrixElem{T}) where {T <: RingElement}
   !issquare(A) && error("Dimensions don't match in charpoly")
   R = base_ring(A)
   base_ring(S) != base_ring(A) && error("Cannot coerce into polynomial ring")
   n = nrows(A)
   if n == 0
      return S(1)
   end
   if n == 1
      return gen(S) - A[1, 1]
   end
   t = R()
   V = Array{T}(undef, n)
   W = Array{T}(undef, n)
   pol = S(1)
   i = 1
   while i < n
      h = A[n - i + 1, n - i]
      while iszero(h)
         k = 1
         while k < n - i && iszero(A[n - i + 1, n - i - k])
            k += 1
         end
         if k == n - i
            b = S()
            fit!(b, i + 1)
            b = setcoeff!(b, i, R(1))
            for kk = 1:i
               b = setcoeff!(b, kk - 1, -A[n - i + 1, n - kk + 1])
            end
            pol *= b
            n -= i
            i = 1
            if n == 1
               pol *= (gen(S) - A[1, 1])
               return pol
            end
         else
            for j = 1:n
               A[n - i - k, j], A[n - i, j] = A[n - i, j], A[n - i - k, j]
            end
            for j = 1:n - i + 1
               A[j, n - i - k], A[j, n - i] = A[j, n - i], A[j, n - i - k]
            end
         end
         h = A[n - i + 1, n - i]
      end
      h = -inv(h)
      for j = 1:n
         V[j] = A[n - i + 1, j]*h
         W[j] = deepcopy(A[n - i + 1, j])
      end
      h = -h
      for j = 1:n - i
         for k = 1:n - i - 1
            t = mul!(t, A[j, n - i], V[k])
            A[j, k] = addeq!(A[j, k], t)
         end
         for k = n - i + 1:n
            t = mul!(t, A[j, n - i], V[k])
            A[j, k] = addeq!(A[j, k], t)
         end
         A[j, n - i] = mul!(A[j, n - i], A[j, n - i], h)
      end
      for j = 1:n - i - 1
         s = R()
         for k = 1:n - i
            s = addmul_delayed_reduction!(s, A[k, j], W[k], t)
         end
         s = reduce!(s)
         A[n - i, j] = s
      end
      for j = n - i:n - 1
         s = R()
         for k = 1:n - i
            s = addmul_delayed_reduction!(s, A[k, j], W[k], t)
         end
         s = addeq!(s, W[j + 1])
         s = reduce!(s)
         A[n - i, j] = s
      end
      s = R()
      for k = 1:n - i
         s = addmul_delayed_reduction!(s, A[k, n], W[k], t)
      end
      s = reduce!(s)
      A[n - i, n] = s
      i += 1
   end
   b = S()
   fit!(b, n + 1)
   b = setcoeff!(b, n, R(1))
   for i = 1:n
      b = setcoeff!(b, i - 1, -A[1, n - i + 1])
   end
   return pol*b
end

@doc Markdown.doc"""
    charpoly(V::Ring, Y::Generic.MatrixElem{T}) where {T <: RingElement}
> Return the characteristic polynomial $p$ of the matrix $M$. The
> polynomial ring $R$ of the resulting polynomial must be supplied
> and the matrix is assumed to be square.
"""
function charpoly(V::Ring, Y::MatrixElem{T}) where {T <: RingElement}
   !issquare(Y) && error("Dimensions don't match in charpoly")
   R = base_ring(Y)
   base_ring(V) != base_ring(Y) && error("Cannot coerce into polynomial ring")
   n = nrows(Y)
   if n == 0
      return V(1)
   end
   F = Array{elem_type(R)}(undef, n)
   A = Array{elem_type(R)}(undef, n)
   M = Array{elem_type(R)}(undef, n - 1, n)
   F[1] = -Y[1, 1]
   for i = 2:n
      F[i] = R()
      for j = 1:i
         M[1, j] = Y[j, i]
      end
      A[1] = Y[i, i]
      p = R()
      for j = 2:i - 1
         for k = 1:i
            s = R()
            for l = 1:i
               s = addmul_delayed_reduction!(s, Y[k, l], M[j - 1, l], p)
            end
            s = reduce!(s)
            M[j, k] = s
         end
         A[j] = M[j, i]
      end
      s = R()
      for j = 1:i
         s = addmul_delayed_reduction!(s, Y[i, j], M[i - 1, j], p)
      end
      s = reduce!(s)
      A[i] = s
      for j = 1:i
         s = -F[j]
         for k = 1:j - 1
            s = addmul_delayed_reduction!(s, A[k], F[j - k], p)
         end
         s = reduce!(s)
         F[j] = -s - A[j]
     end
   end
   z = gen(V)
   f = z^n
   for i = 1:n
      f = setcoeff!(f, n - i, F[i])
   end
   return f
end

###############################################################################
#
#   Minimal polynomial
#
###############################################################################

# We let the rows of A be the Krylov sequence v, Mv, M^2v, ... where v
# is a standard basis vector. We find a minimal polynomial P_v by finding
# a linear relation amongst the rows (rows are added one at a time until
# a relation is found). We then append the rows from A to B and row reduce
# by all the rows already in B. We then look for a new standard basis vector
# v not in the span of the rows in B and repeat the above. The arrays P1 and
# P2 are because we can't efficiently swap rows. The i-th entry encodes the
# row of A (in the case of P1) or B (in the case of P2) which has a pivot in
# column i. The arrays L1 and L2 encode the number of columns that contain
# nonzero entries for each row of A and B respectively. These are required
# by the reduce_row! function. If charpoly_only is set to true, the function
# will return just the polynomial obtained from the first Krylov subspace.
# If the charpoly is the minpoly, this returned polynomial will be the
# charpoly iff it has degree n. Otherwise it is meaningless (but it is
# extremely fast to compute over some fields).

@doc Markdown.doc"""
    minpoly(S::Ring, M::MatElem{T}, charpoly_only::Bool = false) where {T <: FieldElement}
> Return the minimal polynomial $p$ of the matrix $M$. The polynomial ring $S$
> of the resulting polynomial must be supplied and the matrix must be square.
"""
function minpoly(S::Ring, M::MatElem{T}, charpoly_only::Bool = false) where {T <: FieldElement}
   !issquare(M) && error("Not a square matrix in minpoly")
   base_ring(S) != base_ring(M) && error("Unable to coerce polynomial")
   n = nrows(M)
   if n == 0
      return S(1)
   end
   R = base_ring(M)
   p = S(1)
   A = similar(M, n + 1, 2n + 1)
   B = similar(M, n, n)
   L1 = [n + i for i in 1:n + 1]
   L2 = [n for i in 1:n]
   P2 = zeros(Int, n)
   P2[1] = 1
   c2 = 1
   r2 = 1
   first_poly = true
   while r2 <= n
      P1 = [0 for i in 1:2n + 1]
      v = similar(M, n, 1)
      for j = 1:n
         B[r2, j] = v[j, 1]
         A[1, j] = R()
      end
      P1[c2] = 1
      P2[c2] = r2
      v[c2, 1] = R(1)
      B[r2, c2] = v[c2, 1]
      A[1, c2] = R(1)
      A[1, n + 1] = R(1)
      indep = true
      r1 = 1
      c1 = 0
      while c1 <= n && r1 <= n
         r1 += 1
         r2 = indep ? r2 + 1 : r2
         v = M*v
         for j = 1:n
            A[r1, j] = deepcopy(v[j, 1])
         end
         for j = n + 1:n + r1 - 1
            A[r1, j] = R(0)
         end
         A[r1, n + r1] = R(1)
         c1 = reduce_row!(A, P1, L1, r1)
         if indep && r2 <= n && !first_poly
            for j = 1:n
               B[r2, j] = deepcopy(v[j, 1])
            end
            c = reduce_row!(B, P2, L2, r2)
            indep = c != 0
         end
      end
      if first_poly
         for j = 1:n
            P2[j] = P1[j]
         end
         r2 = r1
      end
      c = 0
      for j = c2 + 1:n
         if P2[j] == 0
            c = j
            break
         end
      end
      c2 = c
      b = S()
      fit!(b, r1)
      h = inv(A[r1, n + r1])
      for i = 1:r1
         b = setcoeff!(b, i - 1, A[r1, n + i]*h)
      end
      p = lcm(p, b)
      if charpoly_only == true
         return p
      end
      if first_poly && r2 <= n
         for j = 1:r1 - 1
            for k = 1:n
               B[j, k] = A[j, k]
            end
         end
      end
      first_poly = false
   end
   return p
end

@doc Markdown.doc"""
    minpoly(S::Ring, M::MatElem{T}, charpoly_only::Bool = false) where {T <: RingElement}
> Return the minimal polynomial $p$ of the matrix $M$. The polynomial ring $S$
> of the resulting polynomial must be supplied and the matrix must be square.
"""
function minpoly(S::Ring, M::MatElem{T}, charpoly_only::Bool = false) where {T <: RingElement}
   !issquare(M) && error("Not a square matrix in minpoly")
   base_ring(S) != base_ring(M) && error("Unable to coerce polynomial")
   n = nrows(M)
   if n == 0
      return S(1)
   end
   R = base_ring(M)
   p = S(1)
   A = similar(M, n + 1, 2n + 1)
   B = similar(M, n, n)
   L1 = zeros(Int, n + 1)
   for i in 1:n + 1
      L1[i] = i + n
   end
   L2 = fill(n, n)
   P2 = zeros(Int, n)
   P2[1] = 1
   c2 = 1
   r2 = 1
   first_poly = true
   while r2 <= n
      P1 = [0 for i in 1:2n + 1]
      v = similar(M, n, 1)
      for j = 1:n
         B[r2, j] = v[j, 1]
         A[1, j] = R()
      end
      P1[c2] = 1
      P2[c2] = r2
      v[c2, 1] = R(1)
      B[r2, c2] = v[c2, 1]
      for s = 1:c2 - 1
         if P2[s] != 0
            B[r2, c2] *= B[P2[s], s]
         end
      end
      A[1, c2] = R(1)
      A[1, n + 1] = R(1)
      indep = true
      r1 = 1
      c1 = 0
      while c1 <= n && r1 <= n
         r1 += 1
         r2 = indep ? r2 + 1 : r2
         v = M*v
         for j = 1:n
            A[r1, j] = deepcopy(v[j, 1])
         end
         for j = n + 1:n + r1 - 1
            A[r1, j] = R(0)
         end
         A[r1, n + r1] = R(1)
         c1 = reduce_row!(A, P1, L1, r1)
         if indep && r2 <= n && !first_poly
            for j = 1:n
               B[r2, j] = deepcopy(v[j, 1])
            end
            c = reduce_row!(B, P2, L2, r2)
            indep = c != 0
         end
      end
      if first_poly
         for j = 1:n
            P2[j] = P1[j]
         end
         r2 = r1
      end
      c = 0
      for j = c2 + 1:n
         if P2[j] == 0
            c = j
            break
         end
      end
      c2 = c
      b = S()
      fit!(b, r1)
      for i = 1:r1
         b = setcoeff!(b, i - 1, A[r1, n + i])
      end
      b = reverse(b, r1)
      b = primpart(b)
      b = reverse(b, r1)
      p = lcm(p, b)
      if charpoly_only == true
         return divexact(p, canonical_unit(p))
      end
      if first_poly && r2 <= n
         for j = 1:r1 - 1
            for k = 1:n
               B[j, k] = A[j, k]
            end
         end
      end
      first_poly = false
   end
   return divexact(p, canonical_unit(p))
end

###############################################################################
#
#   Hermite Normal Form
#
###############################################################################

function hnf_cohen(A::MatrixElem{T}) where {T <: RingElement}
   H, U = hnf_cohen_with_transform(A)
   return H
end

function hnf_cohen_with_transform(A::MatrixElem{T}) where {T <: RingElement}
   H = deepcopy(A)
   m = nrows(H)
   U = eye(A, m)
   hnf_cohen!(H, U)
   return H, U
end

function hnf_cohen!(H::MatrixElem{T}, U::MatrixElem{T}) where {T <: RingElement}
   m = nrows(H)
   n = ncols(H)
   l = min(m, n)
   k = 1
   t = base_ring(H)()
   t1 = base_ring(H)()
   t2 = base_ring(H)()
   for i = 1:l
      for j = k+1:m
         if iszero(H[j,i])
            continue
         end
         d, u, v = gcdx(H[k,i], H[j,i])
         a = divexact(H[k,i], d)
         b = -divexact(H[j,i], d)
         for c = i:n
            t = deepcopy(H[j,c])
            t1 = mul_red!(t1, a, H[j, c], false)
            t2 = mul_red!(t2, b, H[k, c], false)
            H[j, c] = add!(H[j, c], t1, t2)
            H[j, c] = reduce!(H[j, c])
            t1 = mul_red!(t1, u, H[k, c], false)
            t2 = mul_red!(t2, v, t, false)
            H[k, c] = add!(H[k, c], t1, t2)
            H[k, c] = reduce!(H[k, c])
         end
         for c = 1:m
            t = deepcopy(U[j,c])
            t1 = mul_red!(t1, a, U[j, c], false)
            t2 = mul_red!(t2, b, U[k, c], false)
            U[j, c] = add!(U[j, c], t1, t2)
            U[j, c] = reduce!(U[j, c])
            t1 = mul_red!(t1, u, U[k, c], false)
            t2 = mul_red!(t2, v, t, false)
            U[k, c] = add!(U[k, c], t1, t2)
            U[k, c] = reduce!(U[k, c])
         end
      end
      if iszero(H[k, i])
         continue
      end
      cu = canonical_unit(H[k, i])
      if cu != 1
         for c = i:n
            H[k, c] = divexact(H[k, c],cu)
        end
         for c = 1:m
            U[k, c] = divexact(U[k, c],cu)
         end
      end
      for j = 1:k-1
         q = -AbstractAlgebra.div(H[j,i], H[k, i])
         for c = i:n
            t = mul!(t, q, H[k, c])
            H[j, c] = addeq!(H[j, c], t)
         end
         for c = 1:m
            t = mul!(t, q, U[k, c])
            U[j, c] = addeq!(U[j, c], t)
         end
      end
      k += 1
   end
   return nothing
end

#  Hermite normal form via Kannan-Bachem for matrices of full column rank
#
#  This is the algorithm of Kannan, Bachem, "Polynomial algorithms for computing
#  the Smith and Hermite normal forms of an integer matrix", Siam J. Comput.,
#  Vol. 8, No. 4, pp. 499-507.

@doc Markdown.doc"""
    hnf_minors(A::Generic.MatrixElem{T}) where {T <: RingElement}
> Compute the upper right row Hermite normal form of $A$ using the algorithm of
> Kannan-Bachem. The input must have full column rank.
"""
function hnf_minors(A::MatrixElem{T}) where {T <: RingElement}
   H = deepcopy(A)
   _hnf_minors!(H, similar(A, 0, 0), Val{false})
   return H
end

@doc Markdown.doc"""
    hnf_minors_with_transform(A::Generic.MatrixElem{T}) where {T <: RingElement}
> Compute the upper right row Hermite normal form $H$ of $A$ and an invertible
> matrix $U$ with $UA = H$ using the algorithm of Kannan-Bachem. The input must
> have full column rank.
"""
function hnf_minors_with_transform(A::MatrixElem{T}) where {T <: RingElement}
   H = deepcopy(A)
   U = similar(A, nrows(A), nrows(A))
   _hnf_minors!(H, U, Val{true})
   return H, U
end

function _hnf_minors!(H::MatrixElem{T}, U::MatrixElem{T}, with_transform::Type{Val{S}} = Val{false}) where {T <: RingElement, S}
   m = nrows(H)
   n = ncols(H)

   l = m
   k = 0

   R = base_ring(H)

   with_trafo = with_transform == Val{true} ? true : false

   if with_trafo
      for i in 1:m
         for j in 1:m
            if j == i
               U[i, j] = one(R)
            else
               U[i, j] = zero(R)
            end
         end
      end
   end

   t = zero(R)
   b = zero(R)
   t2 = zero(R)
   r1d = zero(R)
   r2d = zero(R)
   q = zero(R)
   minus_one = base_ring(H)(-1)

   while k <= n - 1
      k = k + 1
      if k == 1 && !iszero(H[k, k])
         for j2 in 1:n
            H[k, j2] = deepcopy(H[k, j2])
         end
      end

      # We have a matrix of form ( * * * * )
      #                          ( 0 * * * )
      #                          ( 0 0 * * ) <- k - 1
      #                          ( * * * * ) <- k
      # We want to produce zeroes in the first k entries
      # of row k.

      for j in 1:k-1
         # Since we want to mutate the k-th row, first make a copy
         if j == 1
            for j2 in j:n
               H[k, j2] = deepcopy(H[k, j2])
            end
         end

         # Shortcuts
         if iszero(H[k, j])
            continue
         end

         di, q = divides(H[k, j], H[j, j])
         if di
            q = mul!(q, q, minus_one)
            for j2 in j:n
               t = mul!(t, q, H[j, j2])
               H[k, j2] = add!(H[k, j2], H[k, j2], t)
            end
            if with_trafo
               for j2 in 1:m
                  t = mul!(t, q, U[j, j2])
                  U[k, j2] = add!(U[k, j2], U[k, j2], t)
               end
            end
            continue
         end

         # Generic case
         d, u, v = gcdx(H[j, j], H[k, j])

         r1d = divexact(H[j, j], d)
         r2d = divexact(H[k, j], d)

         mul!(r2d, r2d, minus_one)
         for j2 in j:n
            b = mul_red!(b, u, H[j, j2], false)
            t2 = mul_red!(t2, v, H[k, j2], false)
            H[k, j2] = mul_red!(H[k, j2], H[k, j2], r1d, false)
            t = mul_red!(t, r2d, H[j, j2], false)
            H[k, j2] = add!(H[k, j2], H[k, j2], t)
            H[j, j2] = add!(H[j, j2], b, t2)
            H[k, j2] = reduce!(H[k, j2])
            H[j, j2] = reduce!(H[j, j2])
         end
         if with_trafo
            for j2 in 1:m
               b = mul_red!(b, u, U[j, j2], false)
               t2 = mul_red!(t2, v, U[k, j2], false)
               U[k, j2] = mul_red!(U[k, j2], U[k, j2], r1d, false)
               t = mul_red!(t, r2d, U[j, j2], false)
               U[k, j2] = add!(U[k, j2], U[k, j2], t)
               U[j, j2] = add!(U[j, j2], b, t2)
               U[k, j2] = reduce!(U[k, j2])
               U[j, j2] = reduce!(U[j, j2])
            end
         end
      end

      if iszero(H[k, k])
         swap_rows!(H, k, l)
         if with_trafo
            swap_rows!(U, k, l)
         end
         l = l - 1
         k = k - 1
         continue
      end

      u = canonical_unit(H[k, k])
      if !isone(u)
         u = R(inv(u))
         for j in k:n
            H[k, j] = mul!(H[k, j], H[k, j], u)
         end
         if with_trafo
            for j in 1:m
               U[k, j] = mul!(U[k, j], U[k, j], u)
            end
         end
      end

      for i in (k - 1):-1:1
         for j in (i + 1):k
            q = AbstractAlgebra.div(H[i, j], H[j, j])
            if iszero(q)
              continue
            end
            q = mul!(q, q, minus_one)
            for j2 in j:n
               t = mul!(t, q, H[j, j2])
               H[i, j2] = add!(H[i, j2], H[i, j2], t)
            end
            if with_trafo
               for j2 in 1:m
                  t = mul!(t, q, U[j, j2])
                  U[i, j2] = add!(U[i, j2], U[i, j2], t)
               end
            end
         end
      end
      l = m
   end

   # Now the matrix is of form
   # ( * * * )
   # ( 0 * * )
   # ( 0 0 * )
   # ( * * * ) <- n + 1
   # ( * * * )

   for k in (n + 1):m
      for j in 1:n
         if j == 1
            for j2 in 1:n
               H[k, j2] = deepcopy(H[k, j2])
            end
         end

         # We do the same shortcuts as above
         if iszero(H[k, j])
            continue
         end

         di, q = divides(H[k, j], H[j, j])
         if di
            q = mul!(q, q, minus_one)
            for j2 in j:n
               t = mul!(t, q, H[j, j2])
               H[k, j2] = add!(H[k, j2], H[k, j2], t)
            end
            if with_trafo
               for j2 in 1:m
                  t = mul!(t, q, U[j, j2])
                  U[k, j2] = add!(U[k, j2], U[k, j2], t)
               end
            end
            continue
         end

         d, u, v = gcdx(H[j, j], H[k, j])
         r1d = divexact(H[j, j], d)
         r2d = divexact(H[k, j], d)
         mul!(r2d, r2d, minus_one)
         for j2 in j:n
            b = mul_red!(b, u, H[j, j2], false)
            t2 = mul_red!(t2, v, H[k, j2], false)
            H[k, j2] = mul_red!(H[k, j2], H[k, j2], r1d, false)
            t = mul_red!(t, r2d, H[j, j2], false)
            H[k, j2] = add!(H[k, j2], H[k, j2], t)
            H[j, j2] = add!(H[j, j2], b, t2)
            H[k, j2] = reduce!(H[k, j2])
            H[j, j2] = reduce!(H[j, j2])
         end
         if with_trafo
            for j2 in 1:m
               b = mul_red!(b, u, U[j, j2], false)
               t2 = mul_red!(t2, v, U[k, j2], false)
               U[k, j2] = mul_red!(U[k, j2], U[k, j2], r1d, false)
               t = mul_red!(t, r2d, U[j, j2], false)
               U[k, j2] = add!(U[k, j2], U[k, j2], t)
               U[j, j2] = add!(U[j, j2], b, t2)
               U[k, j2] = reduce!(U[k, j2])
               U[j, j2] = reduce!(U[j, j2])
            end
         end
      end
      for i in n:-1:1
         for j in (i + 1):n
            q = AbstractAlgebra.div(H[i, j], H[j, j])
            if iszero(q)
              continue
            end
            q = mul!(q, q, minus_one)
            for j2 in j:n
               t = mul!(t, q, H[j, j2])
               H[i, j2] = add!(H[i, j2], H[i, j2], t)
            end
            if with_trafo
               for j2 in 1:m
                  t = mul!(t, q, U[j, j2])
                  U[i, j2] = add!(U[i, j2], U[i, j2], t)
               end
            end
         end
      end
   end
   return H
end

#  Hermite normal form for arbitrary matrices via a modification of the
#  Kannan-Bachem algorithm

@doc Markdown.doc"""
    hnf_kb(A::Generic.MatrixElem{T}) where {T <: RingElement}
> Compute the upper right row Hermite normal form of $A$ using a modification
> of the algorithm of Kannan-Bachem.
"""
function hnf_kb(A::MatrixElem{T}) where {T <: RingElement}
   return _hnf_kb(A, Val{false})
end

@doc Markdown.doc"""
    hnf_kb_with_transform(A::Generic.MatrixElem{T}) where {T <: RingElement}
> Compute the upper right row Hermite normal form $H$ of $A$ and an invertible
> matrix $U$ with $UA = H$ using a modification of the algorithm of
> Kannan-Bachem.
"""
function hnf_kb_with_transform(A::MatrixElem{T}) where {T <: RingElement}
   return _hnf_kb(A, Val{true})
end

function _hnf_kb(A, trafo::Type{Val{T}} = Val{false}) where T
   H = deepcopy(A)
   m = nrows(H)
   if trafo == Val{true}
      U = eye(A, m)
      hnf_kb!(H, U, true)
      return H, U
   else
      U = similar(A, 0, 0)
      hnf_kb!(H, U, false)
      return H
   end
end

function kb_search_first_pivot(H, start_element::Int = 1)
   for r = start_element:nrows(H)
      for c = start_element:ncols(H)
         if !iszero(H[r,c])
            return r, c
         end
      end
   end
   return 0, 0
end

function kb_reduce_row!(H::MatrixElem{T}, U::MatrixElem{T}, pivot::Array{Int, 1}, c::Int, with_trafo::Bool) where {T <: RingElement}
   r = pivot[c]
   t = base_ring(H)()
   for i = c+1:ncols(H)
      p = pivot[i]
      if p == 0
         continue
      end
      q = -AbstractAlgebra.div(H[r,i], H[p,i])
      for j = i:ncols(H)
         t = mul!(t, q, H[p,j])
         H[r, j] = addeq!(H[r,j], t)
      end
      if with_trafo
         for j = 1:ncols(U)
            t = mul!(t, q, U[p,j])
            U[r, j] = addeq!(U[r,j], t)
         end
      end
   end
   return nothing
end

function kb_reduce_column!(H::MatrixElem{T}, U::MatrixElem{T}, pivot::Array{Int, 1}, c::Int, with_trafo::Bool, start_element::Int = 1) where {T <: RingElement}
   r = pivot[c]
   t = base_ring(H)()
   for i = start_element:c-1
      p = pivot[i]
      if p == 0
         continue
      end
      q = -AbstractAlgebra.div(H[p,c],H[r,c])
      for j = c:ncols(H)
         t = mul!(t, q, H[r,j])
         H[p, j] = addeq!(H[p,j], t)
      end
      if with_trafo
         for j = 1:ncols(U)
            t = mul!(t, q, U[r,j])
            U[p, j] = addeq!(U[p,j], t)
         end
      end
   end
   return nothing
end

function kb_canonical_row!(H, U, r::Int, c::Int, with_trafo::Bool)
   cu = canonical_unit(H[r,c])
   if cu != 1
      for j = c:ncols(H)
         H[r,j] = divexact(H[r,j],cu)
      end
      if with_trafo
         for j = 1:ncols(U)
            U[r,j] = divexact(U[r,j],cu)
         end
      end
   end
   return nothing
end

function kb_sort_rows!(H::MatrixElem{T}, U::MatrixElem{T}, pivot::Array{Int, 1}, with_trafo::Bool, start_element::Int = 1) where {T <:RingElement}
   m = nrows(H)
   n = ncols(H)
   pivot2 = zeros(Int, m)
   for i = 1:n
      if pivot[i] == 0
         continue
      end
      pivot2[pivot[i]] = i
   end

   r1 = start_element
   for i = start_element:n
      r2 = pivot[i]
      if r2 == 0
         continue
      end
      if r1 != r2
         swap_rows!(H, r1, r2)
         with_trafo ? swap_rows!(U, r1, r2) : nothing
         p = pivot2[r1]
         pivot[i] = r1
         if p != 0
            pivot[p] = r2
         end
         pivot2[r1] = i
         pivot2[r2] = p
      end
      r1 += 1
      if r1 == m
         break
      end
   end
   return nothing
end

function hnf_kb!(H, U, with_trafo::Bool = false, start_element::Int = 1)
   m = nrows(H)
   n = ncols(H)
   pivot = zeros(Int, n)
   row1, col1 = kb_search_first_pivot(H, start_element)
   if row1 == 0
      return nothing
   end
   pivot[col1] = row1
   kb_canonical_row!(H, U, row1, col1, with_trafo)
   pivot_max = col1
   t = base_ring(H)()
   t1 = base_ring(H)()
   t2 = base_ring(H)()
   for i=row1:m-1
      new_pivot = false
      for j = start_element:pivot_max
         if iszero(H[i+1,j])
            continue
         end
         if pivot[j] == 0
            pivot[j] = i+1
            kb_canonical_row!(H, U, pivot[j], j, with_trafo)
            kb_reduce_column!(H, U, pivot, j, with_trafo, start_element)
            kb_reduce_row!(H, U, pivot, j, with_trafo)
            pivot_max = max(pivot_max, j)
            new_pivot = true
         else
            p = pivot[j]
            d, u, v = gcdx(H[p, j], H[i + 1, j])
            a = divexact(H[p, j], d)
            b = -divexact(H[i + 1, j], d)
            for c = j:n
               t = deepcopy(H[i + 1, c])
               t1 = mul_red!(t1, a, H[i + 1, c], false)
               t2 = mul_red!(t2, b, H[p, c], false)
               H[i + 1, c] = add!(H[i + 1, c], t1, t2)
               H[i + 1, c] = reduce!(H[i + 1, c])
               t1 = mul_red!(t1, u, H[p, c], false)
               t2 = mul_red!(t2, v, t, false)
               H[p, c] = add!(H[p, c], t1, t2)
               H[p, c] = reduce!(H[p, c])
            end
            if with_trafo
               for c = 1:m
                  t = deepcopy(U[i + 1, c])
                  t1 = mul_red!(t1, a, U[i + 1, c], false)
                  t2 = mul_red!(t2, b, U[p, c], false)
                  U[i + 1, c] = add!(U[i + 1, c], t1, t2)
                  U[i + 1, c] = reduce!(U[i + 1, c])
                  t1 = mul_red!(t1, u, U[p, c], false)
                  t2 = mul_red!(t2, v, t, false)
                  U[p, c] = add!(U[p, c], t1, t2)
                  U[p, c] = reduce!(U[p, c])
               end
            end
            kb_canonical_row!(H, U, pivot[j], j, with_trafo)
            kb_reduce_column!(H, U, pivot, j, with_trafo, start_element)
         end
         if new_pivot
            break
         end
      end
      if !new_pivot
         for c = pivot_max+1:n
            if !iszero(H[i+1,c])
               pivot[c] = i+1
               kb_canonical_row!(H, U, pivot[c], c, with_trafo)
               kb_reduce_column!(H, U, pivot, c, with_trafo, start_element)
               pivot_max = max(pivot_max, c)
               break
            end
         end
      end
   end
   kb_sort_rows!(H, U, pivot, with_trafo, start_element)
   return nothing
end

@doc Markdown.doc"""
    hnf(A::Generic.MatrixElem{T}) where {T <: RingElement}
> Return the upper right row Hermite normal form of $A$.
"""
function hnf(A::MatrixElem{T}) where {T <: RingElement}
  return hnf_kb(A)
end

@doc Markdown.doc"""
    hnf_with_transform(A)
> Return the tuple $H, U$ consisting of the upper right row Hermite normal
> form $H$ of $A$ together with invertible matrix $U$ such that $UA = H$.
"""
function hnf_with_transform(A)
  return hnf_kb_with_transform(A)
end

@doc Markdown.doc"""
    ishnf(M::Generic.MatrixElem{T}) where T <: RingElement
> Return `true` if the matrix is in Hermite normal form.
"""
function ishnf(M::Generic.MatrixElem{T}) where T <: RingElement
   r = nrows(M)
   c = ncols(M)
   row = 1
   col = 1
   pivots = zeros(Int, r)
   # first check the staircase, since it is cheap to do
   while row <= r
      while col <= c && M[row, col] == 0
         for i = row + 1:r
            if M[i, col] != 0
               return false
            end
         end
         col += 1
      end
      if col <= c # found pivot
         pivots[row] = col
         for i = row + 1:r
            if M[i, col] != 0
               return false
            end
         end
      end
      row += 1
      col += 1
   end
   # now check everything above pivots is reduced (expensive)
   row = 1
   while row <= r && pivots[row] != 0
      col = pivots[row]
      p = M[row, col]
      for i = 1:row - 1
         qq, rr = AbstractAlgebra.divrem(M[i, col], p)
         if rr != M[i, col]
            return false
         end
      end
      row += 1
   end
   return true
end

###############################################################################
#
#   Smith Normal Form
#
###############################################################################

function snf_kb(A::MatrixElem{T}) where {T <: RingElement}
   return _snf_kb(A, Val{false})
end

function snf_kb_with_transform(A::MatrixElem{T}) where {T <: RingElement}
   return _snf_kb(A, Val{true})
end

function _snf_kb(A::MatrixElem{T}, trafo::Type{Val{V}} = Val{false}) where {V, T <: RingElement}
   S = deepcopy(A)
   m = nrows(S)
   n = ncols(S)
   if trafo == Val{true}
      U = eye(A, m)
      K = eye(A, n)
      snf_kb!(S, U, K, true)
      return S, U, K
   else
      U = similar(A, 0, 0)
      K = U
      snf_kb!(S, U, K, false)
      return S
   end
end

function kb_clear_row!(S::MatrixElem{T}, K::MatrixElem{T}, i::Int, with_trafo::Bool) where {T <: RingElement}
   m = nrows(S)
   n = ncols(S)
   t = base_ring(S)()
   t1 = base_ring(S)()
   t2 = base_ring(S)()
   for j = i+1:n
      if iszero(S[i, j])
         continue
      end
      d, u, v = gcdx(S[i, i], S[i, j])
      a = divexact(S[i ,i], d)
      b = -divexact(S[i, j], d)
      for r = i:m
         t = deepcopy(S[r, j])
         t1 = mul_red!(t1, a, S[r, j], false)
         t2 = mul_red!(t2, b, S[r, i], false)
         S[r, j] = add!(S[r, j], t1, t2)
         S[r, j] = reduce!(S[r, j])
         t1 = mul_red!(t1, u, S[r, i], false)
         t2 = mul_red!(t2, v, t, false)
         S[r, i] = add!(S[r, i], t1, t2)
         S[r, i] = reduce!(S[r, i])
      end
      if with_trafo
         for r = 1:n
            t = deepcopy(K[r,j])
            t1 = mul_red!(t1, a, K[r, j], false)
            t2 = mul_red!(t2, b, K[r, i], false)
            K[r, j] = add!(K[r, j], t1, t2)
            K[r, j] = reduce!(K[r, j])
            t1 = mul_red!(t1, u, K[r, i], false)
            t2 = mul_red!(t2, v, t, false)
            K[r, i] = add!(K[r, i], t1, t2)
            K[r, i] = reduce!(K[r, i])
         end
      end
   end
   return nothing
end

function snf_kb!(S::MatrixElem{T}, U::MatrixElem{T}, K::MatrixElem{T}, with_trafo::Bool = false) where {T <: RingElement}
   m = nrows(S)
   n = ncols(S)
   l = min(m,n)
   i = 1
   t = base_ring(S)()
   t1 = base_ring(S)()
   t2 = base_ring(S)()
   while i<=l
      kb_clear_row!(S, K, i, with_trafo)
      hnf_kb!(S, U, with_trafo, i)
      c = i+1
      while c <= n && iszero(S[i, c])
         c+=1
      end
      if c != n+1
         continue
      end
      i+=1
   end
   for i = 1:l-1
      for j = i + 1:l
         if isone(S[i, i])
           break
         end
         if iszero(S[i, i]) && iszero(S[j, j])
            continue
         end
         d, u, v = gcdx(S[i, i], S[j, j])
         if with_trafo
            q = -divexact(S[j, j], d)
            t1 = mul!(t1, q, v)
            for c = 1:m
               t = deepcopy(U[i,c])
               U[i, c] = addeq!(U[i, c], U[j,c])
               t2 = mul_red!(t2, t1, U[j, c], false)
               U[j, c] = addeq!(U[j, c], t2)
               t2 = mul_red!(t2, t1, t, false)
               U[j, c] = addeq!(U[j, c], t2)
               U[j, c] = reduce!(U[j, c])
            end
            q1 = -divexact(S[j, j], d)
            q2 = divexact(S[i, i], d)
            for r = 1:n
               t = deepcopy(K[r, i])
               t1 = mul_red!(t1, K[r, i], u, false)
               t2 = mul_red!(t2, K[r, j], v, false)
               K[r, i] = add!(K[r, i], t1, t2)
               K[r, i] = reduce!(K[r, i])
               t1 = mul_red!(t1, t, q1, false)
               t2 = mul_red!(t2, K[r, j], q2, false)
               K[r, j] = add!(K[r, j], t1, t2)
               K[r, j] = reduce!(K[r, j])
            end
         end
         S[j, j] = divexact(S[i, i]*S[j, j],d)
         S[i, i] = d
      end
   end
   return nothing
end

function snf(a::MatrixElem{T}) where {T <: RingElement}
  return snf_kb(a)
end

function snf_with_transform(a::MatrixElem{T}) where {T <: RingElement}
  return snf_kb_with_transform(a)
end

################################################################################
#
#   Popov Form
#
################################################################################

@doc Markdown.doc"""
    weak_popov(A::Mat{T}) where {T <: PolyElem}
> Return the weak Popov form of $A$.
"""
function weak_popov(A::Mat{T}) where {T <: PolyElem}
   return _weak_popov(A, Val{false})
end

@doc Markdown.doc"""
    weak_popov_with_transform(A::Mat{T}) where {T <: PolyElem}
> Compute a tuple $(P, U)$ where $P$ is the weak Popov form of $A$ and $U$
> is a transformation matrix so that $P = UA$.
"""
function weak_popov_with_transform(A::Mat{T}) where {T <: PolyElem}
   return _weak_popov(A, Val{true})
end

function _weak_popov(A::Mat{T}, trafo::Type{Val{S}} = Val{false}) where {T <: PolyElem, S}
   P = deepcopy(A)
   m = nrows(P)
   W = similar(A, 0, 0)
   if trafo == Val{true}
      U = eye(A, m)
      weak_popov!(P, W, U, false, true)
      return P, U
   else
      U = similar(A, 0, 0)
      weak_popov!(P, W, U, false, false)
      return P
   end
end

@doc Markdown.doc"""
    extended_weak_popov(A::Mat{T}, V::Mat{T}) where {T <: PolyElem}
> Compute the weak Popov form $P$ of $A$ by applying simple row transformations
> on $A$ and a vector $W$ by applying the same transformations on the vector $V$.
> Return the tuple $(P, W)$.
"""
function extended_weak_popov(A::Mat{T}, V::Mat{T}) where {T <: PolyElem}
   return _extended_weak_popov(A, V, Val{false})
end

@doc Markdown.doc"""
    extended_weak_popov_with_transform(A::Mat{T}, V::Mat{T}) where {T <: PolyElem}
> Compute the weak Popov form $P$ of $A$ by applying simple row transformations
> on $A$, a vector $W$ by applying the same transformations on the vector $V$,
> and a transformation matrix $U$ so that $P = UA$.
> Return the tuple $(P, W, U)$.
"""
function extended_weak_popov_with_transform(A::Mat{T}, V::Mat{T}) where {T <: PolyElem}
   return _extended_weak_popov(A, V, Val{true})
end

function _extended_weak_popov(A::Mat{T}, V::Mat{T}, trafo::Type{Val{S}} = Val{false}) where {T <: PolyElem, S}
   @assert nrows(V) == nrows(A) && ncols(V) == 1
   P = deepcopy(A)
   W = deepcopy(V)
   m = nrows(P)
   if trafo == Val{true}
      U = eye(A)
      weak_popov!(P, W, U, true, true)
      return P, W, U
   else
      U = similar(A, 0, 0)
      weak_popov!(P, W, U, true, false)
      return P, W
   end
end

function find_pivot_popov(P::Mat{T}, r::Int, last_col::Int = 0) where {T <: PolyElem}
   last_col == 0 ? n = ncols(P) : n = last_col
   pivot = n
   for c = n-1:-1:1
      if degree(P[r,c]) > degree(P[r,pivot])
         pivot = c
      end
   end
   return pivot
end

function init_pivots_popov(P::Mat{T}, last_row::Int = 0, last_col::Int = 0) where {T <: PolyElem}
   last_row == 0 ? m = nrows(P) : m = last_row
   last_col == 0 ? n = ncols(P) : n = last_col
   pivots = Array{Array{Int,1}}(undef, n)
   for i = 1:n
      pivots[i] = zeros(Int, 0)
   end
   # pivots[i] contains the indices of the rows in which the pivot element is in the ith column.
   for r = 1:m
      pivot = find_pivot_popov(P, r, last_col)
      !iszero(P[r,pivot]) ? push!(pivots[pivot], r) : nothing
   end
   return pivots
end

function weak_popov!(P::Mat{T}, W::Mat{T}, U::Mat{T}, extended::Bool = false,
                                       with_trafo::Bool = false, last_row::Int = 0, last_col::Int = 0) where {T <: PolyElem}
   pivots = init_pivots_popov(P, last_row, last_col)
   weak_popov_with_pivots!(P, W, U, pivots, extended, with_trafo, last_row, last_col)
   return nothing
end

#=
The weak Popov form is defined by T. Mulders and A. Storjohann in
"On lattice reduction for polynomial matrices"
=#
function weak_popov_with_pivots!(P::Mat{T}, W::Mat{T}, U::Mat{T}, pivots::Array{Array{Int,1}},
                                                   extended::Bool = false, with_trafo::Bool = false, last_row::Int = 0, last_col::Int = 0) where {T <: PolyElem}
   last_row == 0 ? m = nrows(P) : m = last_row
   last_col == 0 ? n = ncols(P) : n = last_col
   @assert length(pivots) >= n

   t = base_ring(P)()
   change = true
   while change
      change = false
      for i = 1:n
         if length(pivots[i]) <= 1
            continue
         end
         change = true
         # Reduce with the pivot of minimal degree.
         # TODO: Remove the list comprehension as soon as indmin(f, A) works
         #pivotInd = indmin(degree(P[j, i]) for j in pivots[i])
         pivotInd = argmin([degree(P[j, i]) for j in pivots[i]])
         pivot = pivots[i][pivotInd]
         for j = 1:length(pivots[i])
            if j == pivotInd
               continue
            end
            q = -AbstractAlgebra.div(P[pivots[i][j], i], P[pivot, i])
            for c = 1:n
               t = mul!(t, q, P[pivot, c])
               P[pivots[i][j], c] = addeq!(P[pivots[i][j], c], t)
            end
            if with_trafo
               for c = 1:ncols(U)
                  t = mul!(t, q, U[pivot, c])
                  U[pivots[i][j], c] = addeq!(U[pivots[i][j], c], t)
               end
            end
            if extended
               t = mul!(t, q, W[pivot,1])
               W[pivots[i][j], 1] = addeq!(W[pivots[i][j], 1], t)
            end
         end
         old_pivots = pivots[i]
         pivots[i] = [pivot]
         for j = 1:length(old_pivots)
            if j == pivotInd
               continue
            end
            p = find_pivot_popov(P, old_pivots[j], last_col)
            !iszero(P[old_pivots[j],p]) ? push!(pivots[p], old_pivots[j]) : nothing
         end
      end
   end
   return nothing
end

@doc Markdown.doc"""
    rank_profile_popov(A::Mat{T}) where {T <: PolyElem}
> Return an array of $r$ row indices such that these rows of $A$ are linearly
> independent, where $r$ is the rank of $A$.
"""
function rank_profile_popov(A::Mat{T}) where {T <: PolyElem}
   B = deepcopy(A)
   m = nrows(A)
   n = ncols(A)
   U = similar(A, 0, 0)
   V = U
   r = 0
   rank_profile = Array{Int,1}(undef, 0)
   pivots = Array{Array{Int,1}}(undef, n)
   for i = 1:n
      pivots[i] = zeros(Int, 0)
   end
   p = find_pivot_popov(B, 1)
   if !iszero(B[1,p])
      push!(pivots[p], 1)
      r = 1
      push!(rank_profile, 1)
   end
   for i = 2:m
      p = find_pivot_popov(B, i)
      !iszero(B[i,p]) ? push!(pivots[p], i) : nothing
      weak_popov_with_pivots!(B, V, U, pivots, false, false, i)
      s = 0
      for j = 1:n
         # length(pivots[j]) is either 1 or 0, since B is in weak
         # Popov form.
         s += length(pivots[j])
      end
      if s != r
         push!(rank_profile, i)
         r = s
      end
   end
   return rank_profile
end

function det_popov(A::Mat{T}) where {T <: PolyElem}
   nrows(A) != ncols(A) && error("Not a square matrix in det_popov.")
   B = deepcopy(A)
   n = ncols(B)
   R = base_ring(B)
   det = one(R)
   U = similar(A, 0, 0)
   V = U
   t = R()
   pivots1 = init_pivots_popov(B)
   weak_popov_with_pivots!(B, V, U, pivots1, false, false)
   pivots = zeros(Int, n)
   diag_elems = zeros(Int, n)
   for i = 1:n
      if isassigned(pivots1, i) && length(pivots1[i]) > 0
         pivots[i] = pivots1[i][1]
      else
         # If there is no pivot in the ith column, A has not full rank.
         return R(0)
      end
   end
   for i = n-1:-1:1
      # "Remove" the column i+1 and compute a weak Popov Form of the
      # remaining matrix.
      r1 = pivots[i+1]
      c = find_pivot_popov(B, r1, i)
      # If the pivot B[r1, c] is zero then the row is zero.
      while !iszero(B[r1, c])
         r2 = pivots[c]
         if degree(B[r2, c]) > degree(B[r1,c])
            r1, r2 = r2, r1
            pivots[c] = r2
         end
         q = -AbstractAlgebra.div(B[r1, c], B[r2, c])
         for j = 1:i + 1
            t = mul!(t, q, B[r2, j])
            B[r1, j] = addeq!(B[r1, j], t)
         end
         c = find_pivot_popov(B, r1, i)
      end
      if iszero(B[r1, i+1])
         return R(0)
      end
      diag_elems[i+1] = r1
      det = mul!(det, det, B[r1,i+1])
   end
   det = mul!(det, det, B[pivots[1],1])
   diag_elems[1] = pivots[1]
   number_of_swaps = 0
   # Adjust the sign of det by sorting the diagonal elements.
   for i = 1:n
      while diag_elems[i] != i
         r = diag_elems[i]
         diag_elems[i] = diag_elems[r]
         diag_elems[r] = r
         number_of_swaps += 1
      end
   end
   if number_of_swaps%2 == 1
      det = mul!(det, det, R(-1))
   end
   return det
end

@doc Markdown.doc"""
    popov(A::Mat{T}) where {T <: PolyElem}
> Return the Popov form of $A$.
"""
function popov(A::Mat{T}) where {T <: PolyElem}
   return _popov(A, Val{false})
end

@doc Markdown.doc"""
    popov_with_transform(A::Mat{T}) where {T <: PolyElem}
> Compute a tuple $(P, U)$ where $P$ is the Popov form of $A$ and $U$
> is a transformation matrix so that $P = UA$.
"""
function popov_with_transform(A::Mat{T}) where {T <: PolyElem}
   return _popov(A, Val{true})
end

function _popov(A::Mat{T}, trafo::Type{Val{S}} = Val{false}) where {T <: PolyElem, S}
   P = deepcopy(A)
   m = nrows(P)
   if trafo == Val{true}
      U = eye(A, m)
      popov!(P, U, true)
      return P, U
   else
      U = similar(A, 0, 0)
      popov!(P, U, false)
      return P
   end
end

function asc_order_popov!(P::Mat{T}, U::Mat{T}, pivots::Array{Array{Int,1}}, with_trafo::Bool) where {T <: PolyElem}
   m = nrows(P)
   n = ncols(P)
   pivots2 = Array{NTuple{3,Int},1}(undef, m)
   for r = 1:m
      pivots2[r] = (r,n,-1)
   end
   for c = 1:n
      if length(pivots[c]) == 0
         continue
      end
      r = pivots[c][1]
      pivots2[r] = (r, c, degree(P[r,c]))
   end
   sort!(pivots2, lt = (x,y) -> ( x[3] < y[3] || ( x[3] == y[3] && x[2] <= y[2] ) ))
   row_nums = [ i for i = 1:m ]
   for r = 1:m
      if pivots2[r][3] != -1
         c = pivots2[r][2]
         pivots[c] = [r]
      end
      i = pivots2[r][1]
      r2 = row_nums[i]
      if r == r2
         continue
      end
      swap_rows!(P, r, r2)
      with_trafo ? swap_rows!(U, r, r2) : nothing
      j = findfirst(row_nums, r)
      row_nums[i] = r
      row_nums[j] = r2
   end
   return nothing
end

function popov!(P::Mat{T}, U::Mat{T}, with_trafo::Bool = false) where {T <: PolyElem}
   m = nrows(P)
   n = ncols(P)
   W = similar(U, 0, 0)
   pivots = init_pivots_popov(P)
   weak_popov_with_pivots!(P, W, U, pivots, false, with_trafo)
   asc_order_popov!(P, U, pivots, with_trafo)
   t = base_ring(P)()
   for i = 1:n
      if length(pivots[i]) == 0
         continue
      end
      pivot = pivots[i][1]
      d = degree(P[pivot,i])
      for r = 1:pivot-1
         if degree(P[r,i]) < d
            continue
         end
         q = -AbstractAlgebra.div(P[r,i],P[pivot,i])
         for c = 1:n
            t = mul!(t, q, P[pivot,c])
            P[r, c] = addeq!(P[r,c], t)
         end
         if with_trafo
            for c = 1:ncols(U)
               t = mul!(t, q, U[pivot,c])
               U[r, c] = addeq!(U[r,c], t)
            end
         end
      end
   end
   for i = 1:n
      if length(pivots[i]) == 0
         continue
      end
      r = pivots[i][1]
      cu = canonical_unit(P[r,i])
      if cu != 1
         for j = 1:n
            P[r,j] = divexact(P[r,j],cu)
         end
         if with_trafo
            for j = 1:ncols(U)
               U[r,j] = divexact(U[r,j],cu)
            end
         end
      end
   end
   return nothing
end

function hnf_via_popov(A::Mat{T}) where {T <: PolyElem}
   return _hnf_via_popov(A, Val{false})
end

function hnf_via_popov_with_transform(A::Mat{T}) where {T <: PolyElem}
   return _hnf_via_popov(A, Val{true})
end

function _hnf_via_popov(A::Mat{T}, trafo::Type{Val{S}} = Val{false}) where {T <: PolyElem, S}
   H = deepcopy(A)
   m = nrows(H)
   if trafo == Val{true}
      U = eye(A, m)
      hnf_via_popov!(H, U, true)
      return H, U
   else
      U = similar(A, 0, 0)
      hnf_via_popov!(H, U, false)
      return H
   end
end

function hnf_via_popov_reduce_row!(H::Mat{T}, U::Mat{T}, pivots_hermite::Array{Int}, r::Int, with_trafo::Bool) where {T <: PolyElem}
   n = ncols(H)
   t = base_ring(H)()
   for c = 1:n
      if pivots_hermite[c] == 0
         continue
      end
      pivot = pivots_hermite[c]
      q = -AbstractAlgebra.div(H[r, c], H[pivot, c])
      for j = c:n
         t = mul!(t, q, H[pivot, j])
         H[r, j] = addeq!(H[r, j], t)
      end
      if with_trafo
         for j = 1:ncols(U)
            t = mul!(t, q, U[pivot, j])
            U[r, j] = addeq!(U[r, j], t)
         end
      end
   end
   return nothing
end

function hnf_via_popov_reduce_column!(H::Mat{T}, U::Mat{T}, pivots_hermite::Array{Int}, c::Int, with_trafo::Bool) where {T <: PolyElem}
   m = nrows(H)
   n = ncols(H)
   t = base_ring(H)()
   r = pivots_hermite[c]
   for i = 1:m
      if i == r
         continue
      end
      if degree(H[i, c]) < degree(H[r, c])
         continue
      end
      q = -AbstractAlgebra.div(H[i, c], H[r, c])
      for j = 1:n
         t = mul!(t, q, H[r, j])
         H[i, j] = addeq!(H[i, j], t)
      end
      if with_trafo
         for j = 1:ncols(U)
            t = mul!(t, q, U[r, j])
            U[i, j] = addeq!(U[i, j], t)
         end
      end
   end
   return nothing
end

function hnf_via_popov!(H::Mat{T}, U::Mat{T}, with_trafo::Bool = false) where {T <: PolyElem}
   m = nrows(H)
   n = ncols(H)
   R = base_ring(H)
   W = similar(H, 0, 0)
   t = R()
   pivots = init_pivots_popov(H)
   weak_popov_with_pivots!(H, W, U, pivots, false, with_trafo)
   pivots_popov = zeros(Int, n)
   for j = 1:n
      if isassigned(pivots, j) && length(pivots[j]) > 0
         pivots_popov[j] = pivots[j][1]
      else
         error("The rank must be equal to the number of columns.")
      end
   end
   pivots_hermite = zeros(Int, n)
   for i = n-1:-1:1
      # "Remove" the column i+1 and compute a weak Popov Form of the
      # remaining matrix.
      r1 = pivots_popov[i + 1]
      c = find_pivot_popov(H, r1, i)
      new_pivot = true
      # If the pivot H[r1, c] is zero then the row is zero.
      while !iszero(H[r1, c])
         r2 = pivots_popov[c]
         if degree(H[r2, c]) > degree(H[r1,c])
            r1, r2 = r2, r1
            pivots_popov[c] = r2
         end
         q = -AbstractAlgebra.div(H[r1, c], H[r2, c])
         for j = 1:n
            t = mul!(t, q, H[r2, j])
            H[r1, j] = addeq!(H[r1, j], t)
         end
         if with_trafo
            for j = 1:ncols(U)
               t = mul!(t, q, U[r2, j])
               U[r1, j] = addeq!(U[r1, j], t)
            end
         end
         hnf_via_popov_reduce_row!(H, U, pivots_hermite, r1, with_trafo)
         c = find_pivot_popov(H, r1, i)
      end
      new_pivot ? nothing : continue
      pivots_hermite[i+1] = r1
      hnf_via_popov_reduce_column!(H, U, pivots_hermite, i+1, with_trafo)
      l = pivots_popov[i]
      hnf_via_popov_reduce_row!(H, U, pivots_hermite, l, with_trafo)
   end
   pivots_hermite[1] = pivots_popov[1]
   kb_sort_rows!(H, U, pivots_hermite, with_trafo)
   for c = 1:n
      kb_canonical_row!(H, U, pivots_hermite[c], c, with_trafo)
   end
   return nothing
end

###############################################################################
#
#   Transforms
#
###############################################################################

@doc Markdown.doc"""
    similarity!(A::Generic.MatrixElem{T}, r::Int, d::T) where {T <: RingElement}
> Applies a similarity transform to the $n\times n$ matrix $M$ in-place. Let
> $P$ be the $n\times n$ identity matrix that has had all zero entries of row
> $r$ replaced with $d$, then the transform applied is equivalent to
> $M = P^{-1}MP$. We require $M$ to be a square matrix. A similarity transform
> preserves the minimal and characteristic polynomials of a matrix.
"""
function similarity!(A::MatrixElem{T}, r::Int, d::T) where {T <: RingElement}
   n = nrows(A)
   t = base_ring(A)()
   for i = 1:n
      for j = 1:r - 1
         t = mul!(t, A[i, r], d)
         A[i, j] = addeq!(A[i, j], t)
      end
      for j = r + 1:n
         t = mul!(t, A[i, r], d)
         A[i, j] = addeq!(A[i, j], t)
      end
   end
   d = -d
   for i = 1:n
      for j = 1:r - 1
         A[r, i] = addmul_delayed_reduction!(A[r, i], A[j, i], d, t)
      end
      for j = r + 1:n
         A[r, i] = addmul_delayed_reduction!(A[r, i], A[j, i], d, t)
      end
      A[r, i] = reduce!(A[r, i])
   end
end

###############################################################################
#
#   Row swapping
#
###############################################################################

@doc Markdown.doc"""
    swap_rows(a::MatrixElem, i::Int, j::Int)
> Return a matrix $b$ with the entries of $a$, where the $i$th and $j$th
> row are swapped.
"""
function swap_rows(a::MatrixElem, i::Int, j::Int)
   (1 <= i <= nrows(a) && 1 <= j <= nrows(a)) || throw(BoundsError())
   b = deepcopy(a)
   swap_rows!(b, i, j)
   return b
end

@doc Markdown.doc"""
    swap_rows!(a::MatrixElem, i::Int, j::Int)
> Swap the $i$th and $j$th row of $a$.
"""
function swap_rows!(a::MatrixElem, i::Int, j::Int)
   (1 <= i <= nrows(a) && 1 <= j <= nrows(a)) || throw(BoundsError())
   if i != j
      for k = 1:ncols(a)
         x = a[i, k]
         a[i, k] = a[j, k]
         a[j, k] = x
      end
   end
   return a
end

@doc Markdown.doc"""
    swap_cols(a::MatrixElem, i::Int, j::Int)
> Return a matrix $b$ with the entries of $a$, where the $i$th and $j$th
> row are swapped.
"""
function swap_cols(a::MatrixElem, i::Int, j::Int)
   (1 <= i <= ncols(a) && 1 <= j <= ncols(a)) || throw(BoundsError())
   b = deepcopy(a)
   swap_cols!(b, i, j)
   return b
end

@doc Markdown.doc"""
    swap_cols!(a::MatrixElem, i::Int, j::Int)
> Swap the $i$th and $j$th column of $a$.
"""
function swap_cols!(a::MatrixElem, i::Int, j::Int)
   if i != j
      for k = 1:nrows(a)
         x = a[k, i]
         a[k, i] = a[k, j]
         a[k, j] = x
      end
   end
   return a
end

@doc Markdown.doc"""
    invert_rows!(a::MatrixElem)
> Swap the $i$th and $r - i$th row of $a$ for $1 \leq i \leq r/2$,
> where $r$ is the number of rows of $a$.
"""
function invert_rows!(a::MatrixElem)
   k = AbstractAlgebra.div(nrows(a), 2)
   for i in 1:k
      swap_rows!(a, i, nrows(a) - i + 1)
   end
   return a
end

@doc Markdown.doc"""
    invert_rows(a::MatrixElem)
> Return a matrix $b$ with the entries of $a$, where the $i$th and $r - i$th
> row is swapped for $1 \leq i \leq r/2$. Here $r$ is the number of rows of
> $a$.
"""
function invert_rows(a::MatrixElem)
   b = deepcopy(a)
   return invert_rows!(b)
end

@doc Markdown.doc"""
    invert_cols!(a::MatrixElem)
> Swap the $i$th and $r - i$th column of $a$ for $1 \leq i \leq c/2$,
> where $c$ is the number of columns of $a$.
"""
function invert_cols!(a::MatrixElem)
   k = AbstractAlgebra.div(ncols(a), 2)
   for i in 1:k
      swap_cols!(a, i, ncols(a) - i + 1)
   end
   return a
end

@doc Markdown.doc"""
    invert_cols(a::MatrixElem)
> Return a matrix $b$ with the entries of $a$, where the $i$th and $r - i$th
> column is swapped for $1 \leq i \leq c/2$. Here $c$ is the number of columns
> of$a$.
"""
function invert_cols(a::MatrixElem)
   b = deepcopy(a)
   return invert_cols!(b)
end

###############################################################################
#
#   Concatenation
#
###############################################################################

@doc Markdown.doc"""
    hcat(a::AbstractAlgebra.MatElem, b::AbstractAlgebra.MatElem)
> Return the horizontal concatenation of $a$ and $b$. Assumes that the
> number of rows is the same in $a$ and $b$.
"""
function hcat(a::AbstractAlgebra.MatElem, b::AbstractAlgebra.MatElem)
   nrows(a) != nrows(b) && error("Incompatible number of nrows in hcat")
   c = similar(a, nrows(a), ncols(a) + ncols(b))
   n = ncols(a)
   for i = 1:nrows(a)
      for j = 1:ncols(a)
         c[i, j] = a[i, j]
      end
      for j = 1:ncols(b)
         c[i, n + j] = b[i, j]
      end
   end
   return c
end

@doc Markdown.doc"""
    vcat(a::AbstractAlgebra.MatElem, b::AbstractAlgebra.MatElem)
> Return the vertical concatenation of $a$ and $b$. Assumes that the
> number of columns is the same in $a$ and $b$.
"""
function vcat(a::AbstractAlgebra.MatElem, b::AbstractAlgebra.MatElem)
   ncols(a) != ncols(b) && error("Incompatible number of columns in vcat")
   c = similar(a, nrows(a) + nrows(b), ncols(a))
   n = nrows(a)
   for i = 1:nrows(a)
      for j = 1:ncols(a)
         c[i, j] = a[i, j]
      end
   end
   for i = 1:nrows(b)
      for j = 1:ncols(a)
         c[n + i, j] = b[i, j]
      end
   end
   return c
end

###############################################################################
#
#   Random generation
#
###############################################################################

function rand(S::AbstractAlgebra.MatSpace, v...)
   M = S()
   R = base_ring(S)
   for i = 1:nrows(M)
      for j = 1:ncols(M)
         M[i, j] = rand(R, v...)
      end
   end
   return M
end

function randmat_triu(S::AbstractAlgebra.MatSpace, v...)
   M = S()
   R = base_ring(S)
   for i = 1:nrows(M)
      for j = 1:i - 1
         M[i, j] = R()
      end
      for j = i:ncols(M)
         M[i, j] = rand(R, v...)
      end
      while iszero(M[i, i])
         M[i, i] = rand(R, v...)
      end
   end
   return M
end

function randmat_with_rank(S::Generic.MatSpace{T}, rank::Int, v...) where {T <: AbstractAlgebra.RingElement}
   if !isdomain_type(T) && !(T <: ResElem)
      error("Not implemented")
   end
   M = S()
   R = base_ring(S)
   for i = 1:rank
      for j = 1:i - 1
         M[i, j] = R()
      end
      M[i, i] = rand(R, v...)
      while iszero(M[i, i])
         M[i, i] = rand(R, v...)
      end
      for j = i + 1:ncols(M)
         M[i, j] = rand(R, v...)
      end
   end
   for i = rank + 1:nrows(M)
      for j = 1:ncols(M)
         M[i, j] = R()
      end
   end
   m = nrows(M)
   if m > 1
      for i = 1:4*m
         r1 = rand(1:m)
         r2 = rand(1:m - 1)
         r2 = r2 >= r1 ? r2 + 1 : r2
         d = rand(-5:5)
         for j = 1:ncols(M)
            M[r1, j] = M[r1, j] + d*M[r2, j]
         end
      end
   end
   return M
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{S}, ::Type{S}) where {T <: RingElement, S <: Mat{T}} = MatSpaceElem{T}

function promote_rule(::Type{S}, ::Type{U}) where {T <: RingElement, S <: Mat{T}, U <: RingElement}
   promote_rule(T, U) == T ? MatSpaceElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::MatSpace{T})() where {T <: RingElement}
   R = base_ring(a)
   entries = Array{T}(undef, a.nrows, a.ncols)
   for i = 1:a.nrows
      for j = 1:a.ncols
         entries[i, j] = zero(R)
      end
   end
   z = MatSpaceElem{T}(entries)
   z.base_ring = R
   return z
end

function (a::MatSpace{T})(b::S) where {S <: RingElement, T <: RingElement}
   R = base_ring(a)
   entries = Array{T}(undef, a.nrows, a.ncols)
   rb = R(b)
   for i = 1:a.nrows
      for j = 1:a.ncols
         if i != j
            entries[i, j] = zero(R)
         else
            entries[i, j] = rb
         end
      end
   end
   z = MatSpaceElem{T}(entries)
   z.base_ring = R
   return z
end

function (a::MatSpace{T})(b::Mat{T}) where {T <: RingElement}
   parent(b) != a && error("Unable to coerce matrix")
   return b
end

function (a::MatSpace{T})(b::AbstractArray{T, 2}) where T <: RingElement
   R = base_ring(a)
   _check_dim(a.nrows, a.ncols, b)
   for i = 1:a.nrows
      for j = 1:a.ncols
         b[i, j] = R(b[i, j])
      end
   end
   z = MatSpaceElem{T}(b)
   z.base_ring = R
   return z
end

function (a::MatSpace{T})(b::AbstractArray{S, 2}) where {S <: RingElement, T <: RingElement}
   R = base_ring(a)
   _check_dim(a.nrows, a.ncols, b)
   entries = Array{T}(undef, a.nrows, a.ncols)
   for i = 1:a.nrows
      for j = 1:a.ncols
         entries[i, j] = R(b[i, j])
      end
   end
   z = MatSpaceElem{T}(entries)
   z.base_ring = R
   return z
end

function (a::MatSpace{T})(b::AbstractArray{S, 1}) where {S <: RingElement, T <: RingElement}
   _check_dim(a.nrows, a.ncols, b)
   b = Array{S, 2}(transpose(reshape(b, a.ncols, a.nrows)))
   z = a(b)
   return z
end

################################################################################
#
#   Matrix constructors
#
################################################################################

@doc Markdown.doc"""
    matrix(R::Ring, arr::AbstractArray{T, 2}) where {T} -> MatElem{T}

> Constructs the matrix over $R$ with entries as in `arr`.
"""
function matrix(R::Ring, arr::AbstractArray{T, 2}) where {T}
   if elem_type(R) === T
      z = MatSpaceElem{elem_type(R)}(arr)
      z.base_ring = R
      return z
   else
      arr_coerce = convert(Array{elem_type(R), 2}, map(R, arr))::Array{elem_type(R), 2}
      return matrix(R, arr_coerce)
   end
end

@doc Markdown.doc"""
    matrix(R::Ring, r::Int, c::Int, arr::AbstractArray{T, 1}) where {T} -> MatElem{T}

> Constructs the $r \times c$ matrix over $R$, where the entries are taken
> row-wise from `arr`.
"""
function matrix(R::Ring, r::Int, c::Int, arr::AbstractArray{T, 1}) where T
   _check_dim(r, c, arr)
   if elem_type(R) === T
     z = MatSpaceElem{elem_type(R)}(r, c, arr)
     z.base_ring = R
     return z
   else
     arr_coerce = convert(Array{elem_type(R), 1}, map(R, arr))::Array{elem_type(R), 1}
     return matrix(R, r, c, arr_coerce)
   end
end

################################################################################
#
#   Zero matrix
#
################################################################################

@doc Markdown.doc"""
    zero_matrix(R::Ring, r::Int, c::Int) -> MatElem

> Return the $r \times c$ zero matrix over $R$.
"""
function zero_matrix(R::Ring, r::Int, c::Int)
   arr = Array{elem_type(R)}(undef, r, c)
   for i in 1:r
      for j in 1:c
         arr[i, j] = zero(R)
      end
   end
   z = MatSpaceElem{elem_type(R)}(arr)
   z.base_ring = R
   return z
end

################################################################################
#
#   Identity matrix
#
################################################################################

@doc Markdown.doc"""
    identity_matrix(R::Ring, n::Int) -> MatElem

> Return the $n \times n$ identity matrix over $R$.
"""
function identity_matrix(R::Ring, n::Int)
   z = zero_matrix(R, n, n)
   for i in 1:n
      z[i, i] = one(R)
   end
   return z
end

###############################################################################
#
#   MatrixSpace constructor
#
###############################################################################

@doc Markdown.doc"""
    MatrixSpace(R::AbstractAlgebra.Ring, r::Int, c::Int, cached::Bool = true)
> Return parent object corresponding to the space of $r\times c$ matrices over
> the ring $R$. If `cached == true` (the default), the returned parent object
> is cached so that it can returned by future calls to the constructor with the
> same dimensions and base ring.
"""
function MatrixSpace(R::AbstractAlgebra.Ring, r::Int, c::Int, cached::Bool = true)
   T = elem_type(R)
   return MatSpace{T}(R, r, c, cached)
end

@doc Markdown.doc"""
    change_base_ring(M::AbstractAlgebra.MatElem, R::AbstractAlgebra.Ring)
> Return the matrix obtained by coercing each entry into `R`.
"""
function change_base_ring(M::AbstractAlgebra.MatElem, R::AbstractAlgebra.Ring)
   r = nrows(M)
   c = ncols(M)
   N = zero_matrix(R, r, c)
   for i in 1:r
      for j in 1:c
         N[i,j] = R(M[i,j])
      end
   end
   return N
end
