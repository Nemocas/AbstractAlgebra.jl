###############################################################################
#
#   Matrix.jl : matrices over rings
#
###############################################################################

export MatrixSpace, add_column, add_column!, add_row, add_row!, can_solve,
       can_solve_left_reduced_triu, can_solve_with_solution, 
       can_solve_with_solution_interpolation, charpoly,
       charpoly_danilevsky!, charpoly_danilevsky_ff!, charpoly_hessenberg!, 
       _check_dim, _checkbounds, dense_matrix_type, det_popov, diagonal_matrix,
       extended_weak_popov, extended_weak_popov_with_transform, fflu!, fflu,
       find_pivot_popov, gram, hessenberg!, hessenberg, hnf, hnf_cohen,
       hnf_cohen_with_transform, hnf_kb, hnf_kb!, hnf_kb_with_transform,
       hnf_minors, hnf_minors_with_transform, hnf_via_popov,
       hnf_via_popov_with_transform, hnf_with_transform, identity_matrix,
       ishessenberg, ishnf, isinvertible, isinvertible_with_inverse, ispopov,
       isrref, issnf, issquare, istriu, isweak_popov, iszero_column,
       iszero_row, kernel, kronecker_product, left_kernel, lu, lu!,
       map_entries, map_entries!, matrix, minpoly, minors, multiply_column,
       multiply_column!, multiply_row, multiply_row!, nrows, ncols, popov,
       popov_with_transform, powers, pseudo_inv, randmat_triu,
       randmat_with_rank, rank, rank_profile_popov, reverse_cols,
       reverse_cols!, reverse_rows, reverse_rows!, right_kernel, rref, rref!,
       rref_rational, rref_rational!, similarity!, snf, snf_with_transform,
       snf_kb, snf_kb!, snf_kb_with_transform, solve, solve_ff, solve_left,
       solve_rational, solve_triu, solve_with_det, swap_cols, swap_cols!,
       swap_rows, swap_rows!, tr, typed_hvcat, typed_hcat, weak_popov,
       weak_popov_with_transform, zero!, zero_matrix

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring(a::MatSpace{T}) where {T <: RingElement} = a.base_ring::parent_type(T)

base_ring(a::MatrixElem{T}) where {T <: RingElement} = a.base_ring::parent_type(T)

function check_parent(a::MatElem, b::MatElem, throw::Bool = true)
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

_checkbounds(i::Int, j::Int) = 1 <= j <= i

function _checkbounds(A, i::Int, j::Int)
  (_checkbounds(nrows(A), i) && _checkbounds(ncols(A), j)) ||
            Base.throw_boundserror(A, (i, j))
end

function _checkbounds(A, rows::AbstractArray{Int}, cols::AbstractArray{Int})
   Base.has_offset_axes(rows, cols) && throw(ArgumentError("offset arrays are not supported"))
   isempty(rows) || _checkbounds(nrows(A), first(rows)) && _checkbounds(nrows(A), last(rows)) ||
      throw(BoundsError(A, rows))
   isempty(cols) || _checkbounds(ncols(A), first(cols)) && _checkbounds(ncols(A), last(cols)) ||
      throw(BoundsError(A, cols))
end

function check_square(A::MatrixElem)
   issquare(A) || throw(DomainError(A, "matrix must be square"))
   A
end

function check_square(S::MatSpace)
   nrows(S) == ncols(S) || throw(DomainError(S, "matrices must be square"))
   S
end

check_square(S::MatAlgebra) = S

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::MatElem, h::UInt)
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
    nrows(a::MatrixElem)

Return the number of rows of the given matrix.
"""
nrows(::MatrixElem)

@doc Markdown.doc"""
    ncols(a::MatrixElem)

Return the number of columns of the given matrix.
"""
ncols(::MatrixElem)

@doc Markdown.doc"""
    length(a::MatrixElem)

Return the number of entries in the given matrix.
"""
length(a::MatrixElem) = nrows(a) * ncols(a)

@doc Markdown.doc"""
    isempty(a::MatrixElem)

Return `true` if `a` does not contain any entry (i.e. `length(a) == 0`), and `false` otherwise.
"""
isempty(a::MatrixElem) = (nrows(a) == 0) | (ncols(a) == 0)

Base.eltype(::Type{<:MatrixElem{T}}) where {T} = T

function Base.isassigned(a::MatrixElem, i, j)
    try
        a[i, j]
        true
    catch e
        if isa(e, BoundsError) || isa(e, UndefRefError)
            return false
        else
            rethrow()
        end
    end
end

@doc Markdown.doc"""
    zero(a::MatSpace)

Return the zero matrix in the given matrix space.
"""
zero(a::MatSpace) = a()

@doc Markdown.doc"""
    zero(x::MatrixElem, R::Ring, r::Int, c::Int)
    zero(x::MatrixElem, R::Ring=base_ring(x))
    zero(x::MatrixElem, r::Int, c::Int)

Return a zero matrix similar to the given matrix, with optionally different
base ring or dimensions.
"""
zero(x::MatrixElem, R::Ring=base_ring(x)) = zero(x, R, nrows(x), ncols(x))
zero(x::MatrixElem, R::Ring, r::Int, c::Int) = zero!(similar(x, R, r, c))
zero(x::MatrixElem, r::Int, c::Int) = zero(x, base_ring(x), r, c)

function zero!(x::MatrixElem)
   R = base_ring(x)
   for i = 1:nrows(x), j = 1:ncols(x)
      x[i, j] = zero(R)
   end
   x
end

@doc Markdown.doc"""
    one(a::MatSpace)

Return the identity matrix of given matrix space. The matrix space must contain
square matrices or else an error is thrown.
"""
one(a::MatSpace) = check_square(a)(1)

@doc Markdown.doc"""
    one(a::MatrixElem)

Return the identity matrix in the same matrix space as $a$. If the space does
not contain square matrices, an error is thrown.
"""
one(a::MatrixElem) = identity_matrix(a)

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

function isone(a::MatrixElem)
   issquare(a) || return false
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

Return `true` if the $i$-th row of the matrix $M$ is zero.
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

Return `true` if the $i$-th column of the matrix $M$ is zero.
"""
function iszero_column(M::MatrixElem{T}, i::Int) where T <: RingElement
  for j in 1:nrows(M)
    if !iszero(M[j, i])
      return false
    end
  end
  return true
end

###############################################################################
#
#   Similar
#
###############################################################################

function _similar(x::MatrixElem{T}, R::Ring, r::Int, c::Int) where T <: RingElement
   TT = elem_type(R)
   M = Matrix{TT}(undef, (r, c))
   z = x isa MatElem ? Generic.MatSpaceElem{TT}(M) : Generic.MatAlgElem{TT}(M)
   z.base_ring = R
   return z
end

similar(x::MatElem, R::Ring, r::Int, c::Int) = _similar(x, R, r, c)

similar(x::MatElem, R::Ring=base_ring(x)) = similar(x, R, nrows(x), ncols(x))

similar(x::MatElem, r::Int, c::Int) = similar(x, base_ring(x), r, c)

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::MatrixElem) = canonical_unit(a[1, 1])

###############################################################################
#
#   getindex
#
###############################################################################

@doc Markdown.doc"""
    Base.getindex(M::MatElem, rows, cols)

When `rows` and `cols` are specified as an `AbstractVector{Int}`, return a copy of
the submatrix $A$ of $M$ defined by `A[i,j] = M[rows[i], cols[j]]`
for `i=1,...,length(rows)` and `j=1,...,length(cols)`.
Instead of a vector, `rows` and `cols` can also be:
* an integer `i`, which is  interpreted as `i:i`, or
* `:`, which is interpreted as `1:nrows(M)` or `1:ncols(M)` respectively.
"""
function getindex(M::MatElem, rows::AbstractVector{Int}, cols::AbstractVector{Int})
   _checkbounds(M, rows, cols)
   A = similar(M, length(rows), length(cols))
   for i in 1:length(rows)
      for j in 1:length(cols)
         A[i, j] = deepcopy(M[rows[i], cols[j]])
      end
   end
   return A
end

getindex(M::MatElem,
         rows::Union{Int,Colon,AbstractVector{Int}},
         cols::Union{Int,Colon,AbstractVector{Int}}) = M[_to_indices(M, rows, cols)...]

function _to_indices(x, rows, cols)
   if rows isa Integer
      rows = rows:rows
   elseif rows isa Colon
      rows = 1:nrows(x)
   end
   if cols isa Integer
      cols = cols:cols
   elseif cols isa Colon
      cols = 1:ncols(x)
   end
   (rows, cols)
end

function Base.view(M::MatElem{T}, rows::Colon, cols::UnitRange{Int}) where T <: RingElement
   return view(M, 1:nrows(M), cols)
end

function Base.view(M::MatElem{T}, rows::UnitRange{Int}, cols::Colon) where T <: RingElement
   return view(M, rows, 1:ncols(M))
end

function Base.view(M::MatElem{T}, rows::Colon, cols::Colon) where T <: RingElement
   return view(M, 1:nrows(M), 1:ncols(M))
end

Base.firstindex(M::MatrixElem, i::Int) = 1

function Base.lastindex(M::MatrixElem, i::Int)
   if i == 1
      return nrows(M)
   elseif i == 2
      return ncols(M)
   else
      error("Dimension in lastindex must be 1 or 2 (got $i)")
   end
end

###############################################################################
#
#   Array interface
#
###############################################################################

Base.ndims(::MatrixElem) = 2

# Cartesian indexing

Base.eachindex(a::MatrixElem) = CartesianIndices((nrows(a), ncols(a)))

Base.@propagate_inbounds Base.getindex(a::MatrixElem, I::CartesianIndex) =
   a[I[1], I[2]]

Base.@propagate_inbounds function Base.setindex!(a::MatrixElem, x, I::CartesianIndex)
   a[I[1], I[2]] = x
   a
end

# iteration

function Base.iterate(a::MatrixElem, ij=(0, 1))
   i, j = ij
   i += 1
   if i > nrows(a)
      iszero(nrows(a)) && return nothing
      i = 1
      j += 1
   end
   j > ncols(a) && return nothing
   a[i, j], (i, j)
end

Base.IteratorSize(::Type{<:MatrixElem}) = Base.HasShape{2}()
Base.IteratorEltype(::Type{<:MatrixElem}) = Base.HasEltype() # default

###############################################################################
#
#   Block replacement
#
###############################################################################

function setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix}, r::UnitRange{Int}, c::UnitRange{Int}) where T
    _checkbounds(a, r, c)
    size(b) == (length(r), length(c)) || throw(DimensionMismatch("tried to assign a $(size(b, 1))x$(size(b, 2)) matrix to a $(length(r))x$(length(c)) destination"))
    startr = first(r)
    startc = first(c)
    for i in r
        for j in c
            a[i, j] = b[i - startr + 1, j - startc + 1]
        end
    end
end

function setindex!(a::MatrixElem{T}, b::Vector, r::UnitRange{Int}, c::UnitRange{Int}) where T
    _checkbounds(a, r, c)
    if !((length(r) == 1 && length(c) == length(b)) || length(c) == 1 && length(r) == length(b))
      throw(DimensionMismatch("tried to assign vector of length $(length(b)) to a $(length(r))x$(length(c)) destination"))
    end
    startr = first(r)
    startc = first(c)
    for i in r
        for j in c
            a[i, j] = b[i - startr + 1 + j - startc]
        end
    end
end

# UnitRange{Int}, Colon
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::UnitRange{Int}, ::Colon) where T = setindex!(a, b, r, 1:ncols(a))

# Colon, UnitRange{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, ::Colon, c::UnitRange{Int}) where T = setindex!(a, b, 1:nrows(a), c)

# Colon, Colon
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, ::Colon, ::Colon) where T = setindex!(a, b, 1:nrows(a), 1:ncols(a))

# Int, UnitRange{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Int, c::UnitRange{Int}) where T = setindex!(a, b, r:r, c)

# UnitRange{Int}, Int
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::UnitRange{Int}, c::Int) where T = setindex!(a, b, r, c:c)

# Int, Colon
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Int, ::Colon) where T = setindex!(a, b, r:r, 1:ncols(a))

# Colon, Int
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, ::Colon, c::Int) where T = setindex!(a, b, 1:nrows(a), c:c)

function _setindex!(a::MatrixElem{T}, b, r, c) where T
   for (i, i2) in enumerate(r)
      for (j, j2) in enumerate(c)
         a[i2, j2] = b[i, j]
      end
   end
end

function _setindex!(a::MatrixElem{T}, b::Vector, r, c) where T
   for (i, i2) in enumerate(r)
      for (j, j2) in enumerate(c)
         a[i2, j2] = b[i + j - 1]
      end
   end
end

# Vector{Int}, Vector{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Vector{Int}, c::Vector{Int}) where T = _setindex!(a, b, r, c)

# Vector{Int}, UnitRange{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Vector{Int}, c::UnitRange{Int}) where T = _setindex!(a, b, r, c)

# UnitRange{Int}, Vector{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::UnitRange{Int}, c::Vector{Int}) where T = _setindex!(a, b, r, c)

# Vector{Int}, Colon
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Vector{Int}, ::Colon) where T = _setindex!(a, b, r, 1:ncols(a))

# Colon, Vector{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, ::Colon, c::Vector{Int}) where T = _setindex!(a, b, 1:nrows(a), c)

# Int, Vector{Int}
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Int, c::Vector{Int}) where T = setindex!(a, b, r:r, c)

# Vector{Int}, Int
setindex!(a::MatrixElem{T}, b::Union{MatrixElem, Matrix, Vector}, r::Vector{Int}, c::Int) where T = setindex!(a, b, r, c:c)

################################################################################
#
#   Size, axes and issquare
#
################################################################################

size(x::MatrixElem) = (nrows(x), ncols(x))

size(t::MatrixElem, d::Integer) = d <= 2 ? size(t)[d] : 1

axes(t::MatrixElem) = Base.OneTo.(size(t))

axes(t::MatrixElem, d::Integer) = Base.OneTo(size(t, d))

###############################################################################
#
#   Matrix spaces iteration
#
###############################################################################

function Base.iterate(M::MatSpace)
   R = base_ring(M)
   p = ProductIterator(fill(R, nrows(M) * ncols(M)); inplace=true)
   a, st = iterate(p) # R is presumably not empty
   M(a), (p, st)
end

function Base.iterate(M::MatSpace, (p, st))
   a_st = iterate(p, st)
   a_st === nothing && return nothing
   M(first(a_st)), (p, last(a_st))
end

Base.eltype(::Type{M}) where {M<:MatSpace} = elem_type(M)
Base.length(M::MatSpace) = BigInt(length(base_ring(M)))^(nrows(M)*ncols(M))

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::MatrixElem; context = nothing)
   r = nrows(a)
   c = ncols(a)
   isempty(a) && return "$r by $c empty matrix"
   mat = Expr(:vcat)
   for i in 1:r
      row = Expr(:row)
      for j in 1:c
         if isassigned(a, i, j)
            push!(row.args, expressify(a[i, j], context = context))
         else
            push!(row.args,Base.undef_ref_str)
         end
      end
      push!(mat.args, row)
   end
   return mat
end

function Base.show(io::IO, a::MatrixElem)
   show_via_expressify(io, a)
end

function Base.show(io::IO, mi::MIME"text/html", a::MatrixElem)
   show_via_expressify(io, mi, a)
end

function Base.show(io::IO, ::MIME"text/plain", a::MatrixElem)
   r = nrows(a)
   c = ncols(a)

   if isempty(a)
      print(io, "$r by $c empty matrix")
      return
   end

   # preprint each element to know the widths so as to align the columns
   strings = String[sprint(print, isassigned(a, i, j) ? a[i, j] : Base.undef_ref_str,
                           context = :compact => true) for i=1:r, j=1:c]
   maxs = maximum(length, strings, dims=1)

   for i = 1:r
      print(io, "[")
      for j = 1:c
         s = strings[i, j]
         s = ' '^(maxs[j] - length(s)) * s
         print(io, s)
         if j != c
            print(io, "   ")
         end
      end
      print(io, "]")
      if i != r
         println(io, "")
      end
   end
end

function show(io::IO, a::MatSpace)
   print(io, "Matrix Space of ")
   print(io, a.nrows, " rows and ", a.ncols, " columns over ")
   print(IOContext(io, :compact => true), base_ring(a))
end

###############################################################################
#
#   Unary operations
#
###############################################################################

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

function *(x::MatElem{T}, y::MatElem{T}) where {T <: RingElement}
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

function *(x::Union{Integer, Rational, AbstractFloat}, y::MatrixElem)
   z = similar(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         z[i, j] = x*y[i, j]
      end
   end
   return z
end

function *(x::T, y::MatrixElem{T}) where {T <: RingElem}
   z = similar(y)
   for i = 1:nrows(y)
      for j = 1:ncols(y)
         z[i, j] = x*y[i, j]
      end
   end
   return z
end

*(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat}) = y*x

*(x::MatrixElem{T}, y::T) where {T <: RingElem} = y*x

@doc Markdown.doc"""
    +(x::Union{Integer, Rational, AbstractFloat}, y::MatrixElem)

Return $S(x) + y$ where $S$ is the parent of $y$.
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
    +(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat})

Return $x + S(y)$ where $S$ is the parent of $x$.
"""
+(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat}) = y + x

@doc Markdown.doc"""
    +(x::T, y::MatrixElem{T}) where {T <: RingElem}

Return $S(x) + y$ where $S$ is the parent of $y$.
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

Return $x + S(y)$ where $S$ is the parent of $x$.
"""
+(x::MatrixElem{T}, y::T) where {T <: RingElem} = y + x

@doc Markdown.doc"""
    -(x::Union{Integer, Rational, AbstractFloat}, y::MatrixElem)

Return $S(x) - y$ where $S$ is the parent of $y$.
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
    -(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat})

Return $x - S(y)$, where $S$ is the parent of $x$.
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
    -(x::T, y::MatrixElem{T}) where {T <: RingElem}

Return $S(x) - y$ where $S$ is the parent of $y$.
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
    -(x::MatrixElem{T}, y::T) where {T <: RingElem}

Return $x - S(y)$, where $S$ is the parent of $a$.
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

Base.literal_pow(::typeof(^), x::T, ::Val{p}) where {p, T <: MatrixElem} = x^p

@doc Markdown.doc"""
    ^(a::MatrixElem, b::Int)

Return $a^b$. We require that the matrix $a$ is square.
"""
function ^(a::MatrixElem, b::Int)
   !issquare(a) && error("Incompatible matrix dimensions in power")
   if b < 0
      return inv(a)^(-b)
   end
   # special case powers of x for constructing polynomials efficiently
   if b == 0
      return identity_matrix(a)
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

###############################################################################
#
#   Comparisons
#
###############################################################################

@doc Markdown.doc"""
    ==(x::MatrixElem{T}, y::MatrixElem{T}) where {T <: RingElement}

Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
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
    isequal(x::MatrixElem{T}, y::MatrixElem{T}) where {T <: RingElement}

Return `true` if $x == y$ exactly, otherwise return `false`. This function is
useful in cases where the entries of the matrices are inexact, e.g. power
series. Only if the power series are precisely the same, to the same precision,
are they declared equal by this function.
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
    ==(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat})

Return `true` if $x == S(y)$ arithmetically, where $S$ is the parent of $x$,
otherwise return `false`.
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
    ==(x::Union{Integer, Rational, AbstractFloat}, y::MatrixElem)

Return `true` if $S(x) == y$ arithmetically, where $S$ is the parent of $y$,
otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::MatrixElem) = y == x

@doc Markdown.doc"""
    ==(x::MatrixElem{T}, y::T) where {T <: RingElem}

Return `true` if $x == S(y)$ arithmetically, where $S$ is the parent of $x$,
otherwise return `false`.
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
    ==(x::T, y::MatrixElem{T}) where {T <: RingElem}

Return `true` if $S(x) == y$ arithmetically, where $S$ is the parent of $y$,
otherwise return `false`.
"""
==(x::T, y::MatrixElem{T}) where {T <: RingElem} = y == x

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::MatrixElem, y::Union{Integer, Rational, AbstractFloat})
   z = similar(x)
   for i = 1:nrows(x)
      for j = 1:ncols(x)
         z[i, j] = divexact(x[i, j], y)
      end
   end
   return z
end

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
#   Symmetry
#
###############################################################################

@doc Markdown.doc"""
    issymmetric(a::MatrixElem)

Return `true` if the given matrix is symmetric with respect to its main
diagonal, otherwise return `false`.
"""
function issymmetric(a::MatrixElem)
    if !issquare(a)
        return false
    end
    for row in 2:nrows(a)
        for col in 1:(row - 1)
            if a[row, col] != a[col, row]
                return false
            end
        end
    end
    return true
end

###############################################################################
#
#   Kronecker product
#
###############################################################################

function kronecker_product(x::MatrixElem{T}, y::MatrixElem{T}) where {T <: RingElement}
    base_ring(parent(x)) == base_ring(parent(y)) || error("Incompatible matrix spaces in matrix operation")
    z = similar(x, nrows(x)*nrows(y), ncols(x)*ncols(y))
    for ix in 1:nrows(x)
       ixr = (ix - 1)*nrows(y)
       for jx in 1:ncols(x)
          jxc = (jx - 1)*ncols(y)
          for iy in 1:nrows(y)
             for jy in 1:ncols(y)
               z[ixr + iy, jxc + jy] = x[ix, jx]*y[iy, jy]
             end
          end
       end
    end
    return z
end

###############################################################################
#
#   Gram matrix
#
###############################################################################

@doc Markdown.doc"""
    gram(x::MatElem)

Return the Gram matrix of $x$, i.e. if $x$ is an $r\times c$ matrix return
the $r\times r$ matrix whose entries $i, j$ are the dot products of the
$i$-th and $j$-th rows, respectively.
"""
function gram(x::MatElem)
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
    tr(x::MatrixElem)

Return the trace of the matrix $a$, i.e. the sum of the diagonal elements. We
require the matrix to be square.
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
    content(x::MatrixElem)

Return the content of the matrix $a$, i.e. the greatest common divisor of all
its entries, assuming it exists.
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
    *(P::perm, x::MatrixElem)

Apply the pemutation $P$ to the rows of the matrix $x$ and return the result.
"""
function *(P::Perm, x::MatrixElem)
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

function lu!(P::Perm, A::MatrixElem{T}) where {T <: FieldElement}
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
    lu(A::MatrixElem{T}, P = SymmetricGroup(nrows(A))) where {T <: FieldElement}

Return a tuple $r, p, L, U$ consisting of the rank of $A$, a permutation
$p$ of $A$ belonging to $P$, a lower triangular matrix $L$ and an upper
triangular matrix $U$ such that $p(A) = LU$, where $p(A)$ stands for the
matrix whose rows are the given permutation $p$ of the rows of $A$.
"""
function lu(A::MatrixElem{T}, P = SymmetricGroup(nrows(A))) where {T <: FieldElement}
   m = nrows(A)
   n = ncols(A)
   P.n != m && error("Permutation does not match matrix")
   p = one(P)
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
   for i = 1:m
      for j = n + 1:m
         L[i, j] = R()
      end
   end
   return rank, p, L, U
end

function fflu!(P::Perm, A::MatrixElem{T}) where {T <: RingElement}
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

function fflu!(P::Perm, A::MatrixElem{T}) where {T <: Union{FieldElement, ResElem}}
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
    fflu(A::MatrixElem{T}, P = SymmetricGroup(nrows(A))) where {T <: RingElement}

Return a tuple $r, d, p, L, U$ consisting of the rank of $A$, a
denominator $d$, a permutation $p$ of $A$ belonging to $P$, a lower
triangular matrix $L$ and an upper triangular matrix $U$ such that
$p(A) = LDU$, where $p(A)$ stands for the matrix whose rows are the
given permutation $p$ of the rows of $A$ and such that $D$ is the diagonal
matrix diag$(p_1, p_1p_2, \ldots, p_{n-2}p_{n-1}, p_{n-1}p_n)$ where the $p_i$
are the inverses of the diagonal entries of $L$. The denominator $d$ is set to
$\pm \mathrm{det}(S)$ where $S$ is an appropriate submatrix of $A$ ($S = A$ if
$A$ is square and nonsingular) and the sign is decided by the parity of the
permutation.
"""
function fflu(A::MatrixElem{T}, P = SymmetricGroup(nrows(A))) where {T <: RingElement}
   m = nrows(A)
   n = ncols(A)
   P.n != m && error("Permutation does not match matrix")
   p = one(P)
   R = base_ring(A)
   U = deepcopy(A)
   L = similar(A, m, m)
   rank, d = fflu!(p, U)
   i = 1
   j = 1
   k = 1
   while i <= m && j <= n
      if !iszero(U[i, j])
         L[i, k] = U[i, j]
         for l = 1:i - 1
            L[l, k] = 0
         end
         for l = i + 1:m
            L[l, k] = U[l, j]
            U[l, j] = 0
         end
         i += 1
         k += 1
      end
      j += 1
   end

   while k <= m
      for l = 1:k - 1
         L[l, k] = 0
      end
      L[k, k] = 1
      for l = k + 1: m
         L[l, k] = 0
      end
      k += 1
   end

   return rank, d, p, L, U
end

###############################################################################
#
#   Reduced row-echelon form
#
###############################################################################

function rref_rational!(A::MatrixElem{T}) where {T <: RingElement}
   m = nrows(A)
   n = ncols(A)
   R = base_ring(A)
   P = one(SymmetricGroup(m))
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
    rref_rational(M::MatrixElem{T}) where {T <: RingElement}

Return a tuple $(r, A, d)$ consisting of the rank $r$ of $M$ and a
denominator $d$ in the base ring of $M$ and a matrix $A$ such that $A/d$ is
the reduced row echelon form of $M$. Note that the denominator is not usually
minimal.
"""
function rref_rational(M::MatrixElem{T}) where {T <: RingElement}
   A = deepcopy(M)
   r, d = rref_rational!(A)
   return r, A, d
end

function rref!(A::MatrixElem{T}) where {T <: FieldElement}
   m = nrows(A)
   n = ncols(A)
   R = base_ring(A)
   P = one(SymmetricGroup(m))
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
    rref(M::MatrixElem{T}) where {T <: FieldElement}

Return a tuple $(r, A)$ consisting of the rank $r$ of $M$ and a reduced row
echelon form $A$ of $M$.
"""
function rref(M::MatrixElem{T}) where {T <: FieldElement}
   A = deepcopy(M)
   r = rref!(A)
   return r, A
end

@doc Markdown.doc"""
    isrref(M::MatrixElem{T}) where {T <: RingElement}

Return `true` if $M$ is in reduced row echelon form, otherwise return
`false`.
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
    isrref(M::MatrixElem{T}) where {T <: FieldElement}

Return `true` if $M$ is in reduced row echelon form, otherwise return
`false`.
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
   S = PolyRing(R)
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
   P = one(SymmetricGroup(n))
   r, d = fflu!(P, A)
   return r < n ? base_ring(M)() : (parity(P) == 0 ? d : -d)
end

@doc Markdown.doc"""
    det(M::MatrixElem{T}) where {T <: FieldElement}

Return the determinant of the matrix $M$. We assume $M$ is square.
"""
function det(M::MatrixElem{T}) where {T <: FieldElement}
   !issquare(M) && error("Not a square matrix in det")
   if nrows(M) == 0
      return one(base_ring(M))
   end
   return det_fflu(M)
end

@doc Markdown.doc"""
    det(M::MatrixElem{T}) where {T <: RingElement}

Return the determinant of the matrix $M$. We assume $M$ is square.
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
   !isdomain_type(elem_type(typeof(coefficient_ring(base_ring(M))))) &&
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
   RR = coefficient_ring(R)
   x = Array{elem_type(RR)}(undef, bound)
   d = Array{elem_type(RR)}(undef, bound)
   X = zero_matrix(RR, n, n)
   b2 = div(bound, 2)
   pt1 = RR(1 - b2)
   for i = 1:bound
      x[i] = RR(i - b2)
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

Return an array consisting of k-combinations of {1,...,n} as arrays.
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
    minors(A::MatElem, k::Int)

Return an array consisting of the k-minors of A
"""
function minors(A::MatElem, k::Int)
   row_indices = combinations(nrows(A), k)
   col_indices = combinations(ncols(A), k)
   mins = Array{elem_type(base_ring(A)), 1}(undef, 0)
   for ri in row_indices
      for ci in col_indices
         push!(mins, det(A[ri, ci]))
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
    rank(M::MatrixElem{T}) where {T <: RingElement}

Return the rank of the matrix $M$.
"""
function rank(M::MatrixElem{T}) where {T <: RingElement}
   n = nrows(M)
   if n == 0
      return 0
   end
   A = deepcopy(M)
   P = one(SymmetricGroup(n))
   r, d = fflu!(P, A)
   return r
end

@doc Markdown.doc"""
    rank(M::MatrixElem{T}) where {T <: FieldElement}

Return the rank of the matrix $M$.
"""
function rank(M::MatrixElem{T}) where {T <: FieldElement}
   n = nrows(M)
   if n == 0
      return 0
   end
   A = deepcopy(M)
   P = one(SymmetricGroup(n))
   return lu!(P, A)
end

###############################################################################
#
#   Linear solving
#
###############################################################################

# Return `flag, y, d` where flag is set to `true` if `Ay = bd` has a solution
# in the base ring. If not, it returns false. The matrix A can be non-square
# and singular. If a solution exists, `y` is set to one such solution. The
# value `d` is set to an appropriate denominator so that `Ay = bd` is a
# solution over the ring. This implements the fraction free decomposition of
# David J. Jeffrey, see "LU Factoring of non-invertible matrices" in ACM
# Communications in Computer Algebra, July 2010. Note that we handle column
# permutations implicitly and add units along the diagonal of the lower
# triangular matrix L instead of removing columns (and corresponding rows of
# the upper triangular matrix U). We also set free variables to zero.
function can_solve_with_solution_fflu(A::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   base_ring(A) != base_ring(b) && error("Base rings don't match in can_solve_with_solution_fflu")
   nrows(A) != nrows(b) && error("Dimensions don't match in can_solve_with_solution_fflu")
   FFLU = deepcopy(A)
   p = one(SymmetricGroup(nrows(A)))
   rank, d = fflu!(p, FFLU)
   flag, y = solve_fflu_precomp(p, FFLU, b)
   n = nrows(A)
   if flag && rank < n
      b2 = p*b
      A2 = p*A
      A3 = A2[rank + 1:n, :]
      flag = A3*y == b2[rank + 1:n, :]*d
   end
   return flag, y, d
end

# Given a fraction free LU decomposition `LdU` of `p(A)` over an integral
# domain with `L` invertible over fraction field, this function will return a
# pair `flag, y` such that `Ay = bd` and `flag` is set to `true` if a solution
# exists. If not, `flag` may be set to `false`, however it is not required to
# be. If `r` is the rank of `A` then the first `r` rows of `p(A)y = p(b)d` will
# hold iff `flag` is `true`. The remaining rows must be checked by the caller.
function solve_fflu_precomp(p::Perm, FFLU::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   x = p * b
   n = nrows(x)
   m = ncols(x)
   R = base_ring(FFLU)
   c = ncols(FFLU)
   t = base_ring(b)()
   s = base_ring(b)()
   minus_one = R(-1)
   y = similar(x, c, m)
   diag = Array{elem_type(R), 1}(undef, n)
   piv = Array{Int, 1}(undef, n)

   rnk = 0
   l = 0
   for i = 1:n
      l += 1
      while l <= c && iszero(FFLU[i, l])
         l += 1
      end
      piv[i] = l
      if l <= c
         diag[i] = FFLU[i, l]
         rnk += 1
      end
   end

   for k in 1:m
      for i in 1:(n - 1)
         t = mul!(t, x[i, k], minus_one)
         for j in (i + 1):n
            if i <= c && piv[i] <= c
               if i == 1
                  x[j, k] = mul_red!(R(), x[j, k], diag[i], false)
               else
                  x[j, k] = mul_red!(x[j, k], x[j, k], diag[i], false)
               end
            else
               x[j, k] = deepcopy(x[j, k])
            end
            if i <= c && piv[i] <= c
               s = mul_red!(s, FFLU[j, piv[i]], t, false)
               x[j, k] = addeq!(x[j, k], s)
            end
            x[j, k] = reduce!(x[j, k])
            if i > 1
               if i <= c && piv[i - 1] <= c
                  flag, x[j, k] = divides(x[j, k], diag[i - 1])
                  if !flag
                     return false, x
                  end
               end
            end
         end
      end

      l = rnk
      for i in c:-1:1
         if l > 0 && i == piv[l]
            if rnk != 0
               y[i, k] = x[l, k]*diag[rnk]
            else
               y[i, k] = deepcopy(x[l, k])
            end
            for j in (piv[l] + 1):c
               t = mul!(t, y[j, k], FFLU[l, j])
               t = mul!(t, t, minus_one)
               y[i, k] = addeq!(y[i, k], t)
            end

            flag, y[i, k] = divides(y[i, k], diag[l])
            if !flag
               return false, x
            end

            l -= 1
         else
            y[i, k] = R()
         end
      end
   end

   return true, y
end

# Return `flag, y` where flag is set to `true` if `Ay = b` has a solution
# in the base field. If not, it returns false. The matrix A can be non-square
# and singular. If a solution exists, `y` is set to one such solution.
# This implements the LU decomposition, for non-invertible matrices, of
# David J. Jeffrey, see "LU Factoring of non-invertible matrices" in ACM
# Communications in Computer Algebra, July 2010. Note that we handle column
# permutations implicitly and add units along the diagonal of the lower
# triangular matrix L instead of removing columns (and corresponding rows of
# the upper triangular matrix U). We also set free variables to zero.
function can_solve_with_solution_lu(A::MatElem{T}, b::MatElem{T}) where {T <: FieldElement}
   base_ring(A) != base_ring(b) && error("Base rings don't match in can_solve_with_solution_lu")
   nrows(A) != nrows(b) && error("Dimensions don't match in can_solve_with_solution_lu")

   if nrows(A) == 0
      return true, zero_matrix(base_ring(A), ncols(A), ncols(b))
   end

   if ncols(A) == 0
      return iszero(b), zero_matrix(base_ring(A), ncols(A), ncols(b))
   end

   LU = deepcopy(A)
   p = one(SymmetricGroup(nrows(A)))
   rank = lu!(p, LU)

   y = solve_lu_precomp(p, LU, b)

   n = nrows(A)
   flag = true
   if rank < n
      b2 = p*b
      A2 = p*A
      A3 = A2[rank + 1:n, :]
      flag = A3*y == b2[rank + 1:n, :]
   end

   return flag, y
end

# Given an LU decomposition `LU` of `p(A)` over a field with `L` invertible,
# this function will return `y` such that the first `r` rows of `p(A)y = p(b)`
# hold, where `r` is the rank of `A`. The remaining rows must be checked by
# the caller.
function solve_lu_precomp(p::Perm, LU::MatElem{T}, b::MatElem{T}) where {T <: FieldElement}
   x = p * b
   n = nrows(x)
   m = ncols(x)
   R = base_ring(LU)
   c = ncols(LU)
   t = base_ring(b)()
   s = base_ring(b)()
   y = similar(x, c, m)

   diag = Array{elem_type(R), 1}(undef, n)
   piv = Array{Int, 1}(undef, n)

   l = 0
   rnk = 0
   for i = 1:n
      l += 1
      while l <= c && iszero(LU[i, l])
         l += 1
      end
      piv[i] = l
      if l <= c
         diag[i] = LU[i, l]
         rnk += 1
      end
   end

   for k in 1:m
      x[1, k] = deepcopy(x[1, k])
      for i in 2:n
         for j in 1:(i - 1)
            # x[i, k] = x[i, k] - LU[i, j] * x[j, k]
            if j <= c
               t = mul_red!(t, -LU[i, j], x[j, k], false)
               if j == 1
                  x[i, k] = x[i, k] + t # LU[i, j] * x[j, k]
               else
                  x[i, k] = addeq!(x[i, k], t)
               end
            else
               x[i, k] = deepcopy(x[i, k])
            end
         end
         x[i, k] = reduce!(x[i, k])
      end

      # Now every entry of x is a proper copy, so we can change the entries
      # as much as we want.

      l = rnk
      for i in c:-1:1
         if l > 0 && i == piv[l]
            y[i, k] = x[l, k]

            for j in (piv[l] + 1):c
               # x[i, k] = x[i, k] - x[j, k] * LU[l, j]
               t = mul_red!(t, y[j, k], -LU[l, j], false)
               y[i, k] = addeq!(y[i, k], t)
            end

            y[i, k] = reduce!(y[i, k])
            y[i, k] = divexact(y[i, k], diag[l])

            l -= 1
         else
            y[i, k] = R()
         end
      end
   end

   return y
end

function solve_ff(M::MatrixElem{T}, b::MatrixElem{T}) where {T <: FieldElement}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   m = nrows(M)
   flag, x, d = can_solve_with_solution_fflu(M, b)
   !flag && error("System not solvable in solve_ff")
   for i in 1:nrows(x)
      for j in 1:ncols(x)
         x[i, j] = divexact(x[i, j], d)
      end
   end
   return x
end

function can_solve_with_solution_with_det(M::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   # We cannot use solve_fflu directly, since it forgot about the (parity of
   # the) permutation.
   R = base_ring(M)
   FFLU = deepcopy(M)
   p = one(SymmetricGroup(nrows(M)))
   rank, d = fflu!(p, FFLU)
   pivots = zeros(Int, nrows(M))
   c = 1
   for r = 1:nrows(M)
      while c <= ncols(M)
         if !iszero(FFLU[r, c])
            pivots[r] = c
            c += 1
            break
         end
         c += 1
      end
   end
   flag, x = solve_fflu_precomp(p, FFLU, b)
   n = nrows(M)
   if flag && rank < n
      b2 = p*b
      A2 = p*M
      A3 = A2[rank + 1:n, :]
      flag = A3*x == b2[rank + 1:n, :]*d
   end
   if !flag
      return false, rank, p, pivots, x, d
   end
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
   return true, rank, p, pivots, x, d
end

function can_solve_with_solution_with_det(M::MatElem{T}, b::MatElem{T}) where {T <: PolyElem}
   flag, r, p, pivots, x, d = can_solve_with_solution_interpolation_inner(M, b)
   return flag, r, p, pivots, x, d
end

# This can be removed once Nemo implements can_solve_with_solution_with_det
# It's here now only because Nemo overloads it
function solve_with_det(M::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   flag, r, p, piv, x, d = can_solve_with_solution_with_det(M, b)
   !flag && error("System not solvable in solve_with_det")
   return x, d
end

function solve_ff(M::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   m = nrows(M)
   n = ncols(M)
   if m == 0
      return zero_matrix(base_ring(M), ncols(M), ncols(b)), base_ring(M)()
   end
   if n == 0
      b != 0 && error("System not soluble in solve_ff")
      return zero_matrix(base_ring(M), ncols(M), ncols(b)), base_ring(M)()
   end
   flag, S, d = can_solve_with_solution_fflu(M, b)
   !flag && error("System not soluble in solve_ff")
   return S, d
end

function can_solve_with_solution_interpolation_inner(M::MatElem{T}, b::MatElem{T}) where {T <: PolyElem}
   m = nrows(M)
   h = ncols(b)
   c = ncols(M)
   R = base_ring(M)
   prm = one(SymmetricGroup(nrows(M)))
   pivots = zeros(Int, nrows(M))
   if m == 0
      return true, 0, prm, pivots, zero_matrix(R, c, h), one(R)
   end
   R = base_ring(M)
   maxlen = 0
   for i = 1:m
      for j = 1:c
         maxlen = max(maxlen, length(M[i, j]))
      end
   end
   if maxlen == 0
      return iszero(b), 0, prm, pivots, zero_matrix(R, c, h), one(R)
   end
   maxlenb = 0
   for i = 1:m
      for j = 1:h
         maxlenb = max(maxlenb, length(b[i, j]))
      end
   end
   # bound from xd = (M*)b where d is the det
   bound = (maxlen - 1)*(max(m, c) - 1) + max(maxlenb, maxlen)
   RR = coefficient_ring(R)
   tmat = matrix(RR, 0, 0, elem_type(RR)[])
   V = Array{typeof(tmat)}(undef, bound)
   d = Array{elem_type(RR)}(undef, bound)
   y = Array{elem_type(RR)}(undef, bound)
   bj = Array{elem_type(RR)}(undef, bound)
   X = similar(tmat, m, c)
   Y = similar(tmat, m, h)
   x = similar(b, c, h)
   b2 = div(bound, 2)
   pt1 = RR(1 - b2)
   l = 1
   i = 1
   pt = 1
   rnk = -1
   firstprm = true
   failues = 0
   while l <= bound
      y[l] = RR(pt - b2)
      # Running out of interpolation points doesn't imply there is no solution
      (y[l] == pt1 && pt != 1) && error("Not enough interpolation points in ring")
      bad_evaluation = false
      for j = 1:m
         for k = 1:c
            X[j, k] = evaluate(M[j, k], y[l])
            if iszero(X[j, k]) && !iszero(M[j, k])
               bad_evaluation = true
               break
            end
         end
         if bad_evaluation
            break
         end
         for k = 1:h
            Y[j, k] = evaluate(b[j, k], y[l])
            if iszero(Y[j, k]) && !iszero(b[j, k])
               bad_evaluation = true
               break
            end
         end
         if bad_evaluation
            break
         end
      end
      try
         if bad_evaluation
            error("Bad evaluation point")
         end
         flag, r, p, pv, Vl, dl = can_solve_with_solution_with_det(X, Y)
         if !flag
            return flag, r, p, pv, zero(x), zero(R)
         end
         p = inv!(p)
	 # Check that new solution has the same pivots as previous ones
         if r != rnk || p != prm || pv != pivots
            if r < rnk # rank is too low: reject
               pt += 1
               continue
            elseif r > rnk # rank has increased: restart
               l = 1
               rnk = r
               prm = p
               pivots = pv
               firstprm = false
            elseif p != prm || pv != pivots # pivot structure different
               reset_prm = false
               for j = 1:length(p.d)
                  # If earlier pivots or row swaps are encountered, restart
                  if firstprm || (p[j] < prm[j] && pv[j] <= pivots[j]) ||
                                 (p[j] == prm[j] && pv[j] < pivots[j] && pv[j] != 0)
                     prm = p
                     pivots = pv
                     l = 1
                     firstprm = false
                     break
                  elseif p[j] > prm[j] || pv[j] > pivots[j] # worse pivots/row swaps: reject
                     reset_prm = true
                     break
                  end
               end
               if reset_prm
                  pt += 1
                  continue
               end
            end
         end
         V[l] = Vl
         d[l] = dl
         y[l] = RR(pt - b2)
         l += 1
      catch e
         if !(e isa ErrorException)
            rethrow(e)
         end
         i = i + 1
      end

      # We tested bound evaluation points and an impossible inverse was
      # encountered for all the values.

      if i > bound && l == 1
         # impossible inverse doesn't imply no solution
         error("Impossible inverse or too many failures in can_solve_with_solution_interpolation")
      end

      pt = pt + 1
   end
   # Interpolate
   for k = 1:h
      for i = 1:c
         for j = 1:bound
            bj[j] = V[j][i, k]
         end
         try
            x[i, k] = interpolate(R, y, bj)
         catch e
            if !(e isa ErrorException)
               rethrow(e)
            end
            return false, rnk, prm, pivots, zero(x), zero(R)
         end
      end
   end
   return true, rnk, prm, pivots, x, interpolate(R, y, d)
end

function can_solve_with_solution_interpolation(M::MatElem{T}, b::MatElem{T}) where {T <: PolyElem}
   flag, r, p, pv, x, d = can_solve_with_solution_interpolation_inner(M, b)
   return flag, x, d
end

@doc Markdown.doc"""
    solve_rational(M::MatElem{T}, b::MatElem{T}) where T <: RingElement

Given a non-singular $n\times n$ matrix over a ring and an $n\times m$
matrix over the same ring, return a tuple $x, d$ consisting of an
$n\times m$ matrix $x$ and a denominator $d$ such that $Ax = db$. The
denominator will be the determinant of $A$ up to sign. If $A$ is singular an
exception is raised.
"""
function solve_rational(M::MatElem{T}, b::MatElem{T}) where T <: RingElement
   return solve_ringelem(M, b)
end

function solve_ringelem(M::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   return solve_ff(M, b)
end

function solve_rational(M::MatElem{T}, b::MatElem{T}) where {T <: PolyElem}
   base_ring(M) != base_ring(b) && error("Base rings don't match in solve")
   nrows(M) != nrows(b) && error("Dimensions don't match in solve")
   flag = true
   try
      flag, x, d = can_solve_with_solution_interpolation(M, b)
      !flag && error("No solution in solve_rational")
      return x, d
   catch e
      if !isa(e, ErrorException)
         rethrow(e)
      end
      !flag && error("No solution in solve_rational")
      return solve_ff(M, b)
   end
end

@doc Markdown.doc"""
    solve(a::MatElem{S}, b::MatElem{S}) where {S <: RingElement}

Given an $m\times r$ matrix $a$ over a ring and an $m\times n$ matrix $b$
over the same ring, return an $r\times n$ matrix $x$ such that $ax = b$. If
no such matrix exists, an exception is raised.
See also [`solve_left`](@ref).
"""
function solve(a::MatElem{S}, b::MatElem{S}
               ) where S <: RingElement
   can, X = can_solve_with_solution(a, b, side = :right)
   can || throw(ArgumentError("Unable to solve linear system"))
   return X
end

@doc Markdown.doc"""
    solve_left(a::MatElem{S}, b::MatElem{S}) where S <: RingElement

Given an $r\times n$ matrix $a$ over a ring and an $m\times n$ matrix $b$
over the same ring, return an $m\times r$ matrix $x$ such that $xa = b$. If
no such matrix exists, an exception is raised.
See also [`solve`](@ref).
"""
function solve_left(a::MatElem{S}, b::MatElem{S}
                    ) where S <: RingElement
   (flag, x) = can_solve_with_solution(a, b; side = :left)
   flag || error("Unable to solve linear system")
   return x
end

# Find the pivot columns of an rref matrix
function find_pivot(A::MatElem{T}) where T <: RingElement
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

###############################################################################
#
#   Upper triangular solving
#
###############################################################################

@doc Markdown.doc"""
    function istriu(A::MatrixElem{T}) where T <: RingElement

Return `true` if $A$ is an upper triangular matrix.
"""
function istriu(A::MatrixElem{T}) where T <: RingElement
   m = nrows(A)
   n = ncols(A)
   d = 0
   for c = 1:n
      for r = m:-1:1
         if !iszero(A[r, c])
            if r < d
               return false
            end
            d = r
            break
         end
      end
   end
   return true
end


@doc Markdown.doc"""
    solve_triu(U::MatElem{T}, b::MatElem{T}, unit::Bool = false) where {T <: FieldElement}

Given a non-singular $n\times n$ matrix over a field which is upper
triangular, and an $n\times m$ matrix over the same field, return an
$n\times m$ matrix $x$ such that $Ax = b$. If $A$ is singular an exception
is raised. If unit is true then $U$ is assumed to have ones on its
diagonal, and the diagonal will not be read.
"""
function solve_triu(U::MatElem{T}, b::MatElem{T}, unit::Bool = false) where {T <: FieldElement}
   n = nrows(U)
   m = ncols(b)
   R = base_ring(U)
   X = zero(b)
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
    can_solve_left_reduced_triu(r::MatElem{T},
                          M::MatElem{T}) where T <: RingElement
Return a tuple `flag, x` where `flag` is set to true if $xM = r$ has a
solution, where $M$ is an $m\times n$ matrix in (upper triangular) Hermite
normal form or reduced row echelon form and $r$ and $x$ are row vectors with
$m$ columns. If there is no solution, flag is set to `false` and $x$ is set
to the zero row.
"""
function can_solve_left_reduced_triu(r::MatElem{T},
                          M::MatElem{T}) where T <: RingElement
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
      x[1, j], r[1, i] = divrem(r[1, i], M[j, k])
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

function can_solve_with_solution(a::MatElem{S}, b::MatElem{S}; side::Symbol = :right) where S <: FracElem{T} where T <: PolyElem
   if side == :left
      (f, x) = can_solve_with_solution(a', b'; side=:right)
      return (f, x')
   elseif side == :right
      d = numerator(one(base_ring(a)))
      for i = 1:nrows(a)
         for j = 1:ncols(a)
            d = lcm(d, denominator(a[i, j]))
         end
      end
      for i = 1:nrows(b)
         for j = 1:ncols(b)
            d = lcm(d, denominator(b[i, j]))
         end
      end
      A = matrix(parent(d), nrows(a), ncols(a), [numerator(a[i, j]*d) for i in 1:nrows(a) for j in 1:ncols(a)])
      B = matrix(parent(d), nrows(b), ncols(b), [numerator(b[i, j]*d) for i in 1:nrows(b) for j in 1:ncols(b)])
      flag = false
      x = similar(A, ncols(A), nrows(B))
      den = one(base_ring(a))
      try
         flag, x, den = can_solve_with_solution_interpolation(A, B)
      catch
         flag, x, den = can_solve_with_solution_fflu(A, B)
      end
      X = change_base_ring(base_ring(a), x)
      X = divexact(X, base_ring(a)(den))
      return flag, X
   else
      error("Unsupported argument :$side for side: Must be :left or :right.")
   end
end

# The fflu approach is the fastest over a fraction field (see benchmarks on PR 661)
function can_solve_with_solution(a::MatElem{S}, b::MatElem{S}; side::Symbol = :right) where S <: Union{FracElem, Rational{BigInt}}
   if side == :left
      (f, x) = can_solve_with_solution(a', b'; side=:right)
      return (f, x')
   elseif side == :right
      d = numerator(one(base_ring(a)))
      for i = 1:nrows(a)
         for j = 1:ncols(a)
            d = lcm(d, denominator(a[i, j]))
         end
      end
      for i = 1:nrows(b)
         for j = 1:ncols(b)
            d = lcm(d, denominator(b[i, j]))
         end
      end
      A = matrix(parent(d), nrows(a), ncols(a), [numerator(a[i, j]*d) for i in 1:nrows(a) for j in 1:ncols(a)])
      B = matrix(parent(d), nrows(b), ncols(b), [numerator(b[i, j]*d) for i in 1:nrows(b) for j in 1:ncols(b)])
      flag, x, den = can_solve_with_solution_fflu(A, B)
      X = change_base_ring(base_ring(a), x)
      X = divexact(X, base_ring(a)(den))
      return flag, X
   else
      error("Unsupported argument :$side for side: Must be :left or :right.")
   end
end

@doc Markdown.doc"""
    can_solve_with_solution(a::MatElem{S}, b::MatElem{S}; side::Symbol = :right) where S <: RingElement

Given two matrices $a$ and $b$ over the same ring, try to solve $ax = b$
if `side` is `:right` or $xa = b$ if `side` is `:left`. In either case,
return a tuple `(flag, x)`. If a solution exists, `flag` is set to true and
`x` is a solution. If no solution exists, `flag` is set to false and `x`
is arbitrary. If the dimensions of $a$ and $b$ are incompatible, an exception
is raised.
"""
function can_solve_with_solution(a::MatElem{S}, b::MatElem{S}; side::Symbol = :right) where S <: RingElement
   if side == :right
      (f, x) = can_solve_with_solution(a', b'; side=:left)
      return (f, x')
   elseif side == :left
      @assert ncols(a) == ncols(b)
      H, T = hnf_with_transform(a)
      b = deepcopy(b)
      z = zero(a, nrows(b), nrows(a))
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
            q, r = divrem(b[i, k], H[j, k])
            if r != 0
               return (false, zero(a, 0, 0))
            end
            z[i, j] = q
            q = -q
            for h = k:ncols(H)
               t = mul!(t, q, H[j, h])
               b[i, h] = addeq!(b[i, h], t)
            end
         end
      end
      if b != 0
         return (false, zero(a, 0, 0))
      end
      return (true, z*T)
   else
      error("Unsupported argument :$side for side: Must be :left or :right.")
   end
end

function can_solve_with_solution(A::MatElem{T}, B::MatElem{T}; side::Symbol = :right) where T <: FieldElement
   if side == :right
      (f, x) = can_solve_with_solution(A', B', side = :left)
      return (f, x')
   elseif side == :left
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
         return (false, zero(A, 0, 0))
      end
      sol = zero_matrix(R, nrows(B), nrows(A))
      for i = 1:length(p)
         for j = 1:nrows(B)
            sol[j, p[i]] = mu[i, nrows(A) + j]
         end
      end
      return (true, sol)
   else
      error("Unsupported argument :$side for side: Must be :left or :right.")
   end
end

@doc Markdown.doc"""
    can_solve(a::MatElem{S}, b::MatElem{S}; side::Symbol = :right) where S <: RingElement

Given two matrices $a$ and $b$ over the same ring, check the solubility
of $ax = b$ if `side` is `:right` or $xa = b$ if `side` is `:left`.
Return true if a solution exists, false otherwise. If the dimensions
of $a$ and $b$ are incompatible, an exception is raised. If a solution
should be computed as well, use `can_solve_with_solution` instead.
"""
function can_solve(a::MatElem{S}, b::MatElem{S}; side::Symbol = :right) where S <: RingElement
   return can_solve_with_solution(a, b; side=side)[1]
end

###############################################################################
#
#   Inverse
#
###############################################################################

@doc Markdown.doc"""
    pseudo_inv(M::MatrixElem{T}) where {T <: RingElement}

Given a non-singular $n\times n$ matrix $M$ over a ring return a tuple $X, d$
consisting of an $n\times n$ matrix $X$ and a denominator $d$ such that
$MX = dI_n$, where $I_n$ is the $n\times n$ identity matrix. The denominator
will be the determinant of $M$ up to sign. If $M$ is singular an exception
is raised.
"""
function pseudo_inv(M::MatrixElem{T}) where {T <: RingElement}
   issquare(M) || throw(DomainError(M, "Can not invert non-square Matrix"))
   flag, X, d = can_solve_with_solution_fflu(M, identity_matrix(M))
   !flag && error("Singular matrix in pseudo_inv")
   return X, d
end

function Base.inv(M::MatrixElem{T}) where {T <: FieldElement}
   issquare(M) || throw(DomainError(M, "Can not invert non-square Matrix"))
   flag, A = can_solve_with_solution_lu(M, identity_matrix(M))
   !flag && error("Singular matrix in inv")
   return A
end

@doc Markdown.doc"""
    inv(M::MatrixElem{T}) where {T <: RingElement}

Given a non-singular $n\times n$ matrix over a ring, return an
$n\times n$ matrix $X$ such that $MX = I_n$, where $I_n$ is the $n\times n$
identity matrix. If $M$ is not invertible over the base ring an exception is
raised.
"""
function Base.inv(M::MatrixElem{T}) where {T <: RingElement}
   issquare(M) || throw(DomainError(M, "Can not invert non-square Matrix"))
   X, d = pseudo_inv(M)
   isunit(d) || throw(DomainError(M, "Matrix is not invertible."))
   return divexact(X, d)
end

###############################################################################
#
#   Is invertible
#
###############################################################################

@doc Markdown.doc"""
    isinvertible_with_inverse(A::MatrixElem{T}; side::Symbol = :left) where {T <: RingElement}

Given an $n\times m$ matrix $A$ over a ring, return a tuple `(flag, B)`.
If `side` is `:right` and `flag` is true, $B$ is the right inverse of $A$
i.e. $AB$ is the $n\times n$ unit matrix. If `side` is `:left` and `flag` is
true, $B$ is the left inverse of $A$ i.e. $BA$ is the $m\times m$ unit matrix.
If `flag` is false, no right or left inverse exists.
"""
function isinvertible_with_inverse(A::MatrixElem{T}; side::Symbol = :left) where {T <: RingElement}
   if (side == :left && nrows(A) < ncols(A)) || (side == :right && ncols(A) < nrows(A))
      return (false, zero(A, 0, 0))
   end
   I = (side == :left) ? zero(A, ncols(A), ncols(A)) : zero(A, nrows(A), nrows(A))
   for i = 1:ncols(I)
      I[i, i] = one(base_ring(I))
   end
   return can_solve_with_solution(A, I; side = side)
end

@doc Markdown.doc"""
    isinvertible(A::MatrixElem{T}) where {T <: RingElement}

Return true if a given square matrix is invertible, false otherwise. If
the inverse should also be computed, use `isinvertible_with_inverse`.
"""
isinvertible(A::MatrixElem{T}) where {T <: RingElement} = issquare(A) && isunit(det(A))

isinvertible(A::MatrixElem{T}) where {T <: FieldElement} = nrows(A) == ncols(A) == rank(A)

###############################################################################
#
#   Nullspace
#
###############################################################################

@doc Markdown.doc"""
    nullspace(M::MatElem{T}) where {T <: RingElement}

Return a tuple $(\nu, N)$ consisting of the nullity $\nu$ of $M$ and
a basis $N$ (consisting of column vectors) for the right nullspace of $M$,
i.e. such that $MN$ is the zero matrix. If $M$ is an $m\times n$ matrix
$N$ will be an $n\times \nu$ matrix. Note that the nullspace is taken to be
the vector space kernel over the fraction field of the base ring if the
latter is not a field. In AbstractAlgebra we use the name "kernel" for a
function to compute an integral kernel.
"""
function nullspace(M::MatElem{T}) where {T <: RingElement}
   n = ncols(M)
   rank, A, d = rref_rational(M)
   nullity = n - rank
   R = base_ring(M)
   U = zero(M, n, nullity)
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
    nullspace(M::MatElem{T}) where {T <: FieldElement}

Return a tuple $(\nu, N)$ consisting of the nullity $\nu$ of $M$ and
a basis $N$ (consisting of column vectors) for the right nullspace of $M$,
i.e. such that $MN$ is the zero matrix. If $M$ is an $m\times n$ matrix
$N$ will be an $n\times \nu$ matrix.
"""
function nullspace(M::MatElem{T}) where {T <: FieldElement}
   m = nrows(M)
   n = ncols(M)
   rank, A = rref(M)
   nullity = n - rank
   R = base_ring(M)
   X = zero(M, n, nullity)
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
    left_kernel(a::MatElem{T}) where T <: RingElement

Return a tuple `n, M` where $M$ is a matrix whose rows generate the kernel
of $M$ and $n$ is the rank of the kernel. The transpose of the output of this
function is guaranteed to be in flipped upper triangular format (i.e. upper
triangular format if columns and rows are reversed).
"""
function left_kernel(x::MatElem{T}) where T <: RingElement
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

function left_kernel(M::MatElem{T}) where T <: FieldElement
  n, N = nullspace(transpose(M))
  return n, transpose(N)
end

@doc Markdown.doc"""
    right_kernel(a::MatElem{T}) where T <: RingElement

Return a tuple `n, M` where $M$ is a matrix whose columns generate the
kernel of $M$ and $n$ is the rank of the kernel.
"""
function right_kernel(x::MatElem{T}) where T <: RingElement
   n, M = left_kernel(transpose(x))
   return n, transpose(M)
end

function right_kernel(M::MatElem{T}) where T <: FieldElement
   return nullspace(M)
end

@doc Markdown.doc"""
    kernel(a::MatElem{T}; side::Symbol = :right) where T <: RingElement

Return a tuple $(n, M)$, where n is the rank of the kernel and $M$ is a
basis for it. If side is $:right$ or not specified, the right kernel is
computed, i.e. the matrix of columns whose span gives the right kernel
space. If side is $:left$, the left kernel is computed, i.e. the matrix
of rows whose span is the left kernel space.
"""
function kernel(A::MatElem{T}; side::Symbol = :right) where T <: RingElement
   if side == :right
      return right_kernel(A)
   elseif side == :left
      return left_kernel(A)
   else
      error("Unsupported argument: :$side for side: must be :left or :right")
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
    hessenberg(A::MatrixElem{T}) where {T <: RingElement}

Return the Hessenberg form of $M$, i.e. an upper Hessenberg matrix
which is similar to $M$. The upper Hessenberg form has nonzero entries
above and on the diagonal and in the diagonal line immediately below the
diagonal.
"""
function hessenberg(A::MatrixElem{T}) where {T <: RingElement}
   !issquare(A) && error("Dimensions don't match in hessenberg")
   M = deepcopy(A)
   hessenberg!(M)
   return M
end

@doc Markdown.doc"""
    ishessenberg(A::MatrixElem{T}) where {T <: RingElement}

Return `true` if $M$ is in Hessenberg form, otherwise returns `false`.
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
   coefficient_ring(S) != base_ring(A) && error("Cannot coerce into polynomial ring")
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
   coefficient_ring(S) != base_ring(A) && error("Cannot coerce into polynomial ring")
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
   coefficient_ring(S) != base_ring(A) && error("Cannot coerce into polynomial ring")
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
    charpoly(V::Ring, Y::MatrixElem{T}) where {T <: RingElement}

Return the characteristic polynomial $p$ of the matrix $M$. The
polynomial ring $R$ of the resulting polynomial must be supplied
and the matrix is assumed to be square.
"""
function charpoly(V::Ring, Y::MatrixElem{T}) where {T <: RingElement}
   !issquare(Y) && error("Dimensions don't match in charpoly")
   R = base_ring(Y)
   coefficient_ring(V) != base_ring(Y) && error("Cannot coerce into polynomial ring")
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

Return the minimal polynomial $p$ of the matrix $M$. The polynomial ring $S$
of the resulting polynomial must be supplied and the matrix must be square.
"""
function minpoly(S::Ring, M::MatElem{T}, charpoly_only::Bool = false) where {T <: FieldElement}
   !issquare(M) && error("Not a square matrix in minpoly")
   coefficient_ring(S) != base_ring(M) && error("Unable to coerce polynomial")
   n = nrows(M)
   if n == 0
      return S(1)
   end
   R = base_ring(M)
   p = S(1)
   A = similar(M, n + 1, 2n + 1)
   B = similar(M, n, n)
   L1 = Int[n + i for i in 1:n + 1]
   L2 = Int[n for i in 1:n]
   P2 = zeros(Int, n)
   P2[1] = 1
   c2 = 1
   r2 = 1
   first_poly = true
   while r2 <= n
      P1 = Int[0 for i in 1:2n + 1]
      v = zero(M, n, 1)
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

Return the minimal polynomial $p$ of the matrix $M$. The polynomial ring $S$
of the resulting polynomial must be supplied and the matrix must be square.
"""
function minpoly(S::Ring, M::MatElem{T}, charpoly_only::Bool = false) where {T <: RingElement}
   !issquare(M) && error("Not a square matrix in minpoly")
   coefficient_ring(S) != base_ring(M) && error("Unable to coerce polynomial")
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
      v = zero(M, n, 1)
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
   U = identity_matrix(A, m)
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
      for j = k + 1:m
         if iszero(H[j, i])
            continue
         end
         d, u, v = gcdx(H[k, i], H[j, i])
         a = divexact(H[k, i], d)
         b = -divexact(H[j, i], d)
         for c = i:n
            t = deepcopy(H[j, c])
            t1 = mul_red!(t1, a, H[j, c], false)
            t2 = mul_red!(t2, b, H[k, c], false)
            H[j, c] = t1 + t2
            H[j, c] = reduce!(H[j, c])
            t1 = mul_red!(t1, u, H[k, c], false)
            t2 = mul_red!(t2, v, t, false)
            H[k, c] = t1 + t2
            H[k, c] = reduce!(H[k, c])
         end
         for c = 1:m
            t = deepcopy(U[j,c])
            t1 = mul_red!(t1, a, U[j, c], false)
            t2 = mul_red!(t2, b, U[k, c], false)
            U[j, c] = t1 + t2
            U[j, c] = reduce!(U[j, c])
            t1 = mul_red!(t1, u, U[k, c], false)
            t2 = mul_red!(t2, v, t, false)
            U[k, c] = t1 + t2
            U[k, c] = reduce!(U[k, c])
         end
      end
      if iszero(H[k, i])
         continue
      end
      cu = canonical_unit(H[k, i])
      if cu != 1
         for c = i:n
            H[k, c] = divexact(H[k, c], cu)
        end
         for c = 1:m
            U[k, c] = divexact(U[k, c], cu)
         end
      end
      for j = 1:k-1
         q = -div(H[j,i], H[k, i])
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
    hnf_minors(A::MatrixElem{T}) where {T <: RingElement}

Compute the upper right row Hermite normal form of $A$ using the algorithm of
Kannan-Bachem. The input must have full column rank.
"""
function hnf_minors(A::MatrixElem{T}) where {T <: RingElement}
   H = deepcopy(A)
   _hnf_minors!(H, similar(A, 0, 0), Val{false})
   return H
end

@doc Markdown.doc"""
    hnf_minors_with_transform(A::MatrixElem{T}) where {T <: RingElement}

Compute the upper right row Hermite normal form $H$ of $A$ and an invertible
matrix $U$ with $UA = H$ using the algorithm of Kannan-Bachem. The input must
have full column rank.
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
            q = div(H[i, j], H[j, j])
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
            q = div(H[i, j], H[j, j])
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
    hnf_kb(A::MatrixElem{T}) where {T <: RingElement}

Compute the upper right row Hermite normal form of $A$ using a modification
of the algorithm of Kannan-Bachem.
"""
function hnf_kb(A::MatrixElem{T}) where {T <: RingElement}
   return _hnf_kb(A, Val{false})
end

@doc Markdown.doc"""
    hnf_kb_with_transform(A::MatrixElem{T}) where {T <: RingElement}

Compute the upper right row Hermite normal form $H$ of $A$ and an invertible
matrix $U$ with $UA = H$ using a modification of the algorithm of
Kannan-Bachem.
"""
function hnf_kb_with_transform(A::MatrixElem{T}) where {T <: RingElement}
   return _hnf_kb(A, Val{true})
end

function _hnf_kb(A, trafo::Type{Val{T}} = Val{false}) where T
   H = deepcopy(A)
   m = nrows(H)
   if trafo == Val{true}
      U = identity_matrix(A, m)
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

# Reduces the entries above H[pivot[c], c]
function kb_reduce_column!(H::MatrixElem{T}, U::MatrixElem{T}, pivot::Array{Int, 1}, c::Int, with_trafo::Bool, start_element::Int = 1) where {T <: RingElement}

   # Let c = 4 and pivot[c] = 4. H could look like this:
   # ( 0 . * # * )
   # ( . * * # * )
   # ( 0 0 0 0 . )
   # ( 0 0 0 . * )
   # ( * * * * * )
   #
   # (. are pivots, we want to reduce the entries marked with #)
   # The #'s are in rows whose pivot is in a column left of column c.

   r = pivot[c]
   t = base_ring(H)()
   for i = start_element:c - 1
      p = pivot[i]
      if p == 0
         continue
      end
      # So, the pivot in row p is in a column left of c.
      if iszero(H[p, c])
         continue
      end
      q = -div(H[p, c], H[r, c])
      for j = c:ncols(H)
         t = mul!(t, q, H[r, j])
         H[p, j] += t
      end
      if with_trafo
         for j = 1:ncols(U)
            t = mul!(t, q, U[r, j])
            U[p, j] += t
         end
      end
   end
   return nothing
end

# Multiplies row r by a unit such that the entry H[r, c] is "canonical"
function kb_canonical_row!(H, U, r::Int, c::Int, with_trafo::Bool)
   cu = canonical_unit(H[r, c])
   if cu != 1
      for j = c:ncols(H)
         H[r, j] = divexact(H[r, j], cu)
      end
      if with_trafo
         for j = 1:ncols(U)
            U[r, j] = divexact(U[r, j], cu)
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
   pivot = zeros(Int, n) # pivot[j] == i if the pivot of column j is in row i

   # Find the first non-zero entry of H
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
   for i = row1 + 1:m
      new_pivot = false
      for j = start_element:n
         if iszero(H[i, j])
            continue
         end
         if pivot[j] == 0
            # We found a non-zero entry in a column without a pivot: This is a
            # new pivot
            pivot[j] = i
            pivot_max = max(pivot_max, j)
            new_pivot = true
         else
            # We have a pivot for this column: Use it to write 0 in H[i, j]
            p = pivot[j]
            d, u, v = gcdx(H[p, j], H[i, j])
            a = divexact(H[p, j], d)
            b = -divexact(H[i, j], d)
            for c = j:n
               t = deepcopy(H[i, c])
               t1 = mul_red!(t1, a, H[i, c], false)
               t2 = mul_red!(t2, b, H[p, c], false)
               H[i, c] = reduce!(t1 + t2)
               t1 = mul_red!(t1, u, H[p, c], false)
               t2 = mul_red!(t2, v, t, false)
               H[p, c] = reduce!(t1 + t2)
            end
            if with_trafo
               for c = 1:m
                  t = deepcopy(U[i, c])
                  t1 = mul_red!(t1, a, U[i, c], false)
                  t2 = mul_red!(t2, b, U[p, c], false)
                  U[i, c] = reduce!(t1 + t2)
                  t1 = mul_red!(t1, u, U[p, c], false)
                  t2 = mul_red!(t2, v, t, false)
                  U[p, c] = reduce!(t1 + t2)
               end
            end
         end

         # We changed the pivot of column j (or found a new one).
         # We have do reduce the entries marked with # in
         # ( 0 0 0 . * )
         # ( . # # * * )
         # ( 0 0 . * * )
         # ( 0 . # * * )
         # ( * * * * * )
         # where . are pivots and i = 4, j = 2. (This example is for the
         # "new pivot" case.)
         kb_canonical_row!(H, U, pivot[j], j, with_trafo)
         for c = j:pivot_max
            if pivot[c] == 0
               continue
            end
            kb_reduce_column!(H, U, pivot, c, with_trafo, start_element)
         end
         if new_pivot
            break
         end
      end
   end
   kb_sort_rows!(H, U, pivot, with_trafo, start_element)
   return nothing
end

@doc Markdown.doc"""
    hnf(A::MatrixElem{T}) where {T <: RingElement}

Return the upper right row Hermite normal form of $A$.
"""
function hnf(A::MatrixElem{T}) where {T <: RingElement}
  return hnf_kb(A)
end

@doc Markdown.doc"""
    hnf_with_transform(A)

Return the tuple $H, U$ consisting of the upper right row Hermite normal
form $H$ of $A$ together with invertible matrix $U$ such that $UA = H$.
"""
function hnf_with_transform(A)
  return hnf_kb_with_transform(A)
end

@doc Markdown.doc"""
    ishnf(M::MatrixElem{T}) where T <: RingElement

Return `true` if the matrix is in Hermite normal form.
"""
function ishnf(M::MatrixElem{T}) where T <: RingElement
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
         qq, rr = divrem(M[i, col], p)
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

@doc Markdown.doc"""
    issnf(A::MatrixElem{T}) where T <: RingElement

Return `true` if $A$ is in Smith Normal Form.
"""
function issnf(A::MatrixElem{T}) where T <: RingElement
   m = nrows(A)
   n = ncols(A)
   a = A[1, 1]
   for i = 2:min(m, n)
      q, r = divrem(A[i, i], a)
      if !iszero(r)
         return false
      end
      a = A[i,i]
   end
   for i = 1:n
      for j = 1:m
         if i == j
            continue
         end
         if !iszero(A[j, i])
            return false
         end
      end
   end
   return true
end

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
      U = identity_matrix(A, m)
      K = identity_matrix(A, n)
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
         S[r, j] = reduce!(t1 + t2)
         t1 = mul_red!(t1, u, S[r, i], false)
         t2 = mul_red!(t2, v, t, false)
         S[r, i] = reduce!(t1 + t2)
      end
      if with_trafo
         for r = 1:n
            t = deepcopy(K[r,j])
            t1 = mul_red!(t1, a, K[r, j], false)
            t2 = mul_red!(t2, b, K[r, i], false)
            K[r, j] = reduce!(t1 + t2)
            t1 = mul_red!(t1, u, K[r, i], false)
            t2 = mul_red!(t2, v, t, false)
            K[r, i] = reduce!(t1 + t2)
         end
      end
   end
   return nothing
end

function snf_kb!(S::MatrixElem{T}, U::MatrixElem{T}, K::MatrixElem{T}, with_trafo::Bool = false) where {T <: RingElement}
   m = nrows(S)
   n = ncols(S)
   l = min(m, n)
   i = 1
   t = base_ring(S)()
   t1 = base_ring(S)()
   t2 = base_ring(S)()
   while i <= l
      kb_clear_row!(S, K, i, with_trafo)
      hnf_kb!(S, U, with_trafo, i)
      c = i + 1
      while c <= n && iszero(S[i, c])
         c += 1
      end
      if c != n + 1
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
               t = deepcopy(U[i, c])
               U[i, c] += U[j, c]
               t2 = mul_red!(t2, t1, U[j, c], false)
               U[j, c] += t2
               t2 = mul_red!(t2, t1, t, false)
               U[j, c] = reduce!(U[j, c] + t2)
            end
            q1 = -divexact(S[j, j], d)
            q2 = divexact(S[i, i], d)
            for r = 1:n
               t = deepcopy(K[r, i])
               t1 = mul_red!(t1, K[r, i], u, false)
               t2 = mul_red!(t2, K[r, j], v, false)
               K[r, i] = reduce!(t1 + t2)
               t1 = mul_red!(t1, t, q1, false)
               t2 = mul_red!(t2, K[r, j], q2, false)
               K[r, j] = reduce!(t1 + t2)
            end
         end
         S[j, j] = divexact(S[i, i]*S[j, j], d)
         S[i, i] = d
      end
   end
   return nothing
end

@doc Markdown.doc"""
    snf(A::MatrixElem{T}) where {T <: RingElement}

Return the Smith normal form of $A$.
"""
function snf(a::MatrixElem{T}) where {T <: RingElement}
  return snf_kb(a)
end

@doc Markdown.doc"""
    snf_with_transform(A)

Return the tuple $S, T, U$ consisting of the Smith normal form $S$ of $A$
together with invertible matrices $T$ and $U$ such that $TAU = S$.
"""
function snf_with_transform(a::MatrixElem{T}) where {T <: RingElement}
  return snf_kb_with_transform(a)
end

################################################################################
#
#   Popov Form
#
################################################################################

@doc Markdown.doc"""
    isweak_popov(P::MatrixElem{T}, rank::Int) where T <: PolyElem

Return `true` if $P$ is a matrix in weak Popov form of the given rank.
"""
function isweak_popov(P::MatrixElem{T}, rank::Int) where T <: PolyElem
   zero_rows = 0
   pivots = zeros(ncols(P))
   for r = 1:nrows(P)
      p = find_pivot_popov(P, r)
      if P[r, p] == 0
         zero_rows += 1
         continue
      end
      # There is already a pivot in this column
      if pivots[p] != 0
         return false
      end
      pivots[p] = r
   end
   if zero_rows != nrows(P) - rank
      return false
   end
   return true
end

@doc Markdown.doc"""
    ispopov(P::MatrixElem{T}, rank::Int) where T <: PolyElem

Return `true` if $P$ is a matrix in Popov form with the given rank.
"""
function ispopov(P::MatrixElem{T}, rank::Int) where T <: PolyElem
   zero_rows = 0
   for r = 1:nrows(P)
      p = find_pivot_popov(P, r)
      if P[r, p] != 0
         break
      end
      zero_rows += 1
   end
   # The zero rows must all be on top.
   if zero_rows != nrows(P) - rank
      return false
   end
   pivotscr = zeros(Int, ncols(P)) # pivotscr[i] == j means the pivot of column i is in row j
   pivots = zeros(Int, nrows(P)) # the other way round
   for r = zero_rows + 1:nrows(P)
      p = find_pivot_popov(P, r)
      if P[r, p] == 0
         return false
      end
      if pivotscr[p] != 0
         # There is already a pivot in this column
         return false
      end
      pivotscr[p] = r
      pivots[r] = p
   end
   for r = zero_rows + 1:nrows(P)
      p = pivots[r]
      f = P[r, p]
      if !isone(leading_coefficient(f))
         return false
      end
      for i = 1:nrows(P)
         i == r ? continue : nothing
         if degree(P[i, p]) >= degree(f)
            return false
         end
      end
      if r == nrows(P)
         break
      end
      g = P[r + 1, pivots[r + 1]]
      if degree(f) >= degree(g)
         if pivots[r] >= pivots[r + 1]
            return false
         end
      end
   end
   return true
end

@doc Markdown.doc"""
    weak_popov(A::MatElem{T}) where {T <: PolyElem}

Return the weak Popov form of $A$.
"""
function weak_popov(A::MatElem{T}) where {T <: PolyElem}
   return _weak_popov(A, Val{false})
end

@doc Markdown.doc"""
    weak_popov_with_transform(A::MatElem{T}) where {T <: PolyElem}

Compute a tuple $(P, U)$ where $P$ is the weak Popov form of $A$ and $U$
is a transformation matrix so that $P = UA$.
"""
function weak_popov_with_transform(A::MatElem{T}) where {T <: PolyElem}
   return _weak_popov(A, Val{true})
end

function _weak_popov(A::MatElem{T}, trafo::Type{Val{S}} = Val{false}) where {T <: PolyElem, S}
   P = deepcopy(A)
   m = nrows(P)
   W = similar(A, 0, 0)
   if trafo == Val{true}
      U = identity_matrix(A, m)
      weak_popov!(P, W, U, false, true)
      return P, U
   else
      U = similar(A, 0, 0)
      weak_popov!(P, W, U, false, false)
      return P
   end
end

@doc Markdown.doc"""
    extended_weak_popov(A::MatElem{T}, V::MatElem{T}) where {T <: PolyElem}

Compute the weak Popov form $P$ of $A$ by applying simple row transformations
on $A$ and a vector $W$ by applying the same transformations on the vector $V$.
Return the tuple $(P, W)$.
"""
function extended_weak_popov(A::MatElem{T}, V::MatElem{T}) where {T <: PolyElem}
   return _extended_weak_popov(A, V, Val{false})
end

@doc Markdown.doc"""
    extended_weak_popov_with_transform(A::MatElem{T}, V::MatElem{T}) where {T <: PolyElem}

Compute the weak Popov form $P$ of $A$ by applying simple row transformations
on $A$, a vector $W$ by applying the same transformations on the vector $V$,
and a transformation matrix $U$ so that $P = UA$.
Return the tuple $(P, W, U)$.
"""
function extended_weak_popov_with_transform(A::MatElem{T}, V::MatElem{T}) where {T <: PolyElem}
   return _extended_weak_popov(A, V, Val{true})
end

function _extended_weak_popov(A::MatElem{T}, V::MatElem{T}, trafo::Type{Val{S}} = Val{false}) where {T <: PolyElem, S}
   @assert nrows(V) == nrows(A) && ncols(V) == 1
   P = deepcopy(A)
   W = deepcopy(V)
   m = nrows(P)
   if trafo == Val{true}
      U = identity_matrix(A)
      weak_popov!(P, W, U, true, true)
      return P, W, U
   else
      U = similar(A, 0, 0)
      weak_popov!(P, W, U, true, false)
      return P, W
   end
end

function find_pivot_popov(P::MatElem{T}, r::Int, last_col::Int = 0) where {T <: PolyElem}
   last_col == 0 ? n = ncols(P) : n = last_col
   pivot = n
   for c = n-1:-1:1
      if degree(P[r,c]) > degree(P[r,pivot])
         pivot = c
      end
   end
   return pivot
end

function init_pivots_popov(P::MatElem{T}, last_row::Int = 0, last_col::Int = 0) where {T <: PolyElem}
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

function weak_popov!(P::MatElem{T}, W::MatElem{T}, U::MatElem{T}, extended::Bool = false,
                                       with_trafo::Bool = false, last_row::Int = 0, last_col::Int = 0) where {T <: PolyElem}
   pivots = init_pivots_popov(P, last_row, last_col)
   weak_popov_with_pivots!(P, W, U, pivots, extended, with_trafo, last_row, last_col)
   return nothing
end

#=
The weak Popov form is defined by T. Mulders and A. Storjohann in
"On lattice reduction for polynomial matrices"
=#
function weak_popov_with_pivots!(P::MatElem{T}, W::MatElem{T}, U::MatElem{T}, pivots::Array{Array{Int,1}},
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
            q = -div(P[pivots[i][j], i], P[pivot, i])
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
    rank_profile_popov(A::MatElem{T}) where {T <: PolyElem}

Return an array of $r$ row indices such that these rows of $A$ are linearly
independent, where $r$ is the rank of $A$.
"""
function rank_profile_popov(A::MatElem{T}) where {T <: PolyElem}
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

function det_popov(A::MatElem{T}) where {T <: PolyElem}
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
         q = -div(B[r1, c], B[r2, c])
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
    popov(A::MatElem{T}) where {T <: PolyElem}

Return the Popov form of $A$.
"""
function popov(A::MatElem{T}) where {T <: PolyElem}
   return _popov(A, Val{false})
end

@doc Markdown.doc"""
    popov_with_transform(A::MatElem{T}) where {T <: PolyElem}

Compute a tuple $(P, U)$ where $P$ is the Popov form of $A$ and $U$
is a transformation matrix so that $P = UA$.
"""
function popov_with_transform(A::MatElem{T}) where {T <: PolyElem}
   return _popov(A, Val{true})
end

function _popov(A::MatElem{T}, trafo::Type{Val{S}} = Val{false}) where {T <: PolyElem, S}
   P = deepcopy(A)
   m = nrows(P)
   if trafo == Val{true}
      U = identity_matrix(A, m)
      popov!(P, U, true)
      return P, U
   else
      U = similar(A, 0, 0)
      popov!(P, U, false)
      return P
   end
end

function asc_order_popov!(P::MatElem{T}, U::MatElem{T}, pivots::Array{Array{Int,1}}, with_trafo::Bool) where {T <: PolyElem}
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
      j = findfirst(isequal(r), row_nums)
      row_nums[i] = r
      row_nums[j] = r2
   end
   return nothing
end

# Mulders, Storjohann: "On lattice reduction for polynomial matrices", Section 7
function popov!(P::MatElem{T}, U::MatElem{T}, with_trafo::Bool = false) where {T <: PolyElem}
   m = nrows(P)
   n = ncols(P)
   W = similar(U, 0, 0)
   pivots = init_pivots_popov(P)
   weak_popov_with_pivots!(P, W, U, pivots, false, with_trafo)
   asc_order_popov!(P, U, pivots, with_trafo)
   pivotColumns = zeros(Int, m)
   for c = 1:n
      if length(pivots[c]) == 0
         continue
      end
      # P is in weak Popov form, so any non-zero row contains exactly one pivot
      pivotColumns[pivots[c][1]] = c
   end

   t = base_ring(P)()
   for r = 1:m
      # Reduce row r assuming that rows 1, ..., r - 1 are reduced
      if iszero(pivotColumns[r])
         continue
      end
      r2 = 1
      while !iszero(r2)
         r2 = 0
         maxDegreeDiff = -1
         # Find the column c2 for which the difference
         # degree(P[r, c2]) - degree(P[r2, c2]) is maximal, where r2 is the row
         # in which the pivot is in column c2.
         for i = 1:r - 1
            ci = pivotColumns[i]
            if iszero(ci)
               continue
            end
            d = degree(P[r, ci]) - degree(P[i, ci])
            if d <= maxDegreeDiff
               continue
            end
            maxDegreeDiff = d
            r2 = i
         end
         iszero(r2) ? break : nothing
         c2 = pivotColumns[r2]
         q = -div(P[r, c2], P[r2, c2])
         for c = 1:n
            t = mul!(t, q, P[r2, c])
            P[r, c] = addeq!(P[r, c], t)
         end
         if with_trafo
            for c = 1:ncols(U)
               t = mul!(t, q, U[r2, c])
               U[r, c] = addeq!(U[r, c], t)
            end
         end
      end
   end
   # Make the pivots monic
   for i = 1:n
      if length(pivots[i]) == 0
         continue
      end
      r = pivots[i][1]
      cu = canonical_unit(P[r, i])
      if cu != 1
         for j = 1:n
            P[r, j] = divexact(P[r, j], cu)
         end
         if with_trafo
            for j = 1:ncols(U)
               U[r, j] = divexact(U[r, j], cu)
            end
         end
      end
   end
   return nothing
end

function hnf_via_popov(A::MatElem{T}) where {T <: PolyElem}
   return _hnf_via_popov(A, Val{false})
end

function hnf_via_popov_with_transform(A::MatElem{T}) where {T <: PolyElem}
   return _hnf_via_popov(A, Val{true})
end

function _hnf_via_popov(A::MatElem{T}, trafo::Type{Val{S}} = Val{false}) where {T <: PolyElem, S}
   H = deepcopy(A)
   m = nrows(H)
   if trafo == Val{true}
      U = identity_matrix(A, m)
      hnf_via_popov!(H, U, true)
      return H, U
   else
      U = similar(A, 0, 0)
      hnf_via_popov!(H, U, false)
      return H
   end
end

function hnf_via_popov_reduce_row!(H::MatElem{T}, U::MatElem{T}, pivots_hermite::Array{Int}, r::Int, with_trafo::Bool) where {T <: PolyElem}
   n = ncols(H)
   t = base_ring(H)()
   for c = 1:n
      if pivots_hermite[c] == 0
         continue
      end
      pivot = pivots_hermite[c]
      q = -div(H[r, c], H[pivot, c])
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

function hnf_via_popov_reduce_column!(H::MatElem{T}, U::MatElem{T}, pivots_hermite::Array{Int}, c::Int, with_trafo::Bool) where {T <: PolyElem}
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
      q = -div(H[i, c], H[r, c])
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

function hnf_via_popov!(H::MatElem{T}, U::MatElem{T}, with_trafo::Bool = false) where {T <: PolyElem}
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
         q = -div(H[r1, c], H[r2, c])
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
    similarity!(A::MatrixElem{T}, r::Int, d::T) where {T <: RingElement}

Applies a similarity transform to the $n\times n$ matrix $M$ in-place. Let
$P$ be the $n\times n$ identity matrix that has had all zero entries of row
$r$ replaced with $d$, then the transform applied is equivalent to
$M = P^{-1}MP$. We require $M$ to be a square matrix. A similarity transform
preserves the minimal and characteristic polynomials of a matrix.
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

Return a matrix $b$ with the entries of $a$, where the $i$th and $j$th
row are swapped.
"""
function swap_rows(a::MatrixElem, i::Int, j::Int)
   (1 <= i <= nrows(a) && 1 <= j <= nrows(a)) || throw(BoundsError())
   b = deepcopy(a)
   swap_rows!(b, i, j)
   return b
end

@doc Markdown.doc"""
    swap_rows!(a::MatrixElem, i::Int, j::Int)

Swap the $i$th and $j$th row of $a$.
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

Return a matrix $b$ with the entries of $a$, where the $i$th and $j$th
row are swapped.
"""
function swap_cols(a::MatrixElem, i::Int, j::Int)
   (1 <= i <= ncols(a) && 1 <= j <= ncols(a)) || throw(BoundsError())
   b = deepcopy(a)
   swap_cols!(b, i, j)
   return b
end

@doc Markdown.doc"""
    swap_cols!(a::MatrixElem, i::Int, j::Int)

Swap the $i$th and $j$th column of $a$.
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
    reverse_rows!(a::MatrixElem)

Swap the $i$th and $r - i$th row of $a$ for $1 \leq i \leq r/2$,
where $r$ is the number of rows of $a$.
"""
function reverse_rows!(a::MatrixElem)
   k = div(nrows(a), 2)
   for i in 1:k
      swap_rows!(a, i, nrows(a) - i + 1)
   end
   return a
end

@doc Markdown.doc"""
    reverse_rows(a::MatrixElem)

Return a matrix $b$ with the entries of $a$, where the $i$th and $r - i$th
row is swapped for $1 \leq i \leq r/2$. Here $r$ is the number of rows of
$a$.
"""
function reverse_rows(a::MatrixElem)
   b = deepcopy(a)
   return reverse_rows!(b)
end

@doc Markdown.doc"""
    reverse_cols!(a::MatrixElem)

Swap the $i$th and $r - i$th column of $a$ for $1 \leq i \leq c/2$,
where $c$ is the number of columns of $a$.
"""
function reverse_cols!(a::MatrixElem)
   k = div(ncols(a), 2)
   for i in 1:k
      swap_cols!(a, i, ncols(a) - i + 1)
   end
   return a
end

@doc Markdown.doc"""
    reverse_cols(a::MatrixElem)

Return a matrix $b$ with the entries of $a$, where the $i$th and $r - i$th
column is swapped for $1 \leq i \leq c/2$. Here $c$ is the number of columns
of$a$.
"""
function reverse_cols(a::MatrixElem)
   b = deepcopy(a)
   return reverse_cols!(b)
end

################################################################################
#
#  Elementary row/column transformations
#
################################################################################

@doc Markdown.doc"""
    add_column!(a::MatrixElem, s::RingElement, i::Int, j::Int, rows = 1:nrows(a))

Add $s$ times the $i$-th row to the $j$-th row of $a$.

By default, the transformation is applied to all rows of $a$. This can be
changed using the optional `rows` argument.
"""
function add_column!(a::MatrixElem, s::RingElement, i::Int, j::Int, rows = 1:nrows(a))
   v = base_ring(a)(s)
   nc = ncols(a)
   !_checkbounds(nc, i) && error("Column index ($i) must be between 1 and $nc")
   !_checkbounds(nc, j) && error("Column index ($j) must be between 1 and $nc")
   temp = base_ring(a)()
   for r in rows
      temp = mul!(temp, v, a[r, i])
      a[r, j] += temp # cannot mutate matrix entries
   end
   return a
end

@doc Markdown.doc"""
    add_column(a::MatrixElem, s::RingElement, i::Int, j::Int, rows = 1:nrows(a))

Create a copy of $a$ and add $s$ times the $i$-th row to the $j$-th row of $a$.

By default, the transformation is applied to all rows of $a$. This can be
changed using the optional `rows` argument.

"""
function add_column(a::MatrixElem, s::RingElement, i::Int, j::Int, rows = 1:nrows(a))
   b = deepcopy(a)
   return add_column!(b, s, i, j, rows)
end

@doc Markdown.doc"""
    add_row!(a::MatrixElem, s::RingElement, i::Int, j::Int, cols = 1:ncols(a))

Add $s$ times the $i$-th row to the $j$-th row of $a$.

By default, the transformation is applied to all columns of $a$. This can be
changed using the optional `cols` argument.
"""
function add_row!(a::MatrixElem, s::RingElement, i::Int, j::Int, cols = 1:ncols(a))
   v = base_ring(a)(s)
   nr = nrows(a)
   !_checkbounds(nr, i) && error("Row index ($i) must be between 1 and $nr")
   !_checkbounds(nr, j) && error("Row index ($j) must be between 1 and $nr")
   temp = base_ring(a)()
   for c in cols
      temp = mul!(temp, v, a[i, c])
      a[j, c] += temp # cannot mutate matrix entries
   end
   return a
end

@doc Markdown.doc"""
    add_row(a::MatrixElem, s::RingElement, i::Int, j::Int, cols = 1:ncols(a))

Create a copy of $a$ and add $s$ times the $i$-th row to the $j$-th row of $a$.

By default, the transformation is applied to all columns of $a$. This can be
changed using the optional `cols` argument.
"""
function add_row(a::MatrixElem, s::RingElement, i::Int, j::Int, cols = 1:ncols(a))
   b = deepcopy(a)
   return add_row!(b, s, i, j, cols)
end

# Multiply column

@doc Markdown.doc"""
    multiply_column!(a::MatrixElem, s::RingElement, i::Int, rows = 1:nrows(a))

Multiply the $i$th column of $a$ with $s$.

By default, the transformation is applied to all rows of $a$. This can be
changed using the optional `rows` argument.
"""
function multiply_column!(a::MatrixElem, s::RingElement, i::Int, rows = 1:nrows(a))
   c = base_ring(a)(s)
   nc = ncols(a)
   !_checkbounds(nc, i) && error("Row index ($i) must be between 1 and $nc")
   temp = base_ring(a)()
   for r in rows
      a[r, i] = c*a[r, i] # cannot mutate matrix entries
   end
   return a
end

@doc Markdown.doc"""
    multiply_column(a::MatrixElem, s::RingElement, i::Int, rows = 1:nrows(a))

Create a copy of $a$ and multiply the $i$th column of $a$ with $s$.

By default, the transformation is applied to all rows of $a$. This can be
changed using the optional `rows` argument.
"""
function multiply_column(a::MatrixElem, s::RingElement, i::Int, rows = 1:nrows(a))
   b = deepcopy(a)
   return multiply_column!(b, s, i, rows)
end

# Multiply row

@doc Markdown.doc"""
    multiply_row!(a::MatrixElem, s::RingElement, i::Int, cols = 1:ncols(a))

Multiply the $i$th row of $a$ with $s$.

By default, the transformation is applied to all columns of $a$. This can be
changed using the optional `cols` argument.
"""
function multiply_row!(a::MatrixElem, s::RingElement, i::Int, cols = 1:ncols(a))
   c = base_ring(a)(s)
   nr = nrows(a)
   !_checkbounds(nr, i) && error("Row index ($i) must be between 1 and $nr")
   temp = base_ring(a)()
   for r in cols
      a[i, r] = c*a[i, r] # cannot mutate matrix entries
   end
   return a
end

@doc Markdown.doc"""
    multiply_row(a::MatrixElem, s::RingElement, i::Int, cols = 1:ncols(a))

Create a copy of $a$ and multiply  the $i$th row of $a$ with $s$.

By default, the transformation is applied to all columns of $a$. This can be
changed using the optional `cols` argument.
"""
function multiply_row(a::MatrixElem, s::RingElement, i::Int, cols = 1:ncols(a))
   b = deepcopy(a)
   return multiply_row!(b, s, i, cols)
end

###############################################################################
#
#   Concatenation
#
###############################################################################

@doc Markdown.doc"""
    hcat(a::MatElem, b::MatElem)

Return the horizontal concatenation of $a$ and $b$. Assumes that the
number of rows is the same in $a$ and $b$.
"""
function hcat(a::MatElem, b::MatElem)
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
    vcat(a::MatElem, b::MatElem)

Return the vertical concatenation of $a$ and $b$. Assumes that the
number of columns is the same in $a$ and $b$.
"""
function vcat(a::MatElem, b::MatElem)
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

@doc Markdown.doc"""
    vcat(A::MatrixElem...) -> MatrixElem

Return the horizontal concatenation of the matrices $A$.
All component matrices need to have the same base ring and number of columns.
"""
function Base.vcat(A::MatrixElem...)
  return _vcat(A)
end

Base.reduce(::typeof(vcat), A::AbstractVector{<:MatrixElem}) = _vcat(A)

function _vcat(A)
  if length(A) == 0
    error("Number of matrices to concatenate must be positive")
  end

  if any(x -> ncols(x) != ncols(A[1]), A)
    error("Matrices must have the same number of columns")
  end

  if any(x -> base_ring(x) != base_ring(A[1]), A)
    error("Matrices must have the same base ring")
  end

  M = similar(A[1], sum(nrows, A), ncols(A[1]))
  s = 0
  for N in A
    for j in 1:nrows(N)
      for k in 1:ncols(N)
        M[s+j, k] = N[j,k]
      end
    end
    s += nrows(N)
  end
  return M
end

@doc Markdown.doc"""
    hcat(A::MatrixElem...) -> MatrixElem

Return the horizontal concatenating of the matrices $A$.
All component matrices need to have the same base ring and number of rows.
"""
function Base.hcat(A::MatrixElem...)
  return _hcat(A)
end

Base.reduce(::typeof(hcat), A::AbstractVector{<:MatrixElem}) = _hcat(A)

function _hcat(A)
  if length(A) == 0
    error("Number of matrices to concatenate must be positive")
  end

  if any(x -> nrows(x) != nrows(A[1]), A)
    error("Matrices must have the same number of rows")
  end

  if any(x -> base_ring(x) != base_ring(A[1]), A)
    error("Matrices must have the same base ring")
  end

  M = similar(A[1], nrows(A[1]), sum(ncols, A))
  s = 0
  for N in A
    for j in 1:ncols(N)
      for k in 1:nrows(N)
        M[k, s + j] = N[k, j]
      end
    end
    s += ncols(N)
  end
  return M
end

function Base.cat(A::MatrixElem...;dims)
  @assert dims == (1,2) || isa(dims, Int)

  if isa(dims, Int)
    if dims == 1
      return hcat(A...)
    elseif dims == 2
      return vcat(A...)
    else
      error("dims must be 1, 2, or (1,2)")
    end
  end

  local X
  for i in 1:length(A)
    if i == 1
      X = hcat(A[1], zero(A[1], nrows(A[1]), sum(Int[ncols(A[j]) for j=2:length(A)])))
    else
      X = vcat(X, hcat(zero(A[1], nrows(A[i]), sum(ncols(A[j]) for j=1:i-1)), A[i], zero(A[1], nrows(A[i]), sum(Int[ncols(A[j]) for j in (i+1):length(A)]))))
    end
  end
  return X
end

function Base.hvcat(rows::Tuple{Vararg{Int}}, A::MatrixElem...)
  nr = 0
  k = 1
  for i in 1:length(rows)
    nr += nrows(A[k])
    k += rows[i]
  end

  nc = sum(ncols(A[i]) for i in 1:rows[1])

  M = similar(A[1], nr, nc)
  mat_offset = 0
  row_offset = 0
  for j in 1:length(rows)
    s = 0
    for i in 1:rows[j]
      N = A[mat_offset + i]
      for l in 1:ncols(N)
        for k in 1:nrows(N)
          M[row_offset + k, s + l] = N[k, l]
        end
      end
      s += ncols(N)
    end
    row_offset += nrows(A[1+ mat_offset])
    mat_offset += rows[j]
  end

  return M
end

###############################################################################
#
#   Change Base Ring
#
###############################################################################

# like change_base_ring, but without initializing the entries
# this function exists until a better API is implemented
_change_base_ring(R::Ring, a::MatElem) = zero_matrix(R, nrows(a), ncols(a))
_change_base_ring(R::Ring, a::MatAlgElem) = MatrixAlgebra(R, nrows(a))()

@doc Markdown.doc"""
    change_base_ring(R::Ring, M::MatrixElem)

Return the matrix obtained by coercing each entry into `R`.
"""
function change_base_ring(R::Ring, M::MatrixElem)
   N = _change_base_ring(R, M)
   for i = 1:nrows(M), j = 1:ncols(M)
      N[i,j] = R(M[i,j])
   end
   return N
end

###############################################################################
#
#   Map
#
###############################################################################

@doc Markdown.doc"""
    map_entries!(f, dst::MatrixElem, src::MatrixElem)

Like `map_entries`, but stores the result in `dst` rather than a new matrix.
"""
function map_entries!(f, dst::MatrixElem, src::MatrixElem)
   for i = 1:nrows(src), j = 1:ncols(src)
      dst[i, j] = f(src[i, j])
   end
   dst
end

@doc Markdown.doc"""
    map!(f, dst::MatrixElem, src::MatrixElem)

Like `map`, but stores the result in `dst` rather than a new matrix.
This is equivalent to `map_entries!(f, dst, src)`.
"""
Base.map!(f, dst::MatrixElem, src::MatrixElem) = map_entries!(f, dst, src)

@doc Markdown.doc"""
    map_entries(f, a::MatrixElem)

Transform matrix `a` by applying `f` on each element.
"""
function map_entries(f, a::MatrixElem)
   isempty(a) && return _change_base_ring(parent(f(zero(base_ring(a)))), a)
   b11 = f(a[1, 1])
   b = _change_base_ring(parent(b11), a)
   b[1, 1] = b11
   for i = 1:nrows(a), j = 1:ncols(a)
      i == j == 1 && continue
      b[i, j] = f(a[i, j])
   end
   b
end

@doc Markdown.doc"""
    map(f, a::MatrixElem)

Transform matrix `a` by applying `f` on each element.
This is equivalent to `map_entries(f, a)`.
"""
Base.map(f, a::MatrixElem) = map_entries(f, a)

###############################################################################
#
#   Random generation
#
###############################################################################

RandomExtensions.maketype(S::MatSpace, _) = elem_type(S)

function RandomExtensions.make(S::MatSpace, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1]) # forward to default Make constructor
   else
      make(S, make(R, vs...))
   end
end

# Sampler for a MatSpace not needing arguments (passed via make)
# this allows to obtain the Sampler in simple cases without having to know about make
# (when one can do `rand(M)`, one can expect to be able to do `rand(Sampler(rng, M))`)
Random.Sampler(::Type{RNG}, S::MatSpace, n::Random.Repetition
               ) where {RNG<:AbstractRNG} =
   Random.Sampler(RNG, make(S), n)

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{<:MatElem,
                                         <:MatSpace}})
   S, v = sp[][1:end]
   M = S()
   R = base_ring(S)
   for i = 1:nrows(M)
      for j = 1:ncols(M)
         M[i, j] = rand(rng, v)
      end
   end
   return M
end

rand(rng::AbstractRNG, S::MatSpace, v...) = rand(rng, make(S, v...))

rand(S::MatSpace, v...) = rand(Random.GLOBAL_RNG, S, v...)

# resolve ambiguities
rand(rng::AbstractRNG, S::MatSpace, dims::Integer...) =
   rand(rng, make(S), dims...)

rand(S::MatSpace, dims::Integer...) = rand(Random.GLOBAL_RNG, S, dims...)

function randmat_triu(rng::AbstractRNG, S::MatSpace, v...)
   M = S()
   R = base_ring(S)
   for i = 1:nrows(M)
      for j = 1:i - 1
         M[i, j] = R()
      end
      for j = i:ncols(M)
         M[i, j] = rand(rng, R, v...)
      end
      while iszero(M[i, i])
         M[i, i] = rand(rng, R, v...)
      end
   end
   return M
end

randmat_triu(S::MatSpace, v...) = randmat_triu(Random.GLOBAL_RNG, S, v...)

function randmat_with_rank(rng::AbstractRNG, S::MatSpace{T}, rank::Int, v...) where {T <: RingElement}
   if !isdomain_type(T) && !(T <: ResElem)
      error("Not implemented")
   end
   M = S()
   R = base_ring(S)
   for i = 1:rank
      for j = 1:i - 1
         M[i, j] = R()
      end
      M[i, i] = rand(rng, R, v...)
      while iszero(M[i, i])
         M[i, i] = rand(rng, R, v...)
      end
      for j = i + 1:ncols(M)
         M[i, j] = rand(rng, R, v...)
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
         r1 = rand(rng, 1:m)
         r2 = rand(rng, 1:m - 1)
         r2 = r2 >= r1 ? r2 + 1 : r2
         d = rand(rng, -5:5)
         for j = 1:ncols(M)
            M[r1, j] = M[r1, j] + d*M[r2, j]
         end
      end
   end
   return M
end

randmat_with_rank(S::MatSpace{T}, rank::Int, v...) where {T <: RingElement} =
   randmat_with_rank(Random.GLOBAL_RNG, S, rank, v...)

################################################################################
#
#   Matrix constructors
#
################################################################################

@doc Markdown.doc"""
    matrix(R::Ring, arr::AbstractArray{T, 2}) where {T}

Constructs the matrix over $R$ with entries as in `arr`.
"""
function matrix(R::Ring, arr::AbstractArray{T, 2}) where {T}
   if elem_type(R) === T
      z = Generic.MatSpaceElem{elem_type(R)}(arr)
      z.base_ring = R
      return z
   else
      arr_coerce = convert(Array{elem_type(R), 2}, map(R, arr))::Array{elem_type(R), 2}
      return matrix(R, arr_coerce)
   end
end

@doc Markdown.doc"""
    matrix(R::Ring, r::Int, c::Int, arr::AbstractArray{T, 1}) where {T}

Constructs the $r \times c$ matrix over $R$, where the entries are taken
row-wise from `arr`.
"""
function matrix(R::Ring, r::Int, c::Int, arr::AbstractVecOrMat{T}) where T
   _check_dim(r, c, arr)
   ndims(arr) == 2 && return matrix(R, arr)
   if elem_type(R) === T
     z = Generic.MatSpaceElem{elem_type(R)}(r, c, arr)
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
    zero_matrix(R::Ring, r::Int, c::Int)

Return the $r \times c$ zero matrix over $R$.
"""
function zero_matrix(R::Ring, r::Int, c::Int)
   arr = Array{elem_type(R)}(undef, r, c)
   for i in 1:r
      for j in 1:c
         arr[i, j] = zero(R)
      end
   end
   z = Generic.MatSpaceElem{elem_type(R)}(arr)
   z.base_ring = R
   return z
end

################################################################################
#
#   Identity matrix
#
################################################################################

@doc Markdown.doc"""
    identity_matrix(R::Ring, n::Int)

Return the $n \times n$ identity matrix over $R$.
"""
identity_matrix(R::Ring, n::Int) = diagonal_matrix(one(R), n)

################################################################################
#
#   Identity matrix
#
################################################################################

@doc Markdown.doc"""
    identity_matrix(M::MatElem{T}) where T <: RingElement

Construct the identity matrix in the same matrix space as `M`, i.e.
with ones down the diagonal and zeroes elsewhere. `M` must be square.
This is an alias for `one(M)`.
"""
function identity_matrix(M::MatElem{T}) where T <: RingElement
   identity_matrix(check_square(M), nrows(M))
end

function identity_matrix(M::MatElem{T}, n::Int) where T <: RingElement
   z = zero(M, n, n)
   R = base_ring(M)
   for i = 1:n
      z[i, i] = one(R)
   end
   z
end

################################################################################
#
#   Diagonal matrix
#
################################################################################

@doc Markdown.doc"""
    diagonal_matrix(x::RingElement, m::Int, [n::Int])

Return the $m \times n$ matrix over $R$ with `x` along the main diagonal and
zeroes elsewhere. If `n` is not specified, it defaults to `m`.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> diagonal_matrix(ZZ(2), 2, 3)
[2   0   0]
[0   2   0]

julia> diagonal_matrix(QQ(-1), 3)
[-1//1    0//1    0//1]
[ 0//1   -1//1    0//1]
[ 0//1    0//1   -1//1]
```
"""
function diagonal_matrix(x::RingElement, m::Int, n::Int)
   z = zero_matrix(parent(x), m, n)
   for i in 1:min(m, n)
      z[i, i] = x
   end
   return z
end

diagonal_matrix(x::RingElement, m::Int) = diagonal_matrix(x, m, m)

###############################################################################
#
#   MatrixSpace constructor
#
###############################################################################

@doc Markdown.doc"""
    MatrixSpace(R::Ring, r::Int, c::Int; cached::Bool = true)

Return parent object corresponding to the space of $r\times c$ matrices over
the ring $R$. If `cached == true` (the default), the returned parent object
is cached so that it can returned by future calls to the constructor with the
same dimensions and base ring.
"""
function MatrixSpace(R::Ring, r::Int, c::Int; cached::Bool = true)
   return Generic.MatrixSpace(R, r, c, cached=cached)
end

###############################################################################
#
#   Conversion to Array
#
###############################################################################

"""
    Matrix(A::MatrixElem)

Convert `A` to a Julia `Matrix` of the same dimensions with the same elements.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> A = ZZ[1 2 3; 4 5 6]
[1   2   3]
[4   5   6]

julia> Matrix(A)
2×3 Array{BigInt,2}:
 1  2  3
 4  5  6
```
"""
Matrix(M::MatrixElem) = eltype(M)[M[i, j] for i = 1:nrows(M), j = 1:ncols(M)]

"""
    Array(A::MatrixElem)

Convert `A` to a Julia `Matrix` of the same dimensions with the same elements.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> R, x = ZZ["x"]; A = R[x^0 x^1; x^2 x^3]
[  1     x]
[x^2   x^3]

julia> Array(A)
2×2 Array{AbstractAlgebra.Generic.Poly{BigInt},2}:
 1    x
 x^2  x^3
```
"""
Array(M::MatrixElem) = Matrix(M)
