################################################################################
#
#  Matrix normal form over fields (echelon_form)
#
################################################################################

@doc raw"""
    echelon_form(A::MatElem{<:FieldElement}; reduced::Bool = true, shape::Symbol = :upper, trim::Bool = false)

Return a row echelon form $R$ of $A$.

# Keyword arguments
* `reduced`: Whether the columns of $R$ which contain a pivot are reduced. The
  default is `true`.
* `shape`: Whether $R$ is an upper-right (`:upper`, default) or lower-left
  (`:lower`) echelon form.
* `trim`: By default, $R$ will have as many rows as $A$ and potentially involve
  rows only containing zeros. If `trim = true`, these rows are removed, so that
  the number of rows of $R$ coincides with the rank of $A$.

See also [`echelon_form_with_transformation`](@ref).
"""
function echelon_form(A::MatElem{<:FieldElement}; reduced::Bool = true, shape::Symbol = :upper, trim::Bool = false)
  R = deepcopy(A)
  r = echelon_form!(R, reduced = reduced, shape = shape)
  if trim
    if shape === :upper
      return sub(R, 1:r, 1:ncols(R))
    else
      return sub(R, nrows(R) - r + 1:nrows(R), 1:ncols(R))
    end
  end
  return R
end

@doc raw"""
    echelon_form_with_transformation(A::MatElem{<:FieldElement}; reduced::Bool = true, shape::Symbol = :upper)

Return a row echelon form $R$ of $A$ and an invertible matrix $U$ with $R = UA$.

See [`echelon_form`](@ref) for the keyword arguments.
"""
function echelon_form_with_transformation(A::MatElem{<:FieldElement}; reduced::Bool = true, shape::Symbol = :upper)
  if shape === :upper
    R = hcat(A, identity_matrix(base_ring(A), nrows(A)))
  else
    R = hcat(identity_matrix(base_ring(A), nrows(A)), A)
  end
  echelon_form!(R, reduced = reduced, shape = shape)
  if shape === :upper
    return sub(R, 1:nrows(A), 1:ncols(A)), sub(R, 1:nrows(A), ncols(A) + 1:ncols(R))
  else
    return sub(R, 1:nrows(A), nrows(A) + 1:ncols(R)), sub(R, 1:nrows(A), 1:nrows(A))
  end
end

# Return the rank of A and put A into echelon form
# So far, the `reduced` keyword is ignored (that is, the result is always reduced)
function echelon_form!(A::MatElem{<:FieldElement}; reduced::Bool = true, shape::Symbol = :upper)
  if shape !== :upper && shape !== :lower
    throw(ArgumentError("Unsupported argument :$shape for shape: Must be :upper or :lower."))
  end

  if shape === :lower
    reverse_cols!(A)
    r = echelon_form!(A, reduced = reduced, shape = :upper)
    reverse_cols!(A)
    reverse_rows!(A)
    return r
  end

  return rref!(A)
end

################################################################################
#
#  Matrix normal form over (euclidean) domains (hermite_form)
#
################################################################################

@doc raw"""
    hermite_form(A::MatElem{<:RingElement}; reduced::Bool = true, shape::Symbol = :upper, trim::Bool = false)

Return a Hermite normal form $H$ of $A$.
It is assumed that `base_ring(A)` is euclidean.

# Keyword arguments
* `reduced`: Whether the columns of $H$ which contain a pivot are reduced. The
  default is `true`.
* `shape`: Whether $H$ is an upper-right (`:upper`, default) or lower-left
  (`:lower`) echelon form.
* `trim`: By default, $H$ will have as many rows as $A$ and potentially involve
  rows only containing zeros. If `trim = true`, these rows are removed, so that
  the number of rows of $H$ coincides with the rank of $A$.

See also [`hermite_form_with_transformation`](@ref).
"""
function hermite_form(A::MatElem{<:RingElement}; reduced::Bool = true, shape::Symbol = :upper, trim::Bool = false)
  if shape !== :upper && shape !== :lower
    throw(ArgumentError("Unsupported argument :$shape for shape: Must be :upper or :lower."))
  end

  if shape === :lower
    A = reverse_cols(A)
  end
  H = hnf(A)
  r = nrows(H)
  # Compute the rank (if necessary)
  if trim
    while r > 0 && is_zero_row(H, r)
      r -= 1
    end
    H = sub(H, 1:r, 1:ncols(H))
  end
  if shape === :lower
    reverse_cols!(H)
    reverse_rows!(H)
  end
  return H
end

@doc raw"""
    hermite_form_with_transformation(A::MatElem{<:RingElement}; reduced::Bool = true, shape::Symbol = :upper)

Return a Hermite normal form $H$ of $A$ and an invertible matrix $U$ with $H = UA$.
It is assumed that `base_ring(A)` is euclidean.

See [`hermite_form`](@ref) for the keyword arguments.
"""
function hermite_form_with_transformation(A::MatElem{<:RingElement}; reduced::Bool = true, shape::Symbol = :upper)
  if shape !== :upper && shape !== :lower
    throw(ArgumentError("Unsupported argument :$shape for shape: Must be :upper or :lower."))
  end

  if shape === :lower
    A = reverse_cols(A)
  end
  H, U = hnf_with_transform(A)
  if shape === :lower
    reverse_cols!(H)
    reverse_rows!(H)
    reverse_rows!(U)
  end
  return H, U
end

################################################################################
#
#  Matrix normal form over principal ideal rings (howell_form)
#
################################################################################

# Works in theory over any principal ideal ring; internally we require functions
# annihilator, gcdxx and _div_for_howell_form

# Swap rows so that there is a non-zero entry in A[start_row, col].
# Return 0 if this is not possible, 1 if no swapping was necessary and -1
# if rows were swapped.
function _pivot(A::MatElem, start_row::Int, col::Int)
  if !is_zero_entry(A, start_row, col)
    return 1
  end

  for j in start_row + 1:nrows(A)
    if !is_zero_entry(A, j, col)
      swap_rows!(A, j, start_row)
      return -1
    end
  end

  return 0
end

function triangularize!(A::MatElem{<:RingElement})
  n = nrows(A)
  m = ncols(A)
  d = one(base_ring(A))
  row = 1
  col = 1
  while row <= nrows(A) && col <= ncols(A)
    t = _pivot(A, row, col)
    if iszero(t)
      col = col + 1
      continue
    end
    d = d*t
    for i in (row + 1):nrows(A)
      if is_zero_entry(A, i, col)
        continue
      end

      b, q = divides(A[i, col], A[row, col])

      if b
        for k in col:m
          A[i, k] = A[i, k] - q*A[row, k]
        end
      else
        g, s, t, u, v = gcdxx(A[row, col], A[i, col])

        for k in col:m
          t1 = s*A[row, k] + t*A[i, k]
          t2 = u*A[row, k] + v*A[i, k]
          A[row, k] = t1
          A[i, k] = t2
        end
      end
    end
    row = row + 1
    col = col + 1
  end
  return d
end

# Naive version of inplace strong echelon form
# It is assumed that A has more rows then columns.
function strong_echelon_form_naive!(A::MatElem{<:RingElement})
  n = nrows(A)
  m = ncols(A)

  @assert n >= m

  triangularize!(A)

  T = zero_matrix(base_ring(A), 1, ncols(A))
  for j in 1:m
    if !is_zero_entry(A, j, j)
      # Normalize/canonicalize the pivot
      u = canonical_unit(A[j, j])
      if !is_one(u)
        uinv = inv(u)
        for i in j:ncols(A)
          A[j, i] = uinv*A[j, i]
        end
      end

      # This is the reduction
      for i in 1:j - 1
        if is_zero_entry(A, i, j)
          continue
        end
        q = _div_for_howell_form(A[i, j], A[j, j])
        for l in i:m
          A[i, l] = A[i, l] - q*A[j, l]
        end
      end

      a = annihilator(A[j, j])

      for k in 1:m
        T[1, k] = a*A[j, k]
      end
    else
      for k in 1:m
        T[1, k] = A[j, k]
      end
    end

    for i in j+1:m

      if is_zero_entry(T, 1, i)
        continue
      end

      if is_zero_entry(A, i, i)
        for k in i:m
          T[1, k], A[i, k] = A[i, k], T[1, k]
        end
      else
        b, q = divides(T[1, i], A[i, i])
        if b
          for k in i:m
            T[1, k] = T[1, k] - q*A[i, k]
          end
        else
          g, s, t, u, v = gcdxx(A[i, i], T[1, i])
          for k in i:m
            t1 = s*A[i, k] + t*T[1, k]
            t2 = u*A[i, k] + v*T[1, k]
            A[i, k] = t1
            T[1, k] = t2
          end
        end
      end
    end
  end
  return A
end

function howell_form!(A::MatElem{<:RingElement})
  @assert nrows(A) >= ncols(A)

  strong_echelon_form_naive!(A)

  # Move the zero rows to the bottom
  for i in 1:nrows(A)
    if !is_zero_row(A, i)
      continue
    end

    j = findfirst(l -> !is_zero_row(A, l), i + 1:nrows(A))
    if isnothing(j)
      break
    end
    swap_rows!(A, i, i + j)
  end
  return A
end

@doc raw"""
    howell_form(A::MatElem{<:RingElement}; reduced::Bool = true, shape::Symbol = :upper, trim::Bool = false)

Return a Howell form $H$ of $A$.
It is assumed that `base_ring(A)` is a principal ideal ring.

# Keyword arguments
* `reduced`: Whether the columns of $H$ which contain a pivot are reduced. The
  default is `true`.
* `shape`: Whether $H$ is an upper-right (`:upper`, default) or lower-left
  (`:lower`) echelon form.
* `trim`: By default, $H$ will have at least as many rows as $A$ and potentially
  involve rows only containing zeros. If `trim = true`, these rows are removed.

See also [`howell_form_with_transformation`](@ref).
"""
function howell_form(A::MatElem{<:RingElement}; reduced::Bool = true, shape::Symbol = :upper, trim::Bool = false)
  if shape !== :upper && shape !== :lower
    throw(ArgumentError("Unsupported argument :$shape for shape: Must be :upper or :lower."))
  end

  H = deepcopy(A)
  if shape === :lower
    H = reverse_cols!(H)
  end

  if nrows(H) < ncols(H)
    H = vcat(H, zero_matrix(base_ring(H), ncols(H) - nrows(H), ncols(H)))
  end

  howell_form!(H)

  r = nrows(H)
  # Compute the rank (if necessary)
  if trim
    while r > 0 && is_zero_row(H, r)
      r -= 1
    end
    H = sub(H, 1:r, 1:ncols(H))
  end

  if shape === :lower
    reverse_cols!(H)
    reverse_rows!(H)
  end
  return H
end

@doc raw"""
    howell_form_with_transformation(A::MatElem{<:RingElement}; reduced::Bool = true, shape::Symbol = :upper)

Return a Howell form $H$ of $A$ and a matrix $U$ with $H = UA$.
Notice that $H$ may have more rows than $A$ and hence $U$ may not be invertible.
It is assumed that `base_ring(A)` is a principal ideal ring

See [`hermite_form`](@ref) for the keyword arguments.
"""
function howell_form_with_transformation(A::MatElem{<:RingElement}; reduced::Bool = true, shape::Symbol = :upper)
  if shape !== :upper && shape !== :lower
    throw(ArgumentError("Unsupported argument :$shape for shape: Must be :upper or :lower."))
  end

  if shape === :lower
    B = hcat(reverse_cols(A), identity_matrix(A, nrows(A)))
  else
    B = hcat(A, identity_matrix(A, nrows(A)))
  end
  if nrows(B) < ncols(B)
    B = vcat(B, zero(A, ncols(B) - nrows(B), ncols(B)))
  end

  howell_form!(B)

  m = max(nrows(A), ncols(A))
  H = sub(B, 1:m, 1:ncols(A))
  U = sub(B, 1:m, ncols(A) + 1:ncols(B))

  if shape === :lower
    reverse_cols!(H)
    reverse_rows!(H)
    reverse_rows!(U)
  end
  return H, U
end
