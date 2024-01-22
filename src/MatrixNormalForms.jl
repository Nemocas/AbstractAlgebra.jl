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
#  Matrix normal form over (euclidean) rings (hermite_form)
#
################################################################################

@doc raw"""
    hermite_form(A::MatElem{<:RingElement}; reduced::Bool = true, shape::Symbol = :upper, trim::Bool = false)

Return a Hermite normal form $H$ of $A$.
It is assumed that `base_ring(A)` is euclidean.

# Keyword arguments
* `reduced`: Whether the columns of $H$ which contain a pivot are reduced. The
  default is `true`.
* `shape`: Whether $R$ is an upper-right (`:upper`, default) or lower-left
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
  H = hnf_kb(A)
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
  H, U = hnf_kb_with_transform(A)
  if shape === :lower
    reverse_cols!(H)
    reverse_rows!(H)
    reverse_rows!(U)
  end
  return H, U
end
