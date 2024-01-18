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
end

@doc raw"""
    echelon_form_with_transformation(A::MatElem{<:FieldElement}; reduced::Bool = true, shape::Symbol = :upper)

Return a row echelon form $R$ of $A$ and an invertible matrix $U$ with $R = UA$.

See [`echelon_form`](@ref) for the keyword arguments.
"""
function echelon_form_with_transformation(A::MatElem{<:FieldElement}; reduced::Bool = true, shape::Symbol = :upper)
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
end

@doc raw"""
    hermite_form_with_transformation(A::MatElem{<:RingElement}; reduced::Bool = true, shape::Symbol = :upper)

Return a Hermite normal form $H$ of $A$ and an invertible matrix $U$ with $H = UA$.
It is assumed that `base_ring(A)` is euclidean.

See [`hermite_form`](@ref) for the keyword arguments.
"""
function hermite_form_with_transformation(A::MatElem{<:RingElement}; reduced::Bool = true, shape::Symbol = :upper)
end
