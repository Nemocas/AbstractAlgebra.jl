

# These functions are honestly just triangular solving in disguise.

function solve_fflu_precomp(p::Generic.Perm, FFLU::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   x = p * b
   n = nrows(x)
   m = ncols(x)
   R = base_ring(FFLU)

   t = base_ring(b)()
   s = base_ring(b)()
   minus_one = R(-1)

    # For each column of x
    for k in 1:m
        # Backsolve the lower triangular part. However, things have been rearranged
        # so that a single column is accessed in the inner (j) loop, rather than a row.
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

        # Backsolve the upper triangular part. However, things have been rearranged
        # so that a single column of `x` is accessed in the inner (j) loop, rather than a row.
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

function solve_lu_precomp(p::Generic.Perm, LU::MatElem{T}, b::MatrixElem{T}) where {T <: RingElement}

    n = nrows(LU)
    m = ncols(LU)
    
    rk = begin
        rk=min(n,m)
        for i = min(n,m):-1:1
            iszero(LU[i,i]) ? rk -= 1 : break
        end
        rk
    end

    ncolsb = ncols(b)
    R = base_ring(LU)

    #x = p * b # Not correct dimensions for output.

    # TODO: In general, `zero_matrix` will return a type dependent only on the base_ring,
    # which could in theory be a different concrete type from the input. It will be crucial
    # to implement a zero-matrix constructor that takes into account the type of `LU`.
    x  = zero_matrix(R, m, ncolsb)
    pb = p*b
    
    t = base_ring(b)()
    s = base_ring(b)()
    
    # For each column of b, solve Ax=b[:,j].
    for k in 1:ncolsb

        # The back-substitution, along rows.
        # Note that the lower-triangular part by definition has ones on the
        # diagonal. Thus, a solve is always possible. Moreover, for "tall"
        # solve instances (`n>m`), only the first `m` columns of `L` need to
        # be considered. This follows from the following lemma
        #=
            Lemma: Let D be a (not necessarily commutative) division algebra, and
                  `L,U` be matrices in `GL_n(D)`. Then $A = LU$ is the unique `LU`
                  decomposition for A.

            Proof: Write A = LU = (L')(U'). Now inv(L')*L == (U')*inv(U). Since
                   upper (resp. lower) triangular matrices form a subgroup of M_n(D),
                   which once can check by basic arithmetic, we see that L'==L and U'==U. []

            When `U` is singular, uniqueness fails to hold, but hopefully we can show
            up to the natural ambiguity that we still only need the first `m` columns of `L`.
            If you're reading this, I guess you've found a bug!
        =#

        # Populate the initial values in `x`. There is surely a better way to do this.
        for i=1:rk
            x[i, k] = deepcopy(pb[i, k])
        end
        
        for i in 2:rk
            for j in 1:(i - 1)
                # NOTE: This is the correct order of multiplication in non-commutative rings.
                # x[i, k] = x[i, k] - LU[i, j] * x[j, k]
                t = mul_red!(t, LU[i, j], x[j, k], false)
                t = minus!(t)
                if j == 1
                    # This was to allocate memory before. Now we don't have to.
                    # In fact, it is better to allocate memory in a single request
                    # rather than in small pieces.
                    
                    x[i, k] = x[i, k] + t # LU[i, j] * x[j, k]
                else
                    # This doesn't do what you think it does... There is an implicit
                    # copy in `setindex!`. At least one allocation is avoided, but not both.
                    # That said, the assignment is absolutely necessary for immutable input.
                    x[i, k] = addeq!(x[i, k], t)
                end
            end
            x[i, k] = reduce!(x[i, k])
        end

        # Below the last row where `U` has a pivot, the entries
        # of `LU` are just the entries of `L`.
        check_system_is_consistent(LU, view(x, : , k:k), view(b, :, k:k) , rk)
        
        # Since _ut_pivot_columns only checks entries above the diagonal,
        # it effectively only reads the `U` part of the `LU` object.
        pcols  = _ut_pivot_columns(LU)
        
        # PROOF OF CORRECT USAGE:
        # By use of `deepcopy`, we see `x` is freshly allocated space and the entries are
        # deepcopyed/newly allocated.
        # Thus, the elements of `x` share no references, even in part, aside from parents.
        # Moreover, trivially, `x===x`. Thus, we have fulfilled the CONTRACT for using
        # the `!!` method.        
        x[:,k] = _solve_nonsingular_ut!!_I_agree_to_the_terms_and_conditions_of_this_function(x[:,k], LU[:, pcols], x[:,k], rk, 1)

        #TODO: Replace the implicit `get_index` with a `view`.

    end
    return x
end

@doc Markdown.doc"""
    check_system_is_consistent(A,x,b,rk)

Checks if `A[rk+1:n, :]*x == b`, where `n` is the number of rows of `A`.
An error is raised otherwise.
"""
function check_system_is_consistent(A, x, b, rk = 0::Int)

    n = nrows(A)
    m = ncols(A)
    nvecs = max(ncols(x), ncols(b))

    # Containers
    t = base_ring(b)()
    sum = base_ring(b)()

    # Check all the dot products.
    for k = 1:nvecs
        for i = rk+1:n
            sum = zero!(sum)
            for j=1:m 
                sum = addmul_delayed_reduction!(sum, A[i, j], x[j, k], t)
            end
            sum = reduce!(sum)
            sum != b[i,k] && throw(DomainError((A, x, b), "Solve instance is inconsistent."))
        end
    end
    return true
end


# NOTE: Very dangerous function, and breaks the mutability edicts set out in
# https://github.com/Nemocas/Nemo.jl/issues/278
#
# CONTRACT:
#   (*) The contract does not cover cases where allocation and assignment functions (ex: deepcopy)
# have been maliciously or negligently tampered with, in which case the TYPE owner assumes
# all fault.
#
# By involking this function, you (the calling method) have intended target of mutation
# `x` and can certify that either:
#
# 1. `A,b` shares no references (even in part, aside from parents)
#     with `x`, or
#
# 2. `x===b`, the mutation on `b` (i.e `x`) is intended, and elements of `x` share no
#    references between each other, even in part, aside from parents.
#
# 3. The matrix `A` is in reduced upper-triangular form.
#
# 4. A correct proof that conditions 1-3 have been satisfied, given (*), exists at the call site.
#
#
# Under these conditions, I (the function) agree that
#
# 1. `A` and constituants are not modified in any way, shape, or form,
# 2. either `b===x`, or `b` and its constituants are not modified in any, shape, or form.
# 3. The output `x` will be a solution to `Ax==b`.
#
# Failure to abide by the terms of this contract absolves this function of any liability,
# and usually will result in misery and suffering.
#
# Your SIGNATURE of the contract is the `git-blame` info at the call site. MY signature is
# the `git-blame` info at this contract.
#
# The reason this function exists is that solve methods have the same pattern of
# back-substitution, but generally only want to allocate the memory for the solution
# vector once. `x` is assumed to be such an output container, which by definition of
# the mutability edicts should not share references from `A` or `b` in the caller.
#
function _solve_nonsingular_ut!!_I_agree_to_the_terms_and_conditions_of_this_function(x, U, b, rk, k)

    if iszero(rk)
        return x
    end

    n  = nrows(U)
    m  = ncols(U)
    #rk = min(n,m)::Int

    # Arithmetic container.
    t = base_ring(b)()
    ZERO = base_ring(b)()
    
    # NOTE: This is the correct order of division for non-commutative rings.
    x[rk, k] = divexact_left(U[rk, rk], b[rk, k])

    # The upper triangular back-substitution. Along rows.
    for i in rk-1:-1:1
        
        # Theoretically `x[i,k] = add!(x[i,k], b[i,k], ZERO)` should be sufficient,
        # but I really don't trust that `add!` has uniform behaviour.
        t = zero!(t)
        t = addeq!(t, b[i,k])
        x[i,k] = add!(x[i,k], t, ZERO)
        
        for j in (i + 1):rk
            # x[i, k] = x[i, k] - x[j, k] * U[i, j]
            t = mul_red!(t, x[j, k], -U[i, j], false)
            x[i, k] = addeq!(x[i, k], t)
        end
        x[i, k] = reduce!(x[i, k])
        x[i, k] = divexact_left(U[i, i], x[i, k])
    end

    return x
end

#=
@doc Markdown.doc"""
    _ut_pivot_columns(A::MatElem{<:RingElem})

Return the pivot columns of an upper triangular matrix, which is assumed to be
reduced (see docstring for pivot_columns).

Absolutely no check is made, which is surprisingly useful applied to factorizations
of matrices which are stored in-place.
"""
=#
function _ut_pivot_columns(A::MatElem{T}) where T
  p = Int[]
  j = 0
  for i=1:nrows(A)
    j += 1
    if j > ncols(A)
      return p
    end
    while iszero(A[i,j])
      j += 1
      if j > ncols(A)
        return p
      end
    end
    push!(p, j)
  end
  return p
end

# NOTE: This is a UI function. You actually want to call _pivot_columns in
# library code to skip the check.
@doc Markdown.doc"""
    pivot_columns(A::MatElem{<:RingElem})

Return the pivot columns of a matrix in reduced triangular form. A matrix is
in reduced triangular form if the multiset `{ min{i : A[i,j] != 0} }` has no 
elements of multiplicity 2 or more. (where the minimum over the empty set is empty).
i.e, the matrix looks something like this:

[* * * * * *]
[0 * * * * *]
[0 0 0 * * *]
[0 0 0 0 0 *]
[0 0 0 0 0 0]
"""
function pivot_columns(A::MatElem{T}) where T
    !isreduced_ut(A) && throw(DomainError(A, "Matrix is not in Hermite normal form."))
    return _pivot_columns(A)
end

function isreduced_ut(A::MatElem{T}) where T
    first_nonzeros = Int[findfirst(!iszero, A[j,:]) for j=1:nrows(A)]
    for i=1:length(first_nonzeros-1)
        first_nonzeros[i+1] <= first_nonzeros[i] && return false
    end
    return true
end

function solve_consistency_check()
    error("Not Implemented.")
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

#############################################################################
#
#  Kept for reference.
#
#############################################################################



### Code below is not active ###
if false
    
@doc Markdown.doc"""
    solve_ut(A::MatElem{T}, b::MatElem{T}) -> MatElem{T})

Given an upper triangular $m \times m$ matrix $A$ and a matrix $b$ of size $m
\times 1$, this function computes $x$ such that $Ax = b$.  It is assumed that
the pivots of $A$ are invertible.
"""
function solve_ut(A::MatElem{T}, b::MatElem{T}) where T
  m = nrows(A)
  n = ncols(A)
  @assert m == nrows(b)
  @assert m <= n
  x = zero_matrix(base_ring(A), n, 1)
  pivot_cols = Vector{Int}()
  r = 0
  last_pivot = n + 1
  for i = m:-1:1
    for j = 1:last_pivot - 1
      if iszero(A[i, j])
        continue
      end
      x[j, 1] = b[i, 1]
      for k = 1:r
        x[j, 1] -= A[i, pivot_cols[k]]*x[pivot_cols[k], 1]
      end
      x[j, 1] *= inv(A[i, j])
      last_pivot = j
      r += 1
      push!(pivot_cols, j)
      break
    end
  end
  return x
end

@doc Markdown.doc"""
    solve_lt(A::MatElem{T}, b::MatElem{T}) -> MatElem{T})

Given a lower triangular $m \times m$ matrix $A$ and a matrix $b$ of size
$m \times 1$, this function computes $x$ such that $Ax = b$.  It is assumed
that the pivots of $A$ are invertible.
"""
function solve_lt(A::MatElem{T}, b::MatElem{T}) where T
  m = nrows(A)
  n = ncols(A)
  @assert m == nrows(b)
  @assert m <= n
  x = zero_matrix(base_ring(A), n, 1)
  pivot_cols = Vector{Int}()
  r = 0
  last_pivot = 0
  for i = 1:m
    for j = n:-1:last_pivot + 1
      if iszero(A[i, j])
        continue
      end
      x[j, 1] = b[i, 1]
      for k = 1:r
        x[j, 1] -= A[i, pivot_cols[k]]*x[pivot_cols[k], 1]
      end
      x[j, 1] *= inv(A[i, j])
      last_pivot = j
      r += 1
      push!(pivot_cols, j)
      break
    end
  end
  return x
end
    
end #if-false block.
