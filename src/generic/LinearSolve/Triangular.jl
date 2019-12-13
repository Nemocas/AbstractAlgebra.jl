

# These functions are honestly just triangular solving in disguise.

function _solve_fflu_postcomp(p::Generic.Perm, FFLU::MatElem{T}, b::MatElem{T}) where {T <: RingElement}

    # TODO: Implement the new FFLU code & debug.
    @warn "fflu does not work properly for `long` matrices. Caution is advised." maxlog=1
    @warn "fflu does not fully reduce matrix in singular case. Some care is required." maxlog=1

    n = nrows(FFLU)
    m = ncols(FFLU)
    ncolsb = ncols(b)
    
    rk = begin
        rk=min(n,m)
        for i = min(n,m):-1:1
            iszero(FFLU[i,i]) ? rk -= 1 : nothing
        end
        rk
    end

    R = base_ring(FFLU)

    # TODO: In general, `zero_matrix` will return a type dependent only on the base_ring,
    # which could in theory be a different concrete type from the input. It will be crucial
    # to implement a zero-matrix constructor that takes into account the type of `LU`.
    #
    # This breaks at least one test...
    x  = zero_matrix(R, m, ncolsb)
    pb = p*b

    t = base_ring(b)()
    s = base_ring(b)()
    ZERO = zero(R)

    # For each column of x
    for k in 1:ncolsb
        # In the first backsolve step, Only the part `y[1:rk, :]` is needed to
        # construct the solution to the linear equation `(LD^-1)^(-1) b =: y = Ux`.
        # or to check that the system is consistent. Essentially, we are restricting
        # to the principal (rk x rk)-subblock of `U`.
        #
        # Populate the initial values in `x`. There is surely a better way to do this.
        for i=1:rk
            x[i, k] = deepcopy(pb[i, k])
        end

        # Backsolve the lower triangular part. However, things have been rearranged
        # so that a single column is accessed in the inner (j) loop, rather than a row.
        #
        # The trick here is that the atomic lower triangular operations
        # almost commute with `D`, and that the partial quotients `D[i,i]/D[i-1,i-1]`
        # are integral and telescope.
        #
        # Writing `L = L_1 L_2 .. L_n` as a product of atomic lower triangular matrices
        # (See https://en.wikipedia.org/wiki/Triangular_matrix#Atomic_triangular_matrix)
        # Each loop iteration corresponds to solving ...something... .
        #
        for i in 1:(rk - 1)
            t = add!(t, x[i, k], ZERO)
            t = minus!(t)
            for j in (i + 1):rk
                x[j, k] = mul_red!(x[j, k], FFLU[i, i], x[j, k], false)

                s = mul_red!(s, FFLU[j, i], t, false)
                x[j, k] = addeq!(x[j, k], s)
                x[j, k] = reduce!(x[j, k])
                if i > 1
                    x[j, k] = divexact_left(FFLU[i - 1, i - 1], x[j, k])
                end
            end
        end

        # The very last back-solve step is just multiplying `x[rk, k]` by the last pivot
        # of `LU[rk,rk]`, which
        # will be canceled out by the first step of the next part.
        #
        # TODO: There are some clever trick one can implement to get cancellation to happen
        # through the function call.
        x[rk,k] = mul!(x[rk,k], FFLU[rk,rk], x[rk, k])
        
        # Since _ut_pivot_columns only checks entries above the diagonal,
        # it effectively only reads the `U` part of the `LU` object.
        pcols  = _ut_pivot_columns(FFLU)

        for i in (rk - 1):-1:1
            x[i, k] = mul!(x[i, k], FFLU[rk, rk], x[i, k])
        end

        
        # PROOF OF CORRECT USAGE:
        # By construction of `x` using `zero_matrix`, and assignment using `deepcopy`, we see
        # `x` is freshly allocated space and the entries are deepcopyed/newly allocated.
        #
        # After this, entries of `x` are modified only by `mul_red!`, which under correct
        # functioning, does not allow `x` to share a reference with `FFLU` or other operands,
        # aside from itself.
        #
        # Thus, the elements of `x` share no references between each other or the entries of
        # FFLU, even in part, aside from parents.
        # Moreover, trivially, `x===x`. Thus, we have fulfilled the CONTRACT for using
        # the `!!` method.
        if !iszero(rk)
            x[:,k] = _solve_nonsingular_ut!!_I_agree_to_the_terms_and_conditions_of_this_function(x[:,k], FFLU[:, pcols], x[:,k], rk, 1)
        end
        #TODO: Replace the implicit `get_index` with a `view`.
    end
    return x
end

function solve_lu_precomp(p::Generic.Perm, LU::MatElem{T}, b::MatrixElem{T}, rk) where {T <: RingElement}

    n = nrows(LU)
    m = ncols(LU)

    ncolsb = ncols(b)
    R = base_ring(LU)

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
                x[i, k] = addeq!(x[i, k], t)
            end
            x[i, k] = reduce!(x[i, k])
        end

        # Below the last row where `U` has a pivot, the entries
        # of `LU` are just the entries of `L`.
        check_system_is_consistent(LU, view(x, : , k:k), view(pb, :, k:k) , rk)
        
        # Since _ut_pivot_columns only checks entries above the diagonal,
        # it effectively only reads the `U` part of the `LU` object.
        pcols  = _ut_pivot_columns(LU)
        
        # PROOF OF CORRECT USAGE:
        # By use of `deepcopy`, we see `x` is freshly allocated space and the entries are
        # deepcopyed/newly allocated. 
        # Thus, the elements of `x` share no references, even in part, aside from parents.
        # Additionally, we have that `xview === xview`. The number of pivot columns is equal to
        # the rank, so the upper triangular part of LUview is a non-singular square upper
        # triangular matrix.
        # Thus, we have fulfilled the CONTRACT for using the `!!` method.
        if !iszero(rk)
            LUview = view(LU, 1:rk, pcols)
            xview = view(x, 1:rk, k:k)      # The pivot rows are always the same.

            xview = _solve_nonsingular_ut!!_I_agree_to_the_terms_and_conditions_of_this_function(xview, LUview, xview, rk, 1)

            # The entries of `xview` are only assigned to the first `1:rk` rows of `x`. These
            # are not actually the correct since we have ignored the intermediate columns of
            # `LU`. Thus, we need to permute the entries of `x` to the correct positions.
            #
            # Note that the pivot columns are given in ascending order. Entries in `x` below
            # the rank are simply zero by design.
            for ell = length(pcols):-1:1
                x[pcols[ell], k], x[ell, k] = xview[ell,1], x[pcols[ell], k]
            end
        end
    end
    return x
end

###############################################################################
#
#   Internal Nonsingular triangular solving. (The Danger Zone)
#
###############################################################################


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
# 3. The matrix `A` is square, nonsingular, and in reduced lower-triangular form.
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
function _solve_nonsingular_lt!!_I_agree_to_the_terms_and_conditions_of_this_function(x, L, b, rk, k)

    n  = nrows(L)
    m  = ncols(L)
    #rk = min(n,m)::Int

    # Arithmetic container.
    t = base_ring(b)()
    ZERO = base_ring(b)()

    # NOTE: This is the correct order of division for non-commutative rings.
    x[1, k] = divexact_left(L[1, 1], b[1, k])

    # The lower triangular back-substitution. Along rows.
    for i in 2:rk
        
        # Theoretically `x[i,k] = add!(x[i,k], b[i,k], ZERO)` should be sufficient,
        # but I really don't trust that `add!` has uniform behaviour.
        t = zero!(t)
        t = addeq!(t, b[i,k])
        x[i,k] = add!(x[i,k], t, ZERO)
        
        for j in 1:(i - 1)
            # x[i, k] = x[i, k] - x[j, k] * L[i, j]
            t = mul_red!(t, x[j, k], -L[i, j], false)
            x[i, k] = addeq!(x[i, k], t)
        end
        x[i, k] = reduce!(x[i, k])
        x[i, k] = divexact_left(L[i, i], x[i, k])
    end

    return x
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
# 3. The matrix `A` is square, nonsingular, and in reduced upper-triangular form.
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

###############################################################################
#
#   Pivot column selection.
#
###############################################################################

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
