
@doc Markdown.doc"""
    solve_fflu(A::MatElem{T}, b::MatElem{T}; den = Val(true)) where {T <: RingElement}

Gives a solution to the linear equation `A*x = d*b`, with `x,d` defined over `base_ring(A)`. 
The parameter `den` determines if the denominator `d` is returned, and how it is normalized.

`den = Val(false)`   :  No denominator returned.
`den = Val(:false)`

`den = Val(true)`    :  The denonimator from the fflu algorithm is returned,
`den = Val(:true)`      it is the determinant up to a sign.

`den = Val(:det)`    :  The determinant of `A` is returned.

Default is Val(true) for legacy reasons, which is inconsistant with normal solve syntax.
"""
function solve_fflu(A::MatElem{T}, b::MatElem{T};
                    den=Val(true), side = :right, kernel = Val(false)) where {T <: RingElement}

    if side === :left
        # It might be easier to do things this way, but you could in theory allocate
        # the solution space as a TransposeIndexDual. Just remember to return the
        # correct type afterward.
        Adual = TransposeIndexDual(A)
        bdual = TransposeIndexDual(b)
        sol, d = solve_fflu(Adual, bdual, side=:right)
        actual_sol = transpose(sol)
        return actual_sol, d
    end

    if kernel != Val(false)
        error("Kernel not yet supported in `solve_fflu`.")
    end
    
    check_solve_instance_is_well_defined(A,b)
    isempty(b) && return _solve_empty(A,b), base_ring(b)(1) # Empty determinant is 1.
    
    FFLU = deepcopy(A)
    p = PermGroup(nrows(A))()
    r, d = fflu!(p, FFLU)

    if iszero(r)
        x = zero_matrix(base_ring(A), ncols(A), ncols(b))
        check_system_is_consistent(A, x, d*b, r)
        return x, d
    end

    x = _solve_fflu_postcomp(p, FFLU, b)

    # Because of the implicit `D` factor, we have to check consistency at the end.
    check_system_is_consistent(A, x, d*b, r)

    if den == Val(true) || den == Val(:true)
        return x, d

    elseif den == Val(false) || den == Val(:false)
        return x

    elseif den == Val(:det)
        if iszero(parity(p))
            return x, d
        else
            for i = 1:nrows(x)
                for j = 1:ncols(x)
                    x[i,j] = minus!(x[i,j])
                end
            end
            return x, minus!(d)
        end
    end
    error("Option for `den` not accepted.")
end


# Note that for "tall"
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


function solve_lu(A::MatElem{T}, b::MatElem{T};
                  side=:right, kernel=Val(false)) where {T <: FieldElement}

    if side === :left
        # It might be easier to do things this way, but you could in theory allocate
        # the solution space as a TransposeIndexDual. Just remember to return the
        # correct type afterward.
        Adual = TransposeIndexDual(A)
        bdual = TransposeIndexDual(b)
        actual_sol = transpose(solve_lu(Adual, bdual, side=:right, kernel=kernel))
        return actual_sol
    end

    check_solve_instance_is_well_defined(A,b)

    if kernel == Val(true)
        return _losen_lu_mit_kernel(A,b)
    end
        
    isempty(b) && return _solve_empty(A,b)

    # LU decomposition step.
    LU = deepcopy(A)
    p = PermGroup(nrows(A))()
    rk = lu!(p, LU)

    ncolsb = ncols(b)
    R = base_ring(LU)

    # TODO: In general, `zero_matrix` will return a type dependent only on the base_ring,
    # which could in theory be a different concrete type from the input. It will be crucial
    # to implement a zero-matrix constructor that takes into account the type of `LU`.
    x  = zero_matrix(R, ncols(LU), ncolsb)
    pb = p*b

    # Default rank zero behaviour.
    if iszero(rk)
        check_system_is_consistent(LU, x, pb, rk)
        return x
    end
    
    t = base_ring(b)()
    s = base_ring(b)()
    
    # For each column of b, solve Ax=b[:,j].
    for k in 1:ncolsb
        
        # Populate the initial values in `x`. There is surely a better way to do this.
        for i=1:rk
            x[i, k] = deepcopy(pb[i, k])
        end

        LUview = view(LU, 1:rk, 1:rk)
        xview = view(x, 1:rk, k:k)
        xview = _solve_nonsingular_lt!!_I_agree_to_the_terms_and_conditions_of_this_function(xview, LUview, xview, rk, 1, unit_diagonal=Val(true))

        for i=1:rk
            x[i, k] = xview[i, 1]
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

    return x
end

# This function works essentially the same as above, but with the solving done in the
# opposite order.
function _losen_lu_mit_kernel(A,b)

    (isempty(A) || isempty(b)) && error("empty solving with kernel not yet supported.")

    # LU decomposition step.
    LU = transpose(deepcopy(A))
    p = PermGroup(nrows(LU))()
    rk = lu!(p, LU)
    UtLt = transpose(LU)
    
    ncolsb = ncols(b)
    R = base_ring(UtLt)

    # TODO: In general, `zero_matrix` will return a type dependent only on the base_ring,
    # which could in theory be a different concrete type from the input. It will be crucial
    # to implement a zero-matrix constructor that takes into account the type of `UtLt`.
    x  = zero_matrix(R, ncols(UtLt), ncolsb)

    # Default rank zero behaviour.
    if iszero(rk)
        check_system_is_consistent(A, x, b, rk)
        return x, identity_matrix(R, ncols(UtLt))
    end
    
    t = base_ring(b)()
    s = base_ring(b)()
    
    # For each column of b, solve Ax=b[:,j].
    for k in 1:ncolsb

        # _ut_pivot_columns only checks entries above the diagonal.
        prows  = _ut_pivot_columns(LU)
        
        # PROOF OF CORRECT USAGE:
        # By use of `deepcopy`, we see `x` is freshly allocated space and the entries are
        # deepcopyed/newly allocated. 
        # Thus, the elements of `x` share no references, even in part, aside from parents.
        # Additionally, we have that `xview === xview`. The number of pivot columns is equal to
        # the rank, so the upper triangular part of UtLtview is a non-singular square upper
        # triangular matrix.
        # Thus, we have fulfilled the CONTRACT for using the `!!` method.

        # Populate the initial values in `x`. There is surely a better way to do this.
        for i=1:rk
            x[i, k] = deepcopy(b[prows[i], k])
        end

        UtLtview = view(UtLt, prows, 1:rk)
        xview = view(x, 1:rk, k:k)      # The pivot rows are always the same.

        xview = _solve_nonsingular_lt!!_I_agree_to_the_terms_and_conditions_of_this_function(xview, UtLtview, xview, rk, 1)

        for i=1:rk
            x[i, k] = xview[i, 1]
        end        

        #@info "" xview UtLtview
        
        # Below the last row where `U` has a pivot, the entries
        # of `UtLt` are just the entries of `Ut`.
        check_system_is_consistent(UtLt, view(x, : , k:k), view(b, :, k:k) , rk)

        # Backsolve step
        UtLtview = view(UtLt, 1:rk, 1:rk)
        xview = view(x, 1:rk, k:k)
        xview = _solve_nonsingular_ut!!_I_agree_to_the_terms_and_conditions_of_this_function(xview, UtLtview, xview, rk, 1, unit_diagonal=Val(true))

        for i=1:rk
            x[i, k] = xview[i, 1]
        end        
    end

    # Wir machen das kernel hier.
    N = identity_matrix(R, ncols(UtLt))[:, rk+1:ncols(UtLt)]

    #@info "" UtLt N rk typeof(N)

    for k in 1:ncols(UtLt)-rk
        # N is freshly allocated space, so the usual argument applies here
        # to show we have satisfied the CONTRACT.
        N[1:rk, k] = _solve_nonsingular_ut!!_I_agree_to_the_terms_and_conditions_of_this_function(N[1:rk,k], UtLt, -UtLt[1:rk,k+rk], rk, 1, unit_diagonal=Val(true))
    end

    invp = inv(p)
    return invp*x, invp*N
end
