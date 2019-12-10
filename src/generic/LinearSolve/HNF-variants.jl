
################################################################################
#
#  "Can solve with kernel (Field)". (To be demolished)
#
################################################################################


@doc Markdown.doc"""
    can_solve(A::MatElem{T}, B::MatElem{T}; side = :right) where T <: FieldElem -> Bool, MatElem

Tries to solve $Ax = B$ for $x$ if `side = :right` and $xA = B$ if `side =
:left`.
"""
function can_solve(A::MatElem{T}, B::MatElem{T};
                                  side = :right) where T <: FieldElem
  @assert base_ring(A) == base_ring(B)

  if side === :right
    @assert nrows(A) == nrows(B)
    return _can_solve(A, B)
  elseif side === :left
    @assert ncols(A) == ncols(B)
    b, C = _can_solve(transpose(A), transpose(B))
    if b
      return true, transpose(C)
    else
      return false, C
    end
  else
    error("Unsupported argument :$side for side: Must be :left or :right")
  end
end

function _can_solve(A::MatElem{T}, B::MatElem{T}) where T <: FieldElem
  R = base_ring(A)
  mu = [A B]
  rk, mu = rref(mu) # rref call.
  p = find_pivot(mu)
  if any(i -> i > ncols(A), p)
    return false, B
  end
  sol = zero_matrix(R, ncols(A), ncols(B))
  for i = 1:length(p)
    for j = 1:ncols(B)
      sol[p[i], j] = mu[i, ncols(A) + j]
    end
  end
  return true, sol
end


################################################################################
#
#  solve_hnf
#
################################################################################

@doc Markdown.doc"""
    solve_hnf(A::MatElem{T}, b::MatElem{T}, side = :right) where T <: RingElem -> Bool, MatElem
    
Tries to solve $Ax = b$ for $x$ if `side = :right` or $Ax = B$ if `side = :left`.
The method `hnf_with_transform` must be implemented for the type, and it
is assumed that it returns a correct Hermite normal form. 
"""
function solve_hnf(A::MatElem{T}, b::MatElem{T};
                                  side = :right, kernel=Val(false)) where T <: RingElement

    if side === :left
        # It might be easier to do things this way, but you could in theory allocate
        # the solution space as a TransposeIndexDual. Just remember to return the
        # correct type afterward.
        Adual = TransposeIndexDual(A)
        bdual = TransposeIndexDual(b)
        actual_sol = transpose(solve_hnf(Adual, bdual, side=:right))
        return actual_sol
    end

    if kernel != Val(false)
        error("Kernel not yet supported in `solve_lu`.")
    end
    
    check_solve_instance_is_well_defined(A,b)
    isempty(b) && return _solve_empty(A,b)
    
    # Writing `A^T = UH`, we want to solve Ax = (UH)^T*x = b. To optimize the cubic
    # part of the algorithm, we enforce that access is column-major. Note `H^T` is actually
    # lower-triangular.
    HNFt, g = hnf_with_transform(column_major_access_form(TransposeIndexDual(A)))
    lHNF = TransposeIndexDual(HNFt)
    
    #@info "" lHNF typeof(g)
    
    n = nrows(lHNF)
    m = ncols(lHNF)
    nvecs = ncols(b)
    
    R = base_ring(lHNF)

    prows  = _ut_pivot_columns(HNFt)
    rk = length(prows)

    # Also allocates space for `x`. Since it is unclear to me if `g` has special properties
    # for tall/long matrices, this is the extent to which I can optimize.
    x = zero_matrix(R, m, nvecs)
    
    for k=1:nvecs
        # PROOF OF CORRECT USAGE: 
        # As the elements of `x` are allocated in `zero_matrix`, we see they
        # are distinct and share no references between each other, `A`, or `y` aside from
        # parents. Thus, we have fulfilled the CONTRACT for using the `!!`
        # method.
        if !iszero(rk)
            lHNFview = view(lHNF, prows, :)
            xview = view(x, :, k:k)
            bview = view(b, prows, k:k)

            #@info "" lHNFview xview bview
            
            # TODO: This can be optimized to avoid the allocation, but it requires passing
            # incongruent views to a dangerous function. For now, we keep it simple.

            z = zero_matrix(R, length(prows), 1)
            z = _solve_nonsingular_lt!!_I_agree_to_the_terms_and_conditions_of_this_function(z, lHNFview, bview, rk, 1)

            # Since the setindex is not quite as flexible as with AbstractArrays, we do this
            # manually.
            for ell = 1:length(prows)
                x[ell, k] = z[ell,1]
            end

        end
    end

    # Because `divexact` **Does not throw and error** if an invalid division is performed,
    # we need to perform the consistency check on every row.

    sol = transpose(g)*x
    check_system_is_consistent(A, sol, b, 0)
    
    return sol

    ####
    # With kernel
    # TODO: XXX: Implement the logic.
    
    # Scan backward for the first non-zero row, once found, set the nullspace to
    # be the bottom rows (for a left-solve)
    for i = nrows(HNF):-1:1
        for j = 1:ncols(HNF)
            if !iszero(HNF[i,j])
                N = similar(A, ncols(A), nrows(HNF) - i)
                for k = 1:nrows(N)
                    for l = 1:ncols(N)
                        N[k,l] = T[nrows(T) - l + 1, k]
                    end
                end
                return true, transpose(z*T), N
            end
        end
    end
    N =  similar(A, ncols(A), 0)
end


################################################################################
#
#  "Can solve (Ring)". (To be called HNF solve.)
#
################################################################################


@doc Markdown.doc"""
    can_solve(A::MatElem{T}, B::MatElem{T}, side = :right) where T <: RingElem -> Bool, MatElem
    
Tries to solve $Ax = B$ for $x$ if `side = :right` or $Ax = B$ if `side = :left`
over a euclidean ring.
"""
function can_solve(A::MatElem{T}, B::MatElem{T};
                                  side = :right) where T <: RingElem

    # Just "dual" logic.
    @assert base_ring(A) == base_ring(B)

    if side === :right
        @assert nrows(A) == nrows(B)
        return _can_solve(A, B)
    elseif side === :left
        @assert ncols(A) == ncols(B)
        b, C = _can_solve(transpose(A), transpose(B))
        if b
            return true, transpose(C)
        else
            return false, C
        end
    else
        error("Unsupported argument :$side for side: Must be :left or :right")
    end
end

function _can_solve(a::MatElem{S}, b::MatElem{S}, side = :left) where S <: RingElem
    H, T = hnf_with_transform(transpose(a))
    b = deepcopy(b)
    z = similar(a, ncols(b), ncols(a))
    l = min(nrows(a), ncols(a))

    for i = 1:ncols(b)
        for j = 1:l

            # Pivot location.
            k = 1
            while k <= ncols(H) && iszero(H[j, k])
                k += 1
            end
            if k > ncols(H)
                continue
            end

            # Really, this is divexact with a failure catch.
            q, r = divrem(b[k, i], H[j, k]) 
            if !iszero(r)
                return false, b
            end

            # Upper-triangular substitution loop.
            for h = k:ncols(H)
                b[h, i] -= q*H[j, h]
            end
            z[i, j] = q
        end
    end

    # If the `b` 
    if !iszero(b)
        return false, b
    end
    return true, transpose(z*T)
end


################################################################################
#
#  "Can solve with kernel (Field)" (To be demolished.)
#
################################################################################

@doc Markdown.doc"""
    can_solve_with_kernel(A::MatElem{T}, B::MatElem{T}; side = :right) where T <: FieldElem -> Bool, MatElem, MatElem

Tries to solve $Ax = B$ for $x$ if `side = :right` or $Ax = B$ if `side = :left`.
It returns the solution and the right respectively left kernel of $A$.
"""
function can_solve_with_kernel(A::MatElem{T}, B::MatElem{T}; side = :right) where T <: FieldElem
  @assert base_ring(A) == base_ring(B)
  if side === :right
    @assert nrows(A) == nrows(B)
    return _can_solve_with_kernel(A, B)
  elseif side === :left
    b, C, K = _can_solve_with_kernel(transpose(A), transpose(B))
    @assert ncols(A) == ncols(B)
    if b
      return b, transpose(C), transpose(K)
    else
      return b, C, K
    end
  else
    error("Unsupported argument :$side for side: Must be :left or :right")
  end
end

function _can_solve_with_kernel(A::MatElem{T}, B::MatElem{T}) where T <: FieldElem
  R = base_ring(A)
  mu = [A B]
  rk, mu = rref(mu)
  p = find_pivot(mu)
  if any(i->i>ncols(A), p)
    return false, B, B
  end
  sol = zero_matrix(R, ncols(A), ncols(B))
  for i = 1:length(p)
    for j = 1:ncols(B)
      sol[p[i], j] = mu[i, ncols(A) + j]
    end
  end
  n = zero_matrix(R, ncols(A), ncols(A) - length(p))
  np = sort(setdiff(1:ncols(A), p))
  i = 0
  push!(p, ncols(A)+1)
  for j = 1:length(np)
    if np[j] >= p[i+1]
      i += 1
    end
    if i > 0
      n[p[i], j] = -mu[i, np[j]]
    end
    n[np[j], j] = 1
  end
  return true, sol, n
end
