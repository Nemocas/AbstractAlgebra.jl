
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
                                  side = :right) where T <: RingElement

    if side === :left
        # It might be easier to do things this way, but you could in theory allocate
        # the solution space as a TransposeIndexDual. Just remember to return the
        # correct type afterward.
        Adual = TransposeIndexDual(A)
        bdual = TransposeIndexDual(b)
        actual_sol = transpose(solve_hnf(Adual, bdual, side=:right))
        return actual_sol
    end
    
    check_solve_instance_is_well_defined(A,b)

    @warn "Principal `r x r`-block assumed faithful. Please fix!"
    
    HNF, g = hnf_with_transform(A)

    n = nrows(HNF)
    m = ncols(HNF)
    nvecs = ncols(b)
    
    R = base_ring(HNF)

    pcols  = _ut_pivot_columns(HNF)
    rk = length(pcols)

    # Also allocates space for `x`. Since it is unclear to me if `g` has special properties
    # for tall/long matrices, this is the extent to which I can optimize.
    y = g*b
    x = zero_matrix(R, m, nvecs)
    
    for k=1:nvecs
        # PROOF OF CORRECT USAGE: 
        # As the elements of `x` are allocated in `zero_matrix`, we see they
        # are distinct and share no references between each other, `A`, or `y` aside from
        # parents. Thus, we have fulfilled the CONTRACT for using the `!!`
        # method.
        x[:,k] = _solve_nonsingular_ut!!_I_agree_to_the_terms_and_conditions_of_this_function(x[:,k], HNF[:, pcols], y[:,k], rk, 1)
    end

    # Because `divexact` **Does not throw and error** if an invalid division is performed,
    # we need to perform the consistency check on every row.
    check_system_is_consistent(A, x, b, 0)
    
    return x
        
    #=
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
    =#
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



################################################################################
#
#  "Can solve with kernel (Ring)" (To be demolished.)
#
################################################################################


#=
@doc Markdown.doc"""
    can_solve_with_kernel(A::MatElem{T}, B::MatElem{T}) where T <: RingElem -> Bool, MatElem, MatElem

Tries to solve $Ax = B$ for $x$ if `side = :right` or $xA = B$ if `side = :left`.
It returns the solution and the right respectively left kernel of $A$.
"""
function can_solve_with_kernel(A::MatElem{T}, B::MatElem{T}; side = :right) where T <: RingElem
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

function _can_solve_with_kernel(a::MatElem{S}, b::MatElem{S}) where S <: RingElem
  H, T = hnf_with_transform(transpose(a))
  z = similar(a, ncols(b), ncols(a))
  l = min(nrows(a), ncols(a))
  b = deepcopy(b)
  for i=1:ncols(b)
    for j=1:l
      k = 1
      while k <= ncols(H) && iszero(H[j, k])
        k += 1
      end
      if k > ncols(H)
        continue
      end
      q, r = divrem(b[k, i], H[j, k])
      if !iszero(r)
        return false, b, b
      end
      for h=k:ncols(H)
        b[h, i] -= q*H[j, h]
      end
      z[i, j] = q
    end
  end
  if !iszero(b)
    return false, b, b
  end

  for i = nrows(H):-1:1
    for j = 1:ncols(H)
      if !iszero(H[i,j])
        N = similar(a, ncols(a), nrows(H) - i)
        for k = 1:nrows(N)
          for l = 1:ncols(N)
            N[k,l] = T[nrows(T) - l + 1, k]
          end
        end
        return true, transpose(z*T), N
      end
    end
  end
  N =  similar(a, ncols(a), 0)

  return true, transpose(z*T), N
end
=#
