
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
function solve_fflu(A::MatElem{T}, b::MatElem{T}; den=Val(true)) where {T <: RingElement}
    check_solve_instance_is_well_defined(A,b)

    FFLU = deepcopy(A)
    p = PermGroup(nrows(A))()
    r, d = fflu!(p, FFLU)
    #r < nrows(A) && error("Singular matrix in solve_fflu")

    x = _solve_fflu_postcomp(p, FFLU, b)
    
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
end


function solve_lu(A::MatElem{T}, b::MatElem{T}; side=:right) where {T <: FieldElement}

    if side === :left
        # It might be easier to do things this way, but you could in theory allocate
        # the solution space as a TransposeIndexDual. Just remember to return the
        # correct type afterward.
        Adual = TransposeIndexDual(A)
        bdual = TransposeIndexDual(b)
        actual_sol = transpose(solve_lu(Adual, bdual, side=:right))
        return actual_sol
    end

    check_solve_instance_is_well_defined(A,b)
    
    LU = deepcopy(A)
    p = PermGroup(nrows(A))()
    r = lu!(p, LU)

    return solve_lu_precomp(p, LU, b)
end

