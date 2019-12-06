
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
    base_ring(A) != base_ring(b) && error("Base rings don't match in solve_fflu")
    nrows(A) != nrows(b) && error("Dimensions don't match in solve_fflu")
    FFLU = deepcopy(A)
    p = PermGroup(nrows(A))()
    r, d = fflu!(p, FFLU)
    #r < nrows(A) && error("Singular matrix in solve_fflu")
    #return solve_fflu_precomp(p, FFLU, b), d
    if den == Val(true) || den == Val(:true)
        return solve_fflu_precomp(p, FFLU, b), d

    elseif den == Val(false) || den == Val(:false)
        return solve_fflu_precomp(p, FFLU, b)

    elseif den == Val(:det)
        x = solve_fflu_precomp(p, FFLU)
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


function solve_lu(A::MatElem{T}, b::MatElem{T}) where {T <: FieldElement}
    base_ring(A) != base_ring(b) && error("Base rings don't match in solve_lu")
    nrows(A) != nrows(b) && throw(DimensionMismatch("nrows(A) != nrows(b)"))

    if nrows(A) == 0 || ncols(A) == 0 || ncols(b) == 0
        throw(DimensionMismatch("Solve Ax=b instance with either A=$A or b=$b empty."))
    end

    LU = deepcopy(A)
    p = PermGroup(nrows(A))()
    r = lu!(p, LU)
    #r < nrows(A) && error("Singular matrix in solve_lu")
    return solve_lu_precomp(p, LU, b)
end

