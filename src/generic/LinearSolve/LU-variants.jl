
function solve_fflu(A::MatElem{T}, b::MatElem{T}) where {T <: RingElement}
   base_ring(A) != base_ring(b) && error("Base rings don't match in solve_fflu")
   nrows(A) != nrows(b) && error("Dimensions don't match in solve_fflu")
   FFLU = deepcopy(A)
   p = PermGroup(nrows(A))()
   r, d = fflu!(p, FFLU)
    #r < nrows(A) && error("Singular matrix in solve_fflu")
   return solve_fflu_precomp(p, FFLU, b), d
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

