
function solve_left(a::AbstractAlgebra.MatElem{S}, b::AbstractAlgebra.MatElem{S}) where S <: RingElement
    @warn "Function has been depreciated. Use `solve_hnf(a,b,side=:left)` instead." maxlog=1
    return solve_hnf(a,b,side=:left)
end

function solve_left(a::AbstractAlgebra.MatElem{S}, b::AbstractAlgebra.MatElem{S}) where S <: FieldElement
    @warn "Function has been depreciated. Use `solve_lu(a,b,side=:left)` instead." maxlog=1
    return solve_lu(a,b,side=:left)
end


# Find the pivot columns of an rref matrix
function find_pivot(A::AbstractAlgebra.MatElem{T}) where T <: RingElement
  p = Int[]
  j = 0
  for i = 1:nrows(A)
    j += 1
    if j > ncols(A)
      return p
    end
    while iszero(A[i, j])
      j += 1
      if j > ncols(A)
        return p
      end
    end
    push!(p, j)
  end
  return p
end
