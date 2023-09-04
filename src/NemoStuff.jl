export neg!

sub!(z::Rational{Int}, x::Rational{Int}, y::Int) = x - y

neg!(z::Rational{Int}, x::Rational{Int}) = -x

add!(z::Rational{Int}, x::Rational{Int}, y::Int) = x + y

mul!(z::Rational{Int}, x::Rational{Int}, y::Int) = x * y

################################################################################
#
#  Diagonal (block) matrix creation
#
################################################################################

@doc raw"""
    diagonal_matrix(x::T...) where T <: RingElem -> MatElem{T}
    diagonal_matrix(x::Vector{T}) where T <: RingElem -> MatElem{T}
    diagonal_matrix(Q, x::Vector{T}) where T <: RingElem -> MatElem{T}

Returns a diagonal matrix whose diagonal entries are the elements of $x$.

# Examples

```jldoctest
julia> diagonal_matrix(QQ(1), QQ(2))
[1   0]
[0   2]

julia> diagonal_matrix([QQ(3), QQ(4)])
[3   0]
[0   4]

julia> diagonal_matrix(QQ, [5, 6])
[5   0]
[0   6]
```
"""
function diagonal_matrix(R::Ring, x::Vector{<:RingElement})
    x = R.(x)
    M = zero_matrix(R, length(x), length(x))
    for i = 1:length(x)
        M[i, i] = x[i]
    end
    return M
end

function diagonal_matrix(x::T, xs::T...) where {T<:RingElem}
    return diagonal_matrix(collect((x, xs...)))
end

diagonal_matrix(x::Vector{<:RingElement}) = diagonal_matrix(parent(x[1]), x)

@doc raw"""
    diagonal_matrix(x::Vector{T}) where T <: MatElem -> MatElem

Returns a block diagonal matrix whose diagonal blocks are the matrices in $x$.
"""
function diagonal_matrix(x::Vector{T}) where {T<:MatElem}
    return cat(x..., dims=(1, 2))::T
end

function diagonal_matrix(x::T, xs::T...) where {T<:MatElem}
    return cat(x, xs..., dims=(1, 2))::T
end

function diagonal_matrix(R::Ring, x::Vector{<:MatElem})
    if length(x) == 0
        return zero_matrix(R, 0, 0)
    end
    x = [change_base_ring(R, i) for i in x]
    return diagonal_matrix(x)
end

base_ring(::Vector{Int}) = Int

function is_symmetric(M::MatElem)
    for i in 1:nrows(M)
        for j in i:ncols(M)
            if M[i, j] != M[j, i]
                return false
            end
        end
    end
    return true
end

zero_matrix(::Type{Int}, r, c) = zeros(Int, r, c)

#TODO: should be done in Nemo/AbstractAlgebra s.w.
#      needed by ^ (the generic power in Base using square and multiply)
Base.copy(f::Generic.MPoly) = deepcopy(f)
Base.copy(f::Generic.Poly) = deepcopy(f)

################################################################################
#
#  Minpoly and Charpoly
#
################################################################################

function minpoly(M::MatElem)
    k = base_ring(M)
    kx, x = polynomial_ring(k, cached=false)
    return minpoly(kx, M)
end

function charpoly(M::MatElem)
    k = base_ring(M)
    kx, x = polynomial_ring(k, cached=false)
    return charpoly(kx, M)
end

###############################################################################
#
#  Sub
#
###############################################################################

function sub(M::MatElem, rows::Vector{Int}, cols::Vector{Int})
    N = zero_matrix(base_ring(M), length(rows), length(cols))
    for i = 1:length(rows)
        for j = 1:length(cols)
            N[i, j] = M[rows[i], cols[j]]
        end
    end
    return N
end

function sub(M::MatElem{T}, r::AbstractUnitRange{<:Integer}, c::AbstractUnitRange{<:Integer}) where {T}
    z = similar(M, length(r), length(c))
    for i in 1:length(r)
        for j in 1:length(c)
            z[i, j] = M[r[i], c[j]]
        end
    end
    return z
end

#trivia to make life easier

gens(L::SimpleNumField{T}) where {T} = [gen(L)]

function gen(L::SimpleNumField{T}, i::Int) where {T}
    i == 1 || error("index must be 1")
    return gen(L)
end

function Base.getindex(L::SimpleNumField{T}, i::Int) where {T}
    if i == 0
        return one(L)
    elseif i == 1
        return gen(L)
    else
        error("index has to be 0 or 1")
    end
end

ngens(L::SimpleNumField{T}) where {T} = 1

is_unit(a::NumFieldElem) = !iszero(a)

canonical_unit(a::NumFieldElem) = a
