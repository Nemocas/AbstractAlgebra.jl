# Deprecated in 0.11.*

@deprecate rref(a::MatrixElem{<:RingElement}) rref_rational(a)

@deprecate rref!(a::MatrixElem{<:RingElement}) rref_rational!(a)

# Deprecated in 0.12.*

@deprecate MatrixSpace(R::Ring, n::Int, m::Int, cached::Bool) MatrixSpace(R, n, m, cached = cached)

function (G::Generic.SymmetricGroup)()
    Base.depwarn("(::SymmetricGroup)() to get the group identity is deprecated, use one(::SymmetricGroup) instead", :one)
    return one(G)
end

# Deprecated in 0.15.*

@deprecate lead(f) leading_coefficient(f)

@deprecate coeffs(a::AbstractAlgebra.MPolyElem{T}) where T <: RingElement coefficients(a)
