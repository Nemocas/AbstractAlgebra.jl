# Deprecated in 0.9.*

@deprecate set_prec!(a::AbstractAlgebra.Generic.LaurentSeriesElem, prec::Int) set_precision!(a, prec) 

@deprecate set_prec!(a::AbstractAlgebra.SeriesElem, prec::Int) set_precision!(a, prec)

@deprecate set_val!(a::AbstractAlgebra.Generic.LaurentSeriesElem, val::Int) set_valuation!(a, val)

@deprecate set_val!(a::AbstractAlgebra.RelSeriesElem, val::Int) set_valuation!(a, val)

# Deprecated in 0.11.*

@deprecate rref(a::MatrixElem{<:RingElement}) rref_rational(a)

@deprecate rref!(a::MatrixElem{<:RingElement}) rref_rational!(a)

# Deprecated in 0.12.*

@deprecate MatrixSpace(R::Ring, n::Int, m::Int, cached::Bool) MatrixSpace(R, n, m, cached = cached)

@deprecate MatrixAlgebra(R::Ring, n::Int, cached::Bool) MatrixAlgebra(R, n, cached = cached)

function (G::Generic.SymmetricGroup)()
    Base.depwarn("(::SymmetricGroup)() to get the group identity is deprecated, use one(::SymmetricGroup) instead", :one)
    return one(G)
end
