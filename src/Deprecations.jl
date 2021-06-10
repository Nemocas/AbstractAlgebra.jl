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

# Deprecated in 0.14.*

@deprecate powmod(x, y, z) powermod(x, y, z)

# Deprecated in 0.15.*

@deprecate trail(f) trailing_coefficient(f)

@deprecate lead(f) leading_coefficient(f)

@deprecate lc(f) leading_coefficient(f)

@deprecate lm(f) leading_monomial(f)

@deprecate lt(f) leading_term(f)

@deprecate valence(f::Generic.PolyElem) trailing_coefficient(f)

@deprecate coeffs(a::AbstractAlgebra.MPolyElem{T}) where T <: RingElement coefficients(a)

# Deprecated in 0.16.*

@deprecate map_coeffs(g, p::PolyElem; cached=true, parent=Generic._make_parent(g, p, cached)) map_coefficients(g, p; cached=cached, parent=parent)

@deprecate map_coeffs(f, p::MPolyElem; cached = true, parent = Generic._change_mpoly_ring(parent(f(zero(base_ring(p)))), parent(p), cached)) map_coefficients(f, p; cached = cached, parent = parent)

# Deprecated in 0.18.*

@deprecate involves_at_most_one_variable(p::AbstractAlgebra.MPolyElem) isunivariate(p)

# Deprecated in 0.18.*

@deprecate  base_ring(R::PolyRing{T}) where T <: RingElement coefficient_ring(R)

@deprecate  change_base_ring(R::Ring, p::PolyElem{T}; cached::Bool = true, parent::PolyRing = _change_poly_ring(R, parent(p), cached)) where T <: RingElement change_coefficient_ring(R, p, cached = cached, parent = parent)

