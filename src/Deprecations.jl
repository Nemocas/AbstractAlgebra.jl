# Remember: There are more deprecations in Generic.jl

###############################################################################
#
#   Aliases
#
###############################################################################

# ALL aliases here are only a temporary measure to allow for a smooth transition downstream.
# they will be replaced by deprecations eventually

# renamed in 0.45.1
@alias FactoredFractionField factored_fraction_field

#renamed in 0.47.4
@alias dense_poly_type poly_type
@alias dense_poly_ring_type poly_ring_type

###############################################################################
#
#   Deprecated bindings
#
###############################################################################

# Deprecated bindings don't get reexported automatically in Hecke/Oscar/etc.
# By calling this macro from the respective packages, we can ensure that the deprecated bindings are available there.
macro include_deprecated_bindings()
    return esc(quote
        # deprecated in 0.42.0
        Base.@deprecate_binding FreeModule free_module
        Base.@deprecate_binding VectorSpace vector_space
        Base.@deprecate_binding UniversalPolynomialRing universal_polynomial_ring
        Base.@deprecate_binding is_finiteorder is_finite_order

        # renamed in 0.43.0
        Base.@deprecate_binding FreeAssAlgebra FreeAssociativeAlgebra
        Base.@deprecate_binding FreeAssAlgElem FreeAssociativeAlgebraElem
        Base.@deprecate_binding addeq! add!
    end)
end

@include_deprecated_bindings()

###############################################################################
#
#   Deprecations
#
###############################################################################

# deprecated during 0.40.*
@deprecate force_coerce(a, b, throw_error::Type{Val{T}}) where {T} force_coerce(a, b, Val(T))
@deprecate force_op(op::Function, throw_error::Type{Val{T}}, a...) where {T} force_op(op, Val(T), a...)
@deprecate _hnf_minors!(H::MatrixElem{T}, U::MatrixElem{T}, with_transform::Type{Val{S}}) where {T <: RingElement, S} _hnf_minors!(H, U, Val(S))
@deprecate _hnf_kb(A, trafo::Type{Val{T}}) where T _hnf_kb(A, Val(T))
@deprecate _snf_kb(A::MatrixElem{T}, trafo::Type{Val{V}}) where {V, T <: RingElement} _snf_kb(A, Val(V))
@deprecate _weak_popov(A::MatElem{T}, trafo::Type{Val{S}}) where {T <: PolyRingElem, S} _weak_popov(A, Val(S))
@deprecate _extended_weak_popov(A::MatElem{T}, V::MatElem{T}, trafo::Type{Val{S}}) where {T <: PolyRingElem, S} _extended_weak_popov(A, V, Val(S))
@deprecate _popov(A::MatElem{T}, trafo::Type{Val{S}}) where {T <: PolyRingElem, S} _popov(A, Val(S))
@deprecate _hnf_via_popov(A::MatElem{T}, trafo::Type{Val{S}}) where {T <: PolyRingElem, S} _hnf_via_popov(A, Val(S))
@deprecate gen(a::MPolyRing{T}, i::Int, ::Type{Val{ord}}) where {T <: RingElement, ord} gen(a, i, Val(ord))

# deprecated in 0.43.0
@deprecate change_base_ring(p::MPolyRingElem{T}, g, new_polynomial_ring) where {T<:RingElement} map_coefficients(g, p, parent = new_polynomial_ring)
@deprecate mulmod(a::S, b::S, mod::Vector{S}) where {S <: MPolyRingElem} Base.divrem(a * b, mod)[2]
@deprecate var"@attr"(__source__::LineNumberNode, __module__::Base.Module, expr::Expr) var"@attr"(__source__, __module__, :Any, expr) # delegate `@attr functionexpression` to `@attr Any functionexpression` (macros are just functions with this weird extra syntax)

# deprecated in 0.45.0
@deprecate test_iterate ConformanceTests.test_iterate

# deprecated in 0.45.2
import Base: contains
@deprecate contains(I::Ideal{T}, J::Ideal{T}) where {T <: RingElement} issubset(J, I)

# to be deprecated in the next breaking release
zero_matrix(::Type{Int}, r, c) = zeros(Int, r, c)
Array(R::NCRing, r::Int...) = Array{elem_type(R)}(undef, r)
function zeros(R::NCRing, r::Int...)
    #Base.depwarn("`zeros(R::NCRing, r::Int...)` is deprecated, use `zero_matrix(R, r...)` instead.", :zeros)
    T = elem_type(R)
    A = Array{T}(undef, r)
    for i in eachindex(A)
        A[i] = R()
    end
    return A
end

# to be deprecated in a future minor release
PolyRing(R::NCRing) = poly_ring(R, :x)

polynomial_ring_only(R::T, s::Symbol; cached::Bool=true) where T<:NCRing =
   poly_ring(R, s; cached)
polynomial_ring_only(R::T, s::Vector{Symbol}; internal_ordering::Symbol=:lex, cached::Bool=true) where T<:Ring =
   poly_ring(R, s; internal_ordering, cached)
