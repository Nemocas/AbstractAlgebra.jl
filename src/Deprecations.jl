# Remember: There are more deprecations in Generic.jl

###############################################################################
#
#   Aliases
#
###############################################################################

# ALL aliases here are only a temporary measure to allow for a smooth transition downstream.
# they will be replaced by deprecations eventually

# renamed in 0.42.3
@alias FreeAssAlgebra FreeAssociativeAlgebra
@alias FreeAssAlgElem FreeAssociativeAlgebraElem

# renamed in 0.42.6
@alias addeq! add!

###############################################################################
#
#   Deprecated bindings
#
###############################################################################

# Deprecated bindings don't get reexported automatically in Hecke/Oscar/etc.
# By calling this macro from the respective packages, we can ensure that the deprecated bindings are available there.
macro include_deprecated_bindings()
    return esc(quote
        # deprecated in 0.42
        Base.@deprecate_binding FreeModule free_module
        Base.@deprecate_binding VectorSpace vector_space
        Base.@deprecate_binding UniversalPolynomialRing universal_polynomial_ring
        Base.@deprecate_binding is_finiteorder is_finite_order

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
import .Generic: exponent_vector; @deprecate exponent_vector(a::Generic.MPoly{T}, i::Int, ::Type{Val{ord}}) where {T <: RingElement, ord} exponent_vector(a, i, Val(ord))
import .Generic: exponent; @deprecate exponent(a::Generic.MPoly{T}, i::Int, j::Int, ::Type{Val{ord}}) where {T <: RingElement, ord} exponent(a, i, j, Val(ord))
import .Generic: set_exponent_vector!; @deprecate set_exponent_vector!(a::Generic.MPoly{T}, i::Int, exps::Vector{Int}, ::Type{Val{ord}}) where {T <: RingElement, ord} set_exponent_vector!(a, i, exps, Val(ord))
import .Generic: is_gen; @deprecate is_gen(x::Generic.MPoly{T}, ::Type{Val{ord}}) where {T <: RingElement, ord} is_gen(x, Val(ord))
import .Generic: degree; @deprecate degree(f::Generic.MPoly{T}, i::Int, ::Type{Val{ord}}) where {T <: RingElement, ord} degree(f, i, Val(ord))

# will become deprecated in 0.43.0
change_base_ring(p::MPolyRingElem{T}, g, new_polynomial_ring) where {T<:RingElement} = map_coefficients(g, p, parent = new_polynomial_ring)
mulmod(a::S, b::S, mod::Vector{S}) where {S <: MPolyRingElem} = Base.divrem(a * b, mod)[2]
var"@attr"(__source__::LineNumberNode, __module__::Base.Module, expr::Expr) = var"@attr"(__source__, __module__, :Any, expr) # delegate `@attr functionexpression` to `@attr Any functionexpression` (macros are just functions with this weird extra syntax)
