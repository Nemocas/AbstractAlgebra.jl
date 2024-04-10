# Remember: There are more deprecations in Generic.jl

###############################################################################
#
#   Aliases
#
###############################################################################

# ALL aliases here are only a temporary measure to allow for a smooth transition downstream.
# they will be replaced by deprecations eventually

# added during 0.40.*
@alias FreeModule free_module
@alias VectorSpace vector_space
@alias UniversalPolynomialRing universal_polynomial_ring

@alias is_finiteorder is_finite_order


###############################################################################
#
#   Deprecated bindings
#
###############################################################################

# Deprecated bindings don't get reexported automatically in Hecke/Oscar/etc.
# By calling this macro from the respective packages, we can ensure that the deprecated bindings are available there.
macro include_deprecated_bindings()
    return esc(quote
        #= currently none =#

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
