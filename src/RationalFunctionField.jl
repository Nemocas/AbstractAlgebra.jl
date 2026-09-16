###############################################################################
#
#   rational_function_field.jl : Rational function fields
#
###############################################################################

###############################################################################
#
#   rational_function_field constructor
#
###############################################################################

@doc raw"""
    rational_function_field(k::Field, s::VarName=:t; cached::Bool=true)
    rational_function_field(k::Field, varnames::Vector{Symbol}; cached::Bool=true)
    rational_function_field(k::Field, n::Int, s::VarName=:x; cached::Bool=true)

Given a coefficient field `k` and variable names, return a tuple `(S, x)`
consisting of the rational function field $S = k(x, \dots)$ over `k` and its
generator(s) `x`. The numbered form uses the names `Symbol.(s, 1:n)`.

By default (`cached=true`), the output `S` will be cached, i.e. if
`rational_function_field` is invoked again with the same arguments, the same
(*identical*) field is returned. Setting `cached` to `false` ensures a distinct
new field is returned, and will also prevent it from being cached.

For the many ways to specify `varnames` refer to [`polynomial_ring`](@ref) or
[`AbstractAlgebra.@varnames_interface`](@ref).

# Examples

```jldoctest
julia> S, x = rational_function_field(QQ, :x)
(Rational function field over rationals, x)
```
"""
function rational_function_field(k::Field, s::VarName = :t; cached::Bool=true)
   return Generic.rational_function_field(k, Symbol(s); cached=cached)
end

@varnames_interface Generic.rational_function_field(K::Field, s)
