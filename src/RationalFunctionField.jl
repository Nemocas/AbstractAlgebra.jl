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
    rational_function_field(k::Field, varnames::AbstractArray{<:VarName}; cached::Bool=true)
    rational_function_field(k::Field, n::Int, s::VarName=:x; cached::Bool=true)

Given a coefficient field `k` and one variable name `s`, return the rational
function field $S = k(s)$ and its generator. If an array `varnames` is supplied,
return the multivariate rational function field over `k` and an array of
generators with the same shape. The richer variable-name specifications
described for [`AbstractAlgebra.@varnames_interface`](@ref) are also supported.

The numbered form constructs `n` generators whose names have prefix `s`; for
example, `n=3` and `s=:y` create generators `y1`, `y2`, and `y3`. The default
single variable name is `:t`, while the default prefix for numbered variables
is `:x`.

By default (`cached=true`), the output `S` will be cached, i.e. if
`rational_function_field` is invoked again with the same arguments, the same
(*identical*) field is returned. Setting `cached` to `false` ensures a distinct
new field is returned, and will also prevent it from being cached.

# Examples

```jldoctest
julia> S, x = rational_function_field(QQ, :x)
(Rational function field over rationals, x)

julia> S(123)
123

julia> S(numerator(x + 1, false), numerator(x + 2, false))
(x + 1)//(x + 2)

julia> R, (x, y) = rational_function_field(QQ, [:x, :y])
(Rational function field over rationals, AbstractAlgebra.Generic.RationalFunctionFieldElem{Rational{BigInt}, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}[x, y])

julia> (x + y)//y^2
(x + y)//y^2
```
"""
function rational_function_field(k::Field, s::VarName = :t; cached::Bool=true)
   return Generic.rational_function_field(k, Symbol(s); cached=cached)
end

@varnames_interface Generic.rational_function_field(K::Field, s)
