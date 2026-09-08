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

Given a coefficient field `k` and variable names, return a tuple `(S, x)`
consisting of the rational function field $S = k(x, \dots)$ over `k` and its
generator(s) `x`.

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

"""
    @rational_function_field(k::Field, varnames...; cached=true)

Return the field from [`rational_function_field`](@ref) and introduce the
generators into the current scope.

# Examples

```jldoctest
julia> S = @rational_function_field(QQ, [:s, :t])
Rational function field
  over rationals

julia> (s + t)//t
(s + t)//t
```
"""
:(@rational_function_field)
