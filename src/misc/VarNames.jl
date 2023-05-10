"""
    variable_names(a) -> Vector{Symbol}

Create proper variable names from `a`.

#Examples

```jldoctest
julia> variable_names([:x, :y])
2-element Vector{Symbol}
 :x
 :y

julia> variable_names(:x => 3)
3-element Vector{Symbol}
 Symbol("x[1]")
 Symbol("x[2]")
 Symbol("x[3]")

julia> variable_names(:x => (3,2))
6-element Vector{Symbol}
 Symbol("x[1,1]")
 Symbol("x[2,1]")
 Symbol("x[3,1]")
 Symbol("x[1,2]")
 Symbol("x[2,2]")
 Symbol("x[3,2]")

julia> vars = variable_names(Dict(:x => (1,2), :y => 2, :z => missing))
5-element Vector{Symbol}
 Symbol("x[1,1]")
 Symbol("x[1,2]")
 Symbol("y[1]")
 Symbol("y[2]")
 :z

julia> vars == variable_names((x=(1,2), y=2, z=missing))
true

julia> vars == variable_names([:x => (1,2), :y => 2, :z])
true

```
"""
variable_names(a) = [x for sn in _pairs(a) for x in _variable_names(sn)]

_variable_names(s::VarName) = [Symbol(s)]
_variable_names((s, a)::Pair{<:VarName}) = Symbol.(s, _variable_suffix(a))

_variable_suffix(::Missing) = [""]
_variable_suffix(n::Int) = ["[$i]" for i in 1:n]
_variable_suffix(ns::Tuple{Vararg{Int}}) = ["[$(join(is, ','))]" for is in Iterators.product(Base.OneTo.(ns)...)]

"""
    _pairs(a)

Return iterator of pairs, if possible, otherwise keep `a`.
"""
function _pairs(a)
    try
        return pairs(a)
    catch e
        e isa MethodError && return a
        rethrow()
    end
end

_pairs(a::Pair) = [a]
_pairs(a::Vector) = a
