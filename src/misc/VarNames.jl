using Base.Iterators

const VarName = Union{Symbol, AbstractString, Char}
const VarShape = Union{Missing, Int, Tuple{Vararg{Int}}}
# const VarNamesVector = Vector{Union{<:VarName, Pair{<:VarName, <:VarShape}}}
# const VarNames = Union{
#     Pair{<:VarName, <:VarShape},
#     Dict{<:VarName, <:VarShape},
#     NamedTuple{<:Any, Tuple{Vararg{VarShape}}},
#     VarNamesVector
# }

"""
    variable_names(a) -> Vector{Symbol}

Create proper variable names from `a`.

#Examples

```jldoctest
julia> variable_names([:x, :y])
2-element Vector{Symbol}:
 :x
 :y

julia> variable_names(:x => 3)
3-element Vector{Symbol}:
 Symbol("x[1]")
 Symbol("x[2]")
 Symbol("x[3]")

julia> x = variable_names(:x => (3, 2))
6-element Vector{Symbol}:
 Symbol("x[1,1]")
 Symbol("x[2,1]")
 Symbol("x[3,1]")
 Symbol("x[1,2]")
 Symbol("x[2,2]")
 Symbol("x[3,2]")

julia> x == variable_names(:x => [3, 2])
true

julia> vars = variable_names((x = (1, 2), y = 2, z = missing))
5-element Vector{Symbol}:
 Symbol("x[1,1]")
 Symbol("x[1,2]")
 Symbol("y[1]")
 Symbol("y[2]")
 :z

julia> vars == variable_names([:x => (1, 2), :y => 2, :z])
true

julia> variable_names(Dict(:x => (1, 2), :y => 2, :z => missing))
5-element Vector{Symbol}:
 Symbol("y[1]")
 Symbol("y[2]")
 :z
 Symbol("x[1,1]")
 Symbol("x[1,2]")
```
"""
variable_names(a) = Symbol[Symbol(s, suffix) for (s, shape) in _pairs(a) for suffix in _variable_suffix(shape)]::Vector{Symbol}

_variable_suffix(::Missing) = [""]
_variable_suffix(n::Int) = ["[$i]" for i in Base.OneTo(n)]
_variable_suffix(dims::Tuple{Vararg{Int}}) =
    ["[$(join(is, ','))]" for is in Iterators.product(Base.OneTo.(dims)...)]

"""
    _pairs(a)

Return iterator of pairs, if possible, otherwise keep `a`.
"""
_pairs(a::Pair{<:VarName}) = [a]
_pairs(a::Dict{<:VarName}) = pairs(a)
_pairs(a::NamedTuple) = pairs(a)

_pairs(a::Vector) = _pair.(a)
_pair(a::Pair{<:VarName}) = a
_pair(s::VarName) = s => missing

_pairs(a::T) where T = hasmethod(pairs, Tuple{T}) ? pairs(a) : _pair.(a)

"""
    reshape_to_varnames(vec::Vector, varnames)

Turn `vec` into the shape of `varnames`. Reverse flattening from [`variable_names`](@ref).

It holds `reshape_to_varnames(a, variable_names(a)) == a`.

# Examples

```jldoctest
julia> R, vec = polynomial_ring(ZZ, variable_names((x=2, y=(3,2))));

julia> reshape_to_varnames(vec, (x=2, y=(1,2)))
(x = [x[1], x[2]], y = [y[1,1] y[1,2]])
"""
function reshape_to_varnames(vec::Vector, varnames)
    iter = Iterators.Stateful(vec)
    result = reshape_to_varnames(iter, varnames)
    @assert isempty(iter)
    return result
end

reshape_to_varnames(iter::Iterators.Stateful, (s, shape)::Pair) =
    s => _reshape(iter, shape)

reshape_to_varnames(iter::Iterators.Stateful, a::Dict) =
    Dict(s => _reshape(iter, shape) for (s, shape) in a)

reshape_to_varnames(iter::Iterators.Stateful, a::NamedTuple) =
    NamedTuple(s => _reshape(iter, shape) for (s, shape) in pairs(a))

reshape_to_varnames(iter::Iterators.Stateful, a::Vector) =
    [__reshape(iter, s) for s in a]

_reshape(iter, ::Missing) = popfirst!(iter)
_reshape(iter, n::Int) = collect(Iterators.take(iter, n))
_reshape(iter, dims::Tuple{Vararg{Int}}) =
    reshape(collect(Iterators.take(iter, prod(dims))), dims)

__reshape(iter, s::VarName) = popfirst!(iter)
__reshape(iter, (_, a)::Pair) = _reshape(iter, a)

"""
    @with_variable_names f(..., s::Vector{Symbol}, ...)

Add method `f(..., s, ...)` which calls the underlying method after
transforming `s` via a call to `variable_names`.

Format for `@f` is `@f(..., (x[1,2], y[2], z), ...)` or `@f(..., x[3,2], ...)`.
Calling `@f(..., x, ...)` expects `f` to return a single generator.

# Examples

# as if :D
```jldoctest
julia> f(a, s::Union{Vector{Symbol}, Symbol}) = a, s
f (generic function with 1 method)

julia> @with_variable_names f(a, s::Vector{Symbol})
@f (macro with 2 methods)

julia> f
f (generic function with 2 methods)

julia> f("hello", (x = (1, 2), y = 2, z = missing))
("hello", (x = [Symbol("x[1,1]") Symbol("x[1,2]")], y = [Symbol("y[1]"), Symbol("y[2]")], z = :z))

julia> @f "hello" (x = (1, 2), y = 2, z = missing)
"hello"

julia> x
1×2 Matrix{Symbol}:
 Symbol("x[1,1]")
 Symbol("x[1,2]")

julia> y
2-element Vector{Symbol}:
 Symbol("y[1]")
 Symbol("y[2]")

julia> z
 :z
```
"""

macro with_variable_names(e::Expr)
    @assert Base.isexpr(e, :call)
    f = esc(e.args[1])
    args0, _, args1 = _splice(_is_symbol_vector, e.args[2:end])
    s = gensym("varnames") # `:s` could appear in `args0`
    return quote
        function $f($(args0...), $s, $(args1...))
            X, gens = $f($(args0...), variable_names($s), $(args1...))
            return X, reshape_to_varnames(gens, $s)
        end

        macro $f($(args0...), $s::Symbol, $(args1...))
            quote
                X, gen = $$f($$(args0...), $(QuoteNode($s)), $$(args1...))
                $(esc($s)) = gen
                X
            end
        end

        macro $f($(args0...), $s::Expr, $(args1...))
            xs = _expr_pairs($s)
            varnames = Expr(:tuple, (:($k = $v) for (k, v) in xs)...)
            return Expr(:block,
                :((X, gens) = $$f($$(args0...), $varnames, $$(args1...))),
                (:($(esc(x)) = gens.$x) for (x, _) in xs)...,
                :(X)
            )
        end
    end
end

function _expr_pairs(e::Expr)
    Base.isexpr(e, :ref) && (e = Expr(:tuple, e))
    @assert Base.isexpr(e, :tuple)
    return _expr_pair.(e.args)
end
_expr_pair(e::Expr) = (@assert Base.isexpr(e, :ref); e.args[1] => Tuple{Vararg{Int}}(e.args[2:end]))
_expr_pair(s::Symbol) = s => missing


_is_symbol_vector(a) = Base.isexpr(a, :(::), 2) && a.args[2] == :(Vector{Symbol})

function _splice(f::Function, v::Vector)
    k = findfirst(f, v)
    @assert k !== nothing "@with_variable_names called on function without `Vector{Symbol}`"
    return v[1:k-1], v[k], v[k+1:end]
end
