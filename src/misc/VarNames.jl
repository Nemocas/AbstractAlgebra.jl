using Base.Iterators

# const VarName = Union{Symbol, AbstractString, Char}
# const VarShape = Union{Missing, Int, Tuple{Vararg{Int}}}
# const VarNames = Union{
#     Pair{<:VarName, <:VarShape},
#     Dict{<:VarName, <:VarShape},
#     NamedTuple{<:Any, Tuple{Vararg{VarShape}}},
#     VarNamesVector
#     Vector{Union{<:VarName, Pair{<:VarName, <:VarShape}}}
# }

"""
    variable_names(a) -> Vector{Symbol}

Create proper variable names from `a`.

#Examples

```jldoctest
julia> variable_names([:x, :y])
2-element Vector{String}:
 "x"
 "y"

julia> variable_names(:x => 3)
3-element Vector{String}:
 "x[1]"
 "x[2]"
 "x[3]"

julia> x = variable_names(:x => (3, 2))
6-element Vector{String}:
 "x[1,1]"
 "x[2,1]"
 "x[3,1]"
 "x[1,2]"
 "x[2,2]"
 "x[3,2]"

julia> x == variable_names(:x => [3, 2])
true

julia> vars = variable_names((x = (1, 2), y = 2, z = missing))
5-element Vector{String}:
 "x[1,1]"
 "x[1,2]"
 "y[1]"
 "y[2]"
 "z"

julia> vars == variable_names([:x => (1, 2), :y => 2, :z])
true

julia> variable_names(Dict(:x => (1, 2), :y => 2, :z => missing))
5-element Vector{String}:
 "y[1]"
 "y[2]"
 "z"
 "x[1,1]"
 "x[1,2]"
```
"""
variable_names(a) = Symbol[Symbol(s, suffix) for (s, shape) in _pairs(a) for suffix in _variable_suffix(shape)]::Vector{Symbol}

_variable_suffix(::Missing) = [""]
_variable_suffix(n::Int) = ["[$i]" for i in Base.OneTo(n)]
_variable_suffix(dims) = ["[$(join(is, ','))]" for is in Iterators.product(Base.OneTo.(dims)...)]

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
julia> R, vec = polynomial_ring(ZZ, variable_names((x=(1,2), y=2, z=missing)))
(Multivariate Polynomial Ring in x[1,1], x[1,2], y[1], y[2], z over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x[1,1], x[1,2], y[1], y[2], z])

julia> x, y, z = reshape_to_varnames(vec, (x=(1,2), y=2, z=missing))
(x = AbstractAlgebra.Generic.MPoly{BigInt}[x[1,1] x[1,2]], y = AbstractAlgebra.Generic.MPoly{BigInt}[y[1], y[2]], z = z)
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
_reshape(iter, dims) = reshape(collect(Iterators.take(iter, prod(dims))), Tuple(dims))

__reshape(iter, s::VarName) = popfirst!(iter)
__reshape(iter, (_, a)::Pair) = _reshape(iter, a)

"""
    @add_varnames_interface f(args..., varnames)

Add method `X, vars = f(args..., varnames)` and macro `X = @f args... varnames...` to current scope. Read on.

---

    X, gens::Vector{T} = f(args..., vars::Vector{Symbol})

This underlying method is assumed to exist. The existence is currently not tested.

---

    X, vars = f(args..., varnames)

Compute `X` and `gens` by the underlying `f` method. Then reshape `gens` into the shape defined by `varnames`.

---

    X = @f args... varnames...
    X = @f args... (varnames...)

Compute `X` and `gens` by the underlying `f` method. Then introduce the `varnames` into the current scope.

# Examples

```jldoctest
julia> f(a, s::Vector{Symbol}) = a, String.(s)
f (generic function with 1 method)

julia> @add_varnames_interface f(a, s::Vector{Symbol})
@f (macro with 1 method)

julia> f
f (generic function with 2 methods)

julia> f("hello", (x = (1, 2), y = 2, z = missing))
("hello", (x = ["x[1,1]" "x[1,2]"], y = ["y[1]", "y[2]"], z = "z"))

julia> @f "hello" x[1, 2] y[2] z
"hello"

julia> x
1×2 Matrix{String}:
 "x[1,1]"  "x[1,2]"

julia> y
2-element Vector{String}:
 "y[1]"
 "y[2]"

julia> z
 "z"

julia> v = [1,4,1];

julia> @f "variable dims show case" x[v...]; x
1×4×1 Array{String, 3}:
[:, :, 1] =
 "y[1,1,1]"  "y[1,2,1]"  "y[1,3,1]"  "y[1,4,1]"
```
"""
macro add_varnames_interface(e::Expr)
    @assert Base.isexpr(e, :call)
    f = esc(e.args[1])
    args = e.args[2:end-1]
    s = gensym("varnames")
    kv = gensym("kv")
    # TODO cope with `where` clause and keyvalues
    return quote
        # TODO assert existence of underlying method
        function $f($(args...), $s; $kv...)
            X, gens = $f($(args...), variable_names($s); $kv...)
            return X, reshape_to_varnames(gens, $s)
        end

        macro $f($(args...), $s...)
            xs = _expr_pairs($s)
            varnames = Expr(:tuple, (:($k = $v) for (k, v) in xs)...)
            return Expr(:block,
                :((X, gens) = $$f($$(args...), $varnames)),
                (:($(esc(x)) = gens.$x) for (x, _) in xs)...,
                :(X)
            )
        end
    end
end

_expr_pairs(a::Tuple{Vararg{Union{Expr, Symbol}}}) = _expr_pair.(a)
# for `@f args... (varnames...)` variant:
_expr_pairs((a,)::Tuple{Expr}) = Base.isexpr(a, :tuple) ? _expr_pair.(a.args) : (_expr_pair(a),)
_expr_pair(e::Expr) = (@assert Base.isexpr(e, :ref) "$(e.head) ≠ :ref"; e.args[1] => Expr(:tuple, :($(esc.(e.args[2:end])...))))
_expr_pair(s::Symbol) = s => missing

"""
    @add_varname_interface f(args0..., varname::VarName, args1...)

Add method `X, vars = f(args0..., varname::VarName, args1...)` and macro `X = @f args0... varname::Symbol, args1...` to current scope. Read on.

---

    X, gen::T = f(args0..., var::Symbol, args1...)

This underlying method is assumed to exist. The existence is currently not tested.

---

    X, gen = f(args0..., varname::VarName, args1...) = f(args0..., Symbol(varname), args1...)

---

    X = @f args0... varname::Symbol args1...

Compute `X` and `gen` by the underlying `f` method. Then introduce `varname = gen` into the current scope.

# Examples

```jldoctest
julia> f(a, s::Symbol) = a, s
f (generic function with 1 method)

julia> @add_varname_interface f(a, s::Vector{Symbol})
@f (macro with 2 methods)

julia> f
f (generic function with 2 methods)

julia> f("hello", "x")
("hello", :x)

julia> @f "hello" x
"hello"

julia> x
:x
```
"""
macro add_varname_interface(e::Expr)
    @assert Base.isexpr(e, :call)
    f = esc(e.args[1])
    args = e.args[2:end-1]
    s = gensym("varnames")
    kv = gensym("kv")
    return quote
        $f($(args...), $s::Union{AbstractString, Char}; $kv...) =
            $f($(args...), Symbol($s); $kv...)

        macro $f($(args...), $s::Symbol)
            quote
                X, gen = $$f($$(args...), $(QuoteNode($s)))
                $(esc($s)) = gen
                X
            end
        end
    end
end

@add_varname_interface SparsePolynomialRing(R, s)
@add_varnames_interface LaurentPolynomialRing(R, s)
@add_varnames_interface power_series_ring(R, prec, s)
@add_varnames_interface free_associative_algebra(R, s)
@add_varnames_interface polynomial_ring(R, s)
@add_varname_interface polynomial_ring(R, s)
@add_varname_interface laurent_series_ring(R, prec, s)
@add_varname_interface laurent_series_field(R, prec, s)
@add_varname_interface number_field(p, s) # TODO what about strange `t` parameter
@add_varname_interface FunctionField(p, s)
