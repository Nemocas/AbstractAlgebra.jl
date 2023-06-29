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

raw"""
    variable_names(a...) -> Vector{Symbol}
    variable_names(a::Tuple) -> Vector{Symbol}

Create proper variable names from `a`.

#Examples

```jldoctest
julia> variable_names(:x, :y)
[:x, :y]

julia> variable_names(:x => (1, 2), :y => 2, :z)
[Symbol("x[1,1]"), Symbol("x[1,2]"), Symbol("y[1]"), Symbol("y[2]"), :z]

julia> variable_names(["x$i$j" for i in 0:2, j in 0:1], 'y')
[:x00, :x10, :x20, :x01, :x11, :x21, :y]

```
"""
variable_names(as...) = variable_names(as)
variable_names(as::Tuple) = [x for a in as for x in _variable_names(a)]::Vector{Symbol}

_variable_names(a::AbstractArray{<:VarName}) = Symbol.(a)
_variable_names(s::VarName) = [Symbol(s)]
_variable_names((s, n)::Pair{<:VarName, Int}) = Symbol.(s, '[', Base.OneTo(n), ']')
_variable_names((s, dims)::Pair{<:VarName, <:Tuple{Vararg{Int}}}) = Symbol.(s, '[', join.(Tuple.(CartesianIndices(dims)), ','), ']')
_variable_names((s, _)::Pair{<:VarName, Missing}) = [Symbol(s)]

raw"""
    reshape_to_varnames(vec::Vector{T}, varnames...) :: Tuple{Array{<:Any, T}}
    reshape_to_varnames(vec::Vector{T}, varnames::Tuple) :: Tuple{Array{<:Any, T}}

Turn `vec` into the shape of `varnames`. Reverse flattening from [`variable_names`](@ref).

# Examples

```jldoctest
julia> s = ([:a, :b], :x => (1, 2), :y => 2, :z);

julia> reshape_to_varnames(variable_names(s...), s...)
([:a, :b], [Symbol("x[1,1]"), Symbol("x[1,2]")], [Symbol("y[1]"), Symbol("y[2]")], :z)

julia> R, vec = polynomial_ring(ZZ, variable_names(s...))
(Multivariate Polynomial Ring in a, b, x[1,1], x[1,2], y[1], y[2], z over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[a, b, x[1,1], x[1,2], y[1], y[2], z])

julia> (a, b), x, y, z = reshape_to_varnames(vec, s...)
(AbstractAlgebra.Generic.MPoly{BigInt}[a, b], AbstractAlgebra.Generic.MPoly{BigInt}[x[1,1] x[1,2]], AbstractAlgebra.Generic.MPoly{BigInt}[y[1], y[2]], z)

julia> R, (a, b), x, y, z = polynomial_ring(ZZ, s...)
(Multivariate Polynomial Ring in a, b, x[1,1], x[1,2], y[1], y[2], z over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[a, b], AbstractAlgebra.Generic.MPoly{BigInt}[x[1,1] x[1,2]], AbstractAlgebra.Generic.MPoly{BigInt}[y[1], y[2]], z)

```
"""
reshape_to_varnames(vec::Vector, varnames...) = reshape_to_varnames(vec, varnames)
function reshape_to_varnames(vec::Vector, varnames::Tuple)
    iter = Iterators.Stateful(vec)
    result = Tuple(_reshape_to_varnames(iter, x) for x in varnames)
    @assert isempty(iter)
    return result
end

_reshape_to_varnames(iter::Iterators.Stateful, ::VarName) = popfirst!(iter)
_reshape_to_varnames(iter::Iterators.Stateful, (_, shape)::Pair) = _reshape(iter, shape)
_reshape_to_varnames(iter::Iterators.Stateful, a::AbstractArray{<:VarName}) = _reshape(iter, size(a))

_reshape(iter, n::Int) = collect(Iterators.take(iter, n))
_reshape(iter, dims) = reshape(collect(Iterators.take(iter, prod(dims))), Tuple(dims))
_reshape(iter, ::Missing) = popfirst!(iter)

raw"""
    @add_varnames_interface f(args..., varnames)

Add method `X, vars = f(args..., varnames)` and macro `X = @f args... varnames...` to current scope. Read on.

---

    X, gens::Vector{T} = f(args..., vars::Vector{Symbol})

This underlying method is assumed to exist. The existence is currently not tested.

---

    X, vars... = f(args..., varnames...)

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

julia> f("hello", :x => (1, 2), :y => 2, :z)
("hello", ["x[1,1]" "x[1,2]"], ["y[1]", "y[2]"], "z")

julia> f("projective", ["x$i$j" for i in 0:1, j in 0:1], [:y0, :y1], :z)
("projective", ["x00" "x01"; "x10" "x11"], ["y0", "y1"], "z")

julia> f("fun inputs", 'a':'g', Symbol.('x':'z', [0 1])
("fun inputs", ["a", "b", "c", "d", "e", "f", "g"], ["x0" "x1"; "y0" "y1"; "z0" "z1"])

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

julia> v = [1,4,1]; @f "variable dims" x[v...]; x
1×4×1 Array{String, 3}:
[:, :, 1] =
 "y[1,1,1]"  "y[1,2,1]"  "y[1,3,1]"  "y[1,4,1]"

```

# Caveats

For the macro variant, all to be introduced names have to be given explicitly.
`@add_varnames_interface` can only be called with unqualified names.
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
        function $f($(args...), $s...; $kv...)
            X, gens = $f($(args...), variable_names($s); $kv...)
            return X, reshape_to_varnames(gens, $s)...
        end

        macro $f($(args...), $s...)
            xs = _expr_pairs($s)
            gens = esc.(first.(xs))
            varnames = (:($(QuoteNode(k)) => $v) for (k, v) in xs)
            return quote
                X, $(gens...) = $$f($(esc.($(args...))), $(varnames...))
                X
            end
        end
    end
end

_expr_pairs(a::Tuple{Vararg{Union{Expr, Symbol}}}) = _expr_pair.(a)
_expr_pairs((a,)::Tuple{Expr}) = Base.isexpr(a, :tuple) ? _expr_pair.(a.args) : (_expr_pair(a),) # for `@f args... (varnames...)` variant
_expr_pair(e::Expr) = (@assert Base.isexpr(e, :ref) "$(e.head) ≠ :ref"; e.args[1] => Expr(:tuple, esc.(e.args[2:end])...))
_expr_pair(s::Symbol) = s => missing

raw"""
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
    s = gensym("varname")
    kv = gensym("kv")
    return quote
        $f($(args...), $s::Union{AbstractString, Char}; $kv...) =
            $f($(args...), Symbol($s); $kv...)

        macro $f($(args...), $s::Symbol)
            quote
                X, gen = $$f($(esc.($(args...))), $(QuoteNode($s)))
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
