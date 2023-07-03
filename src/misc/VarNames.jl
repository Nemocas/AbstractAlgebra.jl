import MacroTools as MT

@assert VarName === Union{Symbol, AbstractString, Char}
const VarShape = Union{Missing, Int, Tuple{Vararg{Int}}}
const VarNames = Union{
    VarName,
    AbstractArray{<:VarName},
    Pair{<:VarName, <:VarShape},
}

req(cond, msg) = cond || throw(ArgumentError(msg))

@doc raw"""
    variable_names(a...) -> Vector{Symbol}
    variable_names(a::Tuple) -> Vector{Symbol}

Create proper variable names from `a`.

#Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> AbstractAlgebra.variable_names(:x, :y)
2-element Vector{Symbol}:
 :x
 :y

julia> AbstractAlgebra.variable_names(:x => (1, 2), :y => 2, :z)
5-element Vector{Symbol}:
 Symbol("x[1,1]")
 Symbol("x[1,2]")
 Symbol("y[1]")
 Symbol("y[2]")
 :z

julia> AbstractAlgebra.variable_names(["x$i$j" for i in 0:2, j in 0:1], 'y')
7-element Vector{Symbol}:
 :x00
 :x10
 :x20
 :x01
 :x11
 :x21
 :y

```
"""
variable_names(as::VarNames...) = variable_names(as)
variable_names(as::Tuple{Vararg{VarNames}}) = Symbol[x for a in as for x in _variable_names(a)]

_variable_names(a::AbstractArray{<:VarName}) = Symbol.(a)
_variable_names(s::VarName) = [Symbol(s)]
_variable_names((s, _)::Pair{<:VarName, Missing}) = [Symbol(s)]
_variable_names((s, n)::Pair{<:VarName, Int}) = Symbol.(s, '[', Base.OneTo(n), ']')
_variable_names((s, dims)::Pair{<:VarName, <:Tuple{Vararg{Int}}}) = Symbol.(s, '[', join.(Tuple.(CartesianIndices(dims)), ','), ']')

@doc raw"""
    reshape_to_varnames(vec::Vector{T}, varnames...) :: Tuple{Array{<:Any, T}}
    reshape_to_varnames(vec::Vector{T}, varnames::Tuple) :: Tuple{Array{<:Any, T}}

Turn `vec` into the shape of `varnames`. Reverse flattening from [`variable_names`](@ref).

# Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> s = ([:a, :b], :x => (1, 2), :y => 2, :z);

julia> AbstractAlgebra.reshape_to_varnames(AbstractAlgebra.variable_names(s...), s...)
([:a, :b], [Symbol("x[1,1]") Symbol("x[1,2]")], [Symbol("y[1]"), Symbol("y[2]")], :z)

julia> R, vec = polynomial_ring(ZZ, AbstractAlgebra.variable_names(s...))
(Multivariate polynomial ring in 7 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[a, b, x[1,1], x[1,2], y[1], y[2], z])

julia> (a, b), x, y, z = AbstractAlgebra.reshape_to_varnames(vec, s...)
(AbstractAlgebra.Generic.MPoly{BigInt}[a, b], AbstractAlgebra.Generic.MPoly{BigInt}[x[1,1] x[1,2]], AbstractAlgebra.Generic.MPoly{BigInt}[y[1], y[2]], z)

julia> R, (a, b), x, y, z = polynomial_ring(ZZ, s...)
(Multivariate polynomial ring in 7 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[a, b], AbstractAlgebra.Generic.MPoly{BigInt}[x[1,1] x[1,2]], AbstractAlgebra.Generic.MPoly{BigInt}[y[1], y[2]], z)

```
"""
reshape_to_varnames(vec::Vector, varnames::VarNames...) = reshape_to_varnames(vec, varnames)
function reshape_to_varnames(vec::Vector, varnames::Tuple{Vararg{VarNames}})
    iter = Iterators.Stateful(vec)
    result = Tuple(_reshape_to_varnames(iter, x) for x in varnames)
    @assert isempty(iter)
    return result
end

_reshape_to_varnames(iter::Iterators.Stateful, ::VarName) = popfirst!(iter)
_reshape_to_varnames(iter::Iterators.Stateful, (_, shape)::Pair{<:VarName, <:VarShape}) = _reshape(iter, shape)
_reshape_to_varnames(iter::Iterators.Stateful, a::AbstractArray{<:VarName}) = _reshape(iter, size(a))

_reshape(iter, ::Missing) = popfirst!(iter)
_reshape(iter, n::Int) = collect(Iterators.take(iter, n))
_reshape(iter, dims::Tuple{Vararg{Int}}) = reshape(collect(Iterators.take(iter, prod(dims))), dims)

function _varname_interface(e::Expr, @nospecialize s::Union{Expr, Symbol})
    ex = MT.isexpr(e, :(=), :function) ? e : Expr(:(=), e, :())
    d = MT.splitdef(ex)

    callf = esc(d[:name])
    f = esc(MT.postwalk(x -> MT.@capture(x, a_.b_) ? b : x, d[:name]))
    wheres = esc.(d[:whereparams])

    args = d[:args][begin:end-1]
    splitargs = MacroTools.splitarg.(args)
    args = esc.(args)
    req(all(((_, _, slurp, default),) -> (slurp, default) === (false, nothing), splitargs),
        "Default and slurp arguments currently not supported")
    req(isempty(d[:kwargs]), "Keyword arguments currently not supported")
    argnames = first.(splitargs)
    req(all(!isnothing, argnames), "Nameless arguments currently not supported")
    argnames = esc.(argnames)

    s = esc(s)
    argtypes = esc.(a[2] for a in splitargs)
    argtypes = :(Tuple{$(argtypes...), $s} where {$(wheres...)})
    base = f == callf ?
        :(req(hasmethod($f, $argtypes), "base method of `$($f)` for $($argtypes) missing")) :
        :($f($(args...), s::$s; kv...) where {$(wheres...)} = $callf($(argnames...), s; kv...))

    return f, args, argnames, wheres, base
end

@doc raw"""
    @varnames_interface [M.]f(args..., varnames)

Add method `X, vars = f(args..., varnames...)` and macro `X = @f args... varnames...` to current scope. Read on.

---

    X, gens::Vector{T} = f(args..., varnames::Vector{Symbol})

Base method. If `M` is given, this calls `M.f`, otherwise, it has to exist already.

---

    X, vars... = f(args..., varnames...)

Compute `X` and `gens` via the base method. Then reshape `gens` into the shape defined by `varnames`.
Each `varnames` argument can be either an Array of `VarName`s, or `s::VarName => dims`. The latter means an `Array` of size `dims`.

---

    X, x::Vector{T} = f(args..., n::Int, s::VarName = :x)

Shorthand for `X, x = f(args..., ["$s$i" for i in 1:n])`.

---

    X = @f args... varnames...
    X = @f args... (varnames...)

As `f(args..., varnames...)`, and also introduce the `varnames` into the current scope.

# Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> f(a, s::Vector{Symbol}) = a, String.(s)
f (generic function with 1 method)

julia> AbstractAlgebra.@varnames_interface f(a, s)
@f (macro with 1 method)

julia> f
f (generic function with 5 methods)

julia> f("hello", :x, :y, :z)
("hello", "x", "y", "z")

julia> f("hello", :x => (1, 2), :y => 2, :z)
("hello", ["x[1,1]" "x[1,2]"], ["y[1]", "y[2]"], "z")

julia> f("projective", ["x$i$j" for i in 0:1, j in 0:1], [:y0, :y1], :z)
("projective", ["x00" "x01"; "x10" "x11"], ["y0", "y1"], "z")

julia> f("fun inputs", 'a':'g', Symbol.('x':'z', [0 1]))
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
 "x[1,1,1]"  "x[1,2,1]"  "x[1,3,1]"  "x[1,4,1]"

```
"""
macro varnames_interface(e::Expr, options...)
    f, args, argnames, wheres, base = _varname_interface(e, :(Vector{Symbol}))
    fancy_method = quote
        $f($(args...), s::VarNames...; kv...) where {$(wheres...)} = $f($(argnames...), s; kv...)
        function $f($(args...), s::Tuple{Vararg{VarNames}}; kv...) where {$(wheres...)}
            X, gens = $f($(argnames...), variable_names(s...); kv...)
            return X, reshape_to_varnames(gens, s...)...
        end
    end

    opts = parse_options(options, Dict(:n => :n, :macros => :(:all)), Dict(:macros => QuoteNode.([:no, :tuple, :all])))
    en = n = opts[:n]
    if n isa Symbol
        en = :(Base.OneTo($n))
    elseif n isa Expr
        req(n.head === :call, "Value to option `n` can be `n`, `n+1`, or similar, not `$n`")
        n = only(x -> x isa Symbol, n.args[2:end])
    end

    n isa Symbol || return :($base; $fancy_method)
    fancy_n_method = :($f($(args...), $n::Int, s::VarName=:x; kv...) where {$(wheres...)} = $f($(argnames...), Symbol.(s, $en); kv...))

    opts[:macros] === :(:no) && return :($base; $fancy_method; $fancy_n_method)
    ss, xs = opts[:macros] === :(:all) ? (:(s::Union{Expr, Symbol}...), :(_expr_pairs(s))) :
        (:(s::Expr), :((req(s.head === :tuple, "the final macro argument must be a tuple"); _expr_pair.(s.args))))
    fancy_macro = quote
        macro $f($(argnames...), $ss)
            xs = $xs
            gens = esc.(first.(xs))
            varnames = (:($(QuoteNode(k)) => $v) for (k, v) in xs)
            return quote
                X, $(gens...) = $$f($$(argnames...), $(varnames...))
                X
            end
        end
    end

    return :($base; $fancy_method; $fancy_n_method; $fancy_macro)
end

function parse_options(kvs::Tuple{Vararg{Expr}}, default::Dict{Symbol}, valid::Dict{Symbol, <:Vector} = Dict{Symbol, Vector{Any}}())
    result = Dict{Symbol, Any}(default)
    for o in kvs
        MT.@capture(o, k_ = v_) || error("only key value options allowed")
        req(k in keys(result), "invalid key value option key `$k`")
        k in keys(valid) && req(v in valid[k], "invalid option `$v` to key `$k`")
        result[k] = v
    end
    return result
end

_expr_pairs(es::Tuple{Vararg{Union{Expr, Symbol}}}) = _expr_pair.(es)
_expr_pairs((e,)::Tuple{Expr}) = MT.@capture(e, (es__,)) ? _expr_pair.(es) : (_expr_pair(e),) # for `@f args... (varnames...,)` variant
_expr_pair(e::Expr) = MT.@capture(e, x_[a__]) ? x => :($(esc.(a)...),) : error("variable name must be like `x` or `x[...]`, not `$e`")
_expr_pair(s::Symbol) = s => missing

@doc raw"""
    @varname_interface [M.]f(args..., varname)

Add method `X, vars = f(args..., varname::VarName)` and macro `X = @f args... varname::Symbol` to current scope. Read on.

---

    X, gen::T = f(args..., varname::Symbol)

Base method. If `M` is given, this calls `M.f`, otherwise, it has to exist already.

---

    f(args..., varname::VarName) = f(args..., Symbol(varname))

---

    X = @f args... varname::Symbol

As `f(args..., varname)`, and also introduce `varname` into the current scope.

# Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> f(a, s::Symbol) = a, s
f (generic function with 1 method)

julia> AbstractAlgebra.@varname_interface f(a, s)
@f (macro with 1 method)

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
macro varname_interface(e::Expr)
    f, args, argnames, wheres, base = _varname_interface(e, :Symbol)
    fancy_method = :($f($(args...), s::Union{AbstractString, Char}; kv...) where {$(wheres...)} = $f($(argnames...), Symbol(s); kv...))
    fancy_macro = :(
        macro $f($(argnames...), s::Symbol)
            quote
                X, $(esc(s)) = $$f($$(argnames...), $(QuoteNode(s)))
                X
            end
        end
        )
    return :($base; $fancy_method; $fancy_macro)
end


@varname_interface Generic.SparsePolynomialRing(R::Ring, s)
@varname_interface Generic.laurent_series_ring(R::Ring, prec::Int, s)
@varname_interface Generic.laurent_series_field(R::Field, prec::Int, s)
@varname_interface Generic.number_field(p::PolyRingElem, s)
@varname_interface Generic.FunctionField(p::PolyRingElem, s)

@varnames_interface Generic.free_associative_algebra(R::Ring, s)
@varnames_interface Generic.LaurentPolynomialRing(R::Ring, s)

# The various optional arguments would result in ambiguities
@varname_interface Generic.power_series_ring(R::Ring, prec::Int, s)
@varnames_interface Generic.power_series_ring(R::Ring, prec::Int, s) macros=:tuple
@varnames_interface Generic.power_series_ring(R::Ring, weights::Vector{Int}, prec::Int, s)
@varnames_interface Generic.power_series_ring(R::Ring, prec::Vector{Int}, s) n=0 macros=:no

@varname_interface polynomial_ring(R::NCRing, s)
@varnames_interface polynomial_ring(R::Ring, s)
# With `Ring <: NCRing`, we need to resolve ambiguities of `polynomial_ring(::Ring, s...)`
polynomial_ring(R::Ring, s::Symbol; kv...) = invoke(polynomial_ring, Tuple{NCRing, Symbol}, R, s; kv...)
polynomial_ring(R::Ring, s::Union{AbstractString, Char}; kv...) = polynomial_ring(R, Symbol(s); kv...)

# TODO: weights in `graded_polynomial_ring` and `power_series_ring`
