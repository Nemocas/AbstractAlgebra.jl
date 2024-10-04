import MacroTools as MT

@assert Symbol <: VarName
const VarNames = Union{
    AbstractArray{<:VarName},
    Pair{<:VarName},
}

req(cond, msg) = cond || throw(ArgumentError(msg))

@doc raw"""
    variable_names(a...) -> Vector{Symbol}
    variable_names(a::Tuple) -> Vector{Symbol}

Create a vector of variable names from a variable name specification.

Each argument can be either an Array of `VarName`s,
or of the form `s::VarName => iter`, or of the form `s::VarName => (iter...)`.
Here `iter` is supposed to be any iterable, typically a range like `1:5`.
The `:s => iter` specification is shorthand for `["s[$i]" for i in iter]`.
Similarly `:s => (iter1, iter2)` is shorthand for `["s[$i,$j]" for i in iter1, j in iter2]`,
and likewise for three and more iterables.

As an alternative `"s#" => iter` is shorthand for `["s$i" for i in iter]`.
This also works for multiple iterators in that`"s#" => (iter1, iter2)`
is shorthand for `["s$i$j" for i in iter1, j in iter2]`.

# Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> AbstractAlgebra.variable_names([:x, :y])
2-element Vector{Symbol}:
 :x
 :y

julia> AbstractAlgebra.variable_names(:x => (0:0, 0:1), :y => 0:1, [:z])
5-element Vector{Symbol}:
 Symbol("x[0, 0]")
 Symbol("x[0, 1]")
 Symbol("y[0]")
 Symbol("y[1]")
 :z

julia> AbstractAlgebra.variable_names("x#" => (0:0, 0:1), "y#" => 0:1)
4-element Vector{Symbol}:
 :x00
 :x01
 :y0
 :y1

julia> AbstractAlgebra.variable_names("x#" => 9:11)
3-element Vector{Symbol}:
 :x9
 :x10
 :x11

julia> AbstractAlgebra.variable_names(["x$i$i" for i in 1:3])
3-element Vector{Symbol}:
 :x11
 :x22
 :x33

julia> AbstractAlgebra.variable_names('a':'c', ['z'])
4-element Vector{Symbol}:
 :a
 :b
 :c
 :z
```
"""
variable_names(as::VarNames...) = variable_names(as)
# brackets = Val(false) effectively replaces `:x` and `"x"` by `"x#"`, used by macro
variable_names(as::Tuple{Vararg{VarNames}}, brackets::Val = Val(true)) =
    Symbol[x for a in as for x in _variable_names(a, brackets)]

_variable_names(a::AbstractArray{<:VarName}, ::Any) = Symbol.(a)

_variable_names((s, axe)::Pair{<:Union{Char, Symbol}}, ::Val{true}) = Symbol.(s, '[', axe, ']')
_variable_names((s, axe)::Pair{<:Union{Char, Symbol}}, ::Val{false}) = check_names(Symbol.(s, axe))

_variable_names((s, axes)::Pair{<:Union{Char, Symbol}, <:Tuple}, ::Val{true}) = Symbol.(s, '[', join.(Iterators.product(axes...), ", "), ']')
_variable_names((s, axes)::Pair{<:Union{Char, Symbol}, <:Tuple}, ::Val{false}) = check_names(Symbol.(s, join.(Iterators.product(axes...))))

_variable_names((s, axe)::Pair{<:AbstractString}, val::Val) = _variable_names(s => (axe,), val)
function _variable_names((s, axes)::Pair{<:AbstractString, <:Tuple}, val::Val)
    c = count("#", s)
    req(c <= 1, """Only a single '#' allowed, but "$s" contains $c of them.
        Please communicate your use case to the Oscar community.""")
    return c == 0 ? _variable_names(Symbol(s) => axes, val) :
        _check_names([Symbol(replace(s, '#' => join(i))) for i in Iterators.product(axes...)], val)
end

"""
    check_names(names) -> names

Warn, if any of the `names` is no valid Julia identifier. Return `names`.
"""
function check_names(names)
    badname = _first(!Meta.isidentifier, names)
    if badname !== nothing
        badname = first(x for x in names if !Meta.isidentifier(x))
        badstring = repr(string(badname))
        @warn "The variable name $badstring sadly is no Julia identifier. " *
            "You can still access it as `var$badstring`."
    end
    return names
end
_check_names(x, ::Val{true}) = x
_check_names(x, ::Val{false}) = check_names(x)

function _first(f, a)
    i = iterate(Iterators.filter(f, a))
    return i === nothing ? nothing : first(i)
end

@doc raw"""
    reshape_to_varnames(vec::Vector{T}, varnames...) :: Tuple{Array{<:Any, T}}
    reshape_to_varnames(vec::Vector{T}, varnames::Tuple) :: Tuple{Array{<:Any, T}}

Turn `vec` into the shape of `varnames`. Reverse flattening from [`variable_names`](@ref).

# Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> s = ([:a, :b], "x#" => (1:1, 1:2), "y#" => 1:2, [:z]);

julia> AbstractAlgebra.reshape_to_varnames(AbstractAlgebra.variable_names(s...), s...)
([:a, :b], [:x11 :x12], [:y1, :y2], [:z])

julia> R, v = polynomial_ring(ZZ, AbstractAlgebra.variable_names(s...))
(Multivariate polynomial ring in 7 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[a, b, x11, x12, y1, y2, z])

julia> (a, b), x, y, z = AbstractAlgebra.reshape_to_varnames(v, s...)
(AbstractAlgebra.Generic.MPoly{BigInt}[a, b], AbstractAlgebra.Generic.MPoly{BigInt}[x11 x12], AbstractAlgebra.Generic.MPoly{BigInt}[y1, y2], AbstractAlgebra.Generic.MPoly{BigInt}[z])

julia> R, (a, b), x, y, z = polynomial_ring(ZZ, s...)
(Multivariate polynomial ring in 7 variables over integers, AbstractAlgebra.Generic.MPoly{BigInt}[a, b], AbstractAlgebra.Generic.MPoly{BigInt}[x11 x12], AbstractAlgebra.Generic.MPoly{BigInt}[y1, y2], AbstractAlgebra.Generic.MPoly{BigInt}[z])
```
"""
reshape_to_varnames(vec::Vector, varnames::VarNames...) =
    reshape_to_varnames(vec, varnames)
function reshape_to_varnames(vec::Vector, varnames::Tuple{Vararg{VarNames}})
    iter = Iterators.Stateful(vec)
    result = _unpeel_reshape_to_varnames(iter, varnames)
    @assert isempty(iter)
    return result
end

function _unpeel_reshape_to_varnames(iter, x::Tuple)
    if length(x) === 1
        return (_reshape_to_varnames(iter, x[1]), )
    else
        return tuple(_reshape_to_varnames(iter, x[1]), _unpeel_reshape_to_varnames(iter, Base.tail(x))...)
    end
end

_reshape_to_varnames(iter::Iterators.Stateful, a::AbstractArray{<:VarName}) =
    _reshape(iter, size(a))
_reshape_to_varnames(iter::Iterators.Stateful, (_, shape)::Pair{<:VarName}) =
    __reshape(iter, shape)

__reshape(iter, axes::Tuple) = _reshape(iter, ntuple(i -> size(axes[i])[1], length(axes)))
__reshape(iter, axe) = _reshape(iter, size(axe))

_reshape(iter, dims) = reshape(collect(Iterators.take(iter, prod(dims))), Tuple(dims))

"""
    keyword_arguments((kvs::Expr...), default::Dict, [valid::Dict]) :: Dict

Mimic usual keyword arguments for usage in macros.

* `kvs`: tuple of Expr of form :(k = v)
* `default`: dictionary providing the allowed keys and their default values
* `valid`: optional `Dict{Symbol, <:AbstractVector}` constraining the valid values for some keys

Return a copy of `default` with the key value pairs from `kvs` applied.

# Example
```jldoctest; setup = :(using AbstractAlgebra)
julia> AbstractAlgebra.keyword_arguments((:(a=1), :(b=:yes)),
       Dict(:a=>0, :b=>:no, :c=>0),
       Dict(:b => [:(:yes), :(:no)]))
Dict{Symbol, Any} with 3 entries:
  :a => 1
  :b => :(:yes)
  :c => 0
```
"""
function keyword_arguments(kvs::Tuple{Vararg{Expr}}, default::Dict{Symbol},
        valid::Dict{Symbol, <:AbstractVector} = Dict{Symbol, Vector{Any}}()) ::
        Dict{Symbol}
    result = Dict{Symbol, Any}(default)
    for o in kvs
        req(MT.@capture(o, k_ = v_), "Only key value options allowed")
        req(k in keys(result), "Invalid key value option key `$k`")
        k in keys(valid) && req(v in valid[k], "Invalid option `$v` to key `$k`")
        result[k] = v
    end
    return result
end

raw"""
    normalise_keyword_arguments(args::Union{Tuple, Vector}) :: Vector

Turn argument list like `(1, a=2; b=3).args` into `(1; a=2, b=3).args`.

Intended to let a macro call `@m(1, a=2; b=3)` mimic a usual function call.

```jldoctest; setup = :(using AbstractAlgebra)
julia> args = AbstractAlgebra.normalise_keyword_arguments(:(1, a=2; b=3).args); :($(args...),)
:((1; a = 2, b = 3))

julia> macro nka_test(args...) AbstractAlgebra.normalise_keyword_arguments(args) end
@nka_test (macro with 1 method)

julia> args = @nka_test(1, a=2; b=3); :($(args...),)
:((1; a = 2, b = 3))

julia> args = @nka_test(1+2, a=x; b=x^2, kw...); :($(args...),)
:((1 + 2; a = x, b = x^2, kw...))
```
"""
function normalise_keyword_arguments(args)
    # Keyword arguments after `;`
    if Meta.isexpr(first(args), :parameters)
        kv = first(args)
        args = args[2:end]
    else
        kv = Expr(:parameters)
    end
    # Keyword arguments without previous `;`
    append!(kv.args, (Expr(:kw, e.args...) for e in args if Meta.isexpr(e, :(=))))
    # normal arguments
    args = (e for e in args if !Meta.isexpr(e, :(=)))
    return [kv, args...]
end

function _eval(m::Core.Module, e::Expr)
    try
        Base.eval(m, e)
    catch err
        if isa(err, UndefVarError)
            @error "Inconveniently, you may only use literals and variables " *
                "from the global scope of the current module (`$m`) " *
                "when using variable name constructor macros"
        end
        rethrow()
    end
end

# input is :([M.]f(args..., s) where {wheres} [ = ... ])
function _splitdef(e::Expr)
    Meta.isexpr(e, (:(=), :function)) || (e = Expr(:(=), e, :()))
    d = MT.splitdef(e)

    req(isempty(d[:kwargs]), "Keyword arguments currently not supported")

    args = d[:args][begin:end-1] # the last argument is just a placeholder
    splitargs = MT.splitarg.(args)
    req(all(((_, _, slurp, default),) -> (slurp, default) === (false, nothing), splitargs),
        "Default and slurp arguments currently not supported")

    argnames = first.(splitargs)
    req(all(!isnothing, argnames), "Nameless arguments currently not supported")

    base_f = d[:name]
    return Dict{Symbol, Any}(
        :base_f => esc(base_f),
        :f => esc(unqualified_name(base_f)),
        :wheres => esc.(d[:whereparams]),
        :args => esc.(args),
        :argnames => esc.(argnames),
        :argtypes => (esc(a[2]) for a in splitargs),
    )
end

unqualified_name(name::Symbol) = name
unqualified_name(name::QuoteNode) = name.value
function unqualified_name(name::Expr)
    req(Meta.isexpr(name, :., 2), "Expected a binding, but `$name` given")
    unqualified_name(name.args[2])
end

function base_method(d::Dict{Symbol},
        @nospecialize s_type::Union{Symbol, Expr})
    f, base_f, wheres = d[:f], d[:base_f], d[:wheres]
    if f == base_f
        argtypes = :(Tuple{$(d[:argtypes]...), $s_type} where {$(wheres...)})
        :(req(hasmethod($f, $argtypes),
              "base method of `$($f)` for $($argtypes) missing"))
    else
        :($f($(d[:args]...), s::$s_type; kv...) where {$(wheres...)} =
            $base_f($(d[:argnames]...), s; kv...))
    end
end

function varnames_method(d::Dict{Symbol})
    f, args, argnames, wheres = d[:f], d[:args], d[:argnames], d[:wheres]
    quote
        $f($(args...), s1::VarNames, s::VarNames...; kv...) where {$(wheres...)} =
            $f($(argnames...), (s1, s...); kv...)
        function $f($(args...), s::Tuple{Vararg{VarNames}}; kv...) where {$(wheres...)}
            X, gens = $f($(argnames...), variable_names(s...); kv...)
            return X, reshape_to_varnames(gens, s...)...
        end
    end
end

function n_vars_method(d::Dict{Symbol}, n, range)
    f, args, argnames, wheres = d[:f], d[:args], d[:argnames], d[:wheres]
    n === :(:no) && return :()
    req(n isa Symbol, "Value to option `n` must be `:no` or " *
        "an alternative name like `m`, not `$n`")
    quote
        $f($(args...), $n::Int, s::VarName=:x; kv...) where {$(wheres...)} =
            $f($(argnames...), Symbol.(s, $range); kv...)
    end
end

function varnames_macro(f, args_count, opt_in, n, range)
    opt_in === :(:yes) || return :()
    quote
        macro $f(args...)
            args = normalise_keyword_arguments(args)
            req(length(args) >= $args_count+2, "Not enough arguments")
            base_args = args[1:$args_count+1] # includes keyword arguments

            m = VERSION > v"9999" ? __module__ : $(esc(:__module__)) # julia issue #51602
            s = _eval(m, :($$_varnames_macro($(args[$args_count+2:end]...))))
            $(_n_vars_macro_support(n, range))

            varnames_macro_code($f, esc.(base_args), s)
        end
    end
end

_varnames_macro(arg::VarName) = Symbol(arg)
_varnames_macro(args::VarNames...) = variable_names(args, Val(false))
_varnames_macro(n::Int, s::VarName=:x) = (n, Symbol(s)) # defer work to `_n_vars_macro_support` outside `eval`

function _n_vars_macro_support(n, range)
    if n == :(:no)
        return :(s isa Tuple{Int, Symbol} && throw(ArgumentError("`Int` argument variant not supported")))
    else
        return quote
            if s isa Tuple{Int, Symbol}
                $n = s[1]
                s = Symbol.(s[2], $range)
            end
        end
    end
end

function varnames_macro_code(f, args, s::Symbol)
    quote
        X, $(esc(s)) = $f($(args...), $(QuoteNode(s)))
        X
    end
end

function varnames_macro_code(f, args, s::Vector{Symbol})
    quote
        X, ($(esc.(s)...),) = $f($(args...), $s)
        X
    end
end

@doc raw"""
    @varnames_interface [M.]f(args..., varnames) macros=:yes n=n range=1:n

Add methods `X, vars = f(args..., varnames...)` and macro `X = @f(args..., varnames...`) to current scope.

# Created methods

    X, gens::Vector{T} = f(args..., varnames::Vector{Symbol})

Base method, called by everything else defined below. If a module `M` is
specified, this is implemented as a call to `M.f`. Otherwise, a method `f` with
this signature must already exist.

---

    X, gens... = f(args..., varnames...; kv...)
    X, gens... = f(args..., varnames::Tuple; kv...)

Compute `X` and `gens` via the base method. Then reshape `gens` into the shape defined by `varnames` according to [`variable_names`](@ref).

The vararg `varnames...` method needs at least one argument to avoid confusion.
Moreover a single `VarName` argument will be dispatched to use a univariate method of `f` if it exists (e.g. `polynomial_ring(R, :x)`).
If you need those cases, use the `Tuple` method.

Keyword arguments are passed on to the base method.

---

    X, x::Vector{T} = f(args..., n::Int, s::VarName = :x; kv...)

Shorthand for `X, x = f(args..., "$s#" => 1:n; kv...)`.
The name of the argument `n` can be changed via the `n` option. The range `1:n` is given via the `range` option.

Setting `n=:no` disables creation of this method.

---

    X = @f(args..., varnames...; kv...)
    X = @f(args..., varnames::Tuple; kv...)
    X = @f(args..., n::Int, s::VarName = :x; kv...)
    X = @f(args..., varname::VarName; kv...)

These macros behave like their `f(args..., varnames; kv...)` counterparts but also introduce the indexed `varnames` into the current scope.
The first version needs at least one `varnames` argument to avoid confusion.
The last version calls the univariate base method if it exists (e.g. `polynomial_ring(R, varname)`).

Setting `macros=:no` disables macro creation.

!!! warning
    Turning `varnames` into a vector of symbols happens by evaluating `variable_names(varnames)` in the global scope of the current module.
    For interactive usage in the REPL this is fine, but in general you have no access to local variables and should not use any side effects in `varnames`.

# Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> f(a, s::Vector{Symbol}) = a, String.(s)
f (generic function with 1 method)

julia> AbstractAlgebra.@varnames_interface f(a, s)
@f (macro with 1 method)

julia> f
f (generic function with 5 methods)

julia> f("hello", [:x, :y, :z])
("hello", ["x", "y", "z"])

julia> f("numbered", 3)
("numbered", ["x1", "x2", "x3"])

julia> f("hello", :x => (1:1, 1:2), :y => 1:2, [:z])
("hello", ["x[1, 1]" "x[1, 2]"], ["y[1]", "y[2]"], ["z"])

julia> f("projective", ["x$i$j" for i in 0:1, j in 0:1], [:y0, :y1], [:z])
("projective", ["x00" "x01"; "x10" "x11"], ["y0", "y1"], ["z"])

julia> f("fun inputs", 'a':'g', Symbol.('x':'z', [0 1]))
("fun inputs", ["a", "b", "c", "d", "e", "f", "g"], ["x0" "x1"; "y0" "y1"; "z0" "z1"])

julia> @f("hello", "x#" => (1:1, 1:2), "y#" => (1:2), [:z])
"hello"

julia> (x11, x12, y1, y2, z)
("x11", "x12", "y1", "y2", "z")

julia> g(a, s::Vector{Symbol}; kv...) = (a, kv...), String.(s)
g (generic function with 1 method)

julia> AbstractAlgebra.@varnames_interface g(a, s)
@g (macro with 1 method)

julia> @g("parameters", [:x, :y], a=1, b=2; c=3)
("parameters", :c => 3, :a => 1, :b => 2)
```
"""
macro varnames_interface(e::Expr, options::Expr...)
    d = _splitdef(e)

    opts = keyword_arguments(options,
        Dict(:n => :n, :range => :(1:n), :macros => :(:yes)),
        Dict(:macros => QuoteNode.([:no, :yes])))

    quote
        $(base_method(d, :(Vector{Symbol})))
        $(varnames_method(d))
        $(n_vars_method(d, opts[:n], opts[:range]))
        $(varnames_macro(d[:f], length(d[:argnames]), opts[:macros], opts[:n], opts[:range]))
    end
end
