"""
    @alias alias_name real_name

This macro defines a const variable `alias_name` set to `real_name`, and
adds a docstring to it that explains the relation.

Moreover, `alias_name` is exported if `real_name` is exported from
the current module, and vice-versa. Note that for this to work, the `export`
statement must come before the invocation of this macro.

This can be used to provide aliases like `is_isomorphic` for `is_isomorphic`
(or vice versa); unlike `@deprecate`, methods can be installed under either name,
as the two names really refer to the same object.

The alias also gets a special docstring that points out that it is an alias

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> foo(x::Int) = x

julia> @alias bar foo

julia> foo === bar
true
"""
macro alias(alias_name::Symbol, real_name::Symbol)
    result = quote
        if isdefined($__module__, $(QuoteNode(alias_name)))
            $alias_name === $real_name || error("Alias '$($(string(alias_name)))' is already defined with a different value")
        else
            @doc """
                $($(string(alias_name)))

            Alias for `$($(string(real_name)))`.
            """
            const $alias_name = $real_name
        end
        if $(QuoteNode(real_name)) in names($__module__)
            export $alias_name
        elseif $(QuoteNode(alias_name)) in names($__module__)
            export $real_name
        end
    end
    return esc(result)
end
