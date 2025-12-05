import MacroTools

using Base: ismutabletype

"""
    @attributes typedef

This is a helper macro that ensures that there is storage for attributes in
the type declared in the expression `typedef`, which must be either a `mutable
struct` definition expression, or the name of a `mutable struct` type.

The latter variant is useful to enable attribute storage for types defined in
other packages. Note that `@attributes` is idempotent: when applied to a type
for which attribute storage is already available, it does nothing.

For singleton types, attribute storage is also supported, and in fact always
enabled. Thus it is not necessary to apply this macro to such a type.

!!! note
    When applied to a struct definition this macro adds a new field to the
    struct. For structs without constructor, this will change the signature of
    the default inner constructor, which requires explicit values for every
    field, including the attribute storage field this macro adds. Usually it
    is thus preferable to add an explicit default constructor, as in the
    example below.

# Examples

Applying the macro to a struct definition results in internal storage
of the attributes:
```jldoctest
julia> @attributes mutable struct MyGroup
           order::Int
           MyGroup(order::Int) = new(order)
       end

julia> G = MyGroup(5)
MyGroup(5, #undef)

julia> set_attribute!(G, :isfinite, :true)

julia> get_attribute(G, :isfinite)
true
```

Applying the macro to a typename results in external storage of the
attributes:
```jldoctest
julia> mutable struct MyOtherGroup
           order::Int
           MyOtherGroup(order::Int) = new(order)
       end

julia> @attributes MyOtherGroup

julia> G = MyOtherGroup(5)
MyOtherGroup(5)

julia> set_attribute!(G, :isfinite, :true)

julia> get_attribute(G, :isfinite)
true
```
"""
macro attributes(expr)
   # The following two lines are borrowed from Base.@kwdef
   expr = macroexpand(__module__, expr) # to expand @static
   if expr isa Expr && expr.head === :struct && expr.args[1]
      # Handle the following usage:
      #    @attributes mutable struct Type ... end

      # add member for storing the attributes
      push!(expr.args[3].args, :(__attrs::Dict{Symbol,Any}))
      return quote
        Base.@__doc__($(esc(expr)))
      end
   elseif expr isa Expr && expr.head === :struct && !expr.args[1] && all(x -> x isa LineNumberNode, expr.args[3].args)
      # Ignore application to singleton types:
      #    @attributes struct Singleton end
      return esc(expr)
   elseif expr isa Symbol || (expr isa Expr && expr.head === :. &&
                              length(expr.args) == 2 && expr.args[2] isa QuoteNode) ||
                              (expr isa Expr && expr.head === :curly &&
                                expr.args[1] isa Symbol || (expr.args[1] isa Expr && expr.args[1].head === :. &&
                                length(expr.args[1].args) == 2 && expr.args[1].args[2] isa QuoteNode))
      # Handle the following usage:
      #    @attributes Type
      #    @attributes Module.Type
      #    @attributes Module.Submodule.Type
      #    etc.
      #    @attributes [Module[.Submodule].]Type{T}
      return esc(quote
         # do nothing if the type already has storage for attributes
         if !AbstractAlgebra.is_attribute_storing_type($expr)
           isstructtype($expr) && AbstractAlgebra.ismutabletype($expr) || error("attributes can only be attached to mutable structs")
           let attr = Base.WeakKeyDict{$expr, Dict{Symbol, Any}}()
             AbstractAlgebra._get_attributes(G::$expr) = Base.get(attr, G, nothing)
             AbstractAlgebra._get_attributes!(G::$expr) = Base.get!(() -> Dict{Symbol, Any}(), attr, G)
             AbstractAlgebra._get_attributes(::Type{$expr}) = attr
             AbstractAlgebra.is_attribute_storing_type(::Type{$expr}) = true
           end
         end
      end)
   end
   error("attributes can only be attached to mutable structs")
end

_is_attribute_storing_type(::Type{T}) where T = Base.issingletontype(T) || (isstructtype(T) && ismutabletype(T) && hasfield(T, :__attrs))

# storage for attributes of singletons
const _singleton_attr_storage = Dict{Type, Dict{Symbol, Any}}()

"""
    _get_attributes(G::T) where T

Return the dictionary storing the attributes for `G` if available, otherwise
return `nothing`. If type `T` does not support attribute storage, an exception
is thrown.

This is used to implement `has_attribute` and `get_attribute`.
Alternate attribute storage solutions need only provide custom methods
for `_get_attributes` and `_get_attributes!`.
"""
function _get_attributes(G::T) where T
   if Base.issingletontype(T)
      return Base.get(_singleton_attr_storage, T, nothing)
   end
   is_attribute_storing_type(T) || error("attributes storage not supported for type $T")
   return isdefined(G, :__attrs) ? G.__attrs : nothing
end

"""
    _get_attributes!(G::T) where T

Return the dictionary storing the attributes for `G`, creating it if necessary.
If type `T` does not support attribute storage, an exception is thrown.

This is used to implement `set_attribute!` and `get_attribute!`.
Alternate attribute storage solutions need only provide custom methods
for `_get_attributes` and `_get_attributes!`.
directly.
"""
function _get_attributes!(G::T) where T
   if Base.issingletontype(T)
      return Base.get!(() -> Dict{Symbol, Any}(), _singleton_attr_storage, T)
   end
   is_attribute_storing_type(T) || error("attributes storage not supported for type $T")
   if !isdefined(G, :__attrs)
      G.__attrs = Dict{Symbol, Any}()
   end
   return G.__attrs
end

"""
    is_attribute_storing(G::Any)

Return a boolean indicating whether `G` has the ability to store attributes.
"""
is_attribute_storing(G::T) where T = is_attribute_storing_type(T)

"""
    is_attribute_storing_type(T::Type)

Return a boolean indicating whether instances of type `T` have
the ability to store attributes.
"""
is_attribute_storing_type(::Type{T}) where T = _is_attribute_storing_type(T)

"""
    has_attribute(G::Any, attr::Symbol)

Return a boolean indicating whether `G` has a value stored for the attribute `attr`.
"""
function has_attribute(G::Any, attr::Symbol)
   D = _get_attributes(G)
   return D isa Dict && haskey(D, attr)
end

"""
    get_attribute(f::Function, G::Any, attr::Symbol)

Return the value stored for the attribute `attr`, or if no value has been set,
return `f()`.

This is intended to be called using `do` block syntax.

```julia
get_attribute(obj, attr) do
    # default value calculated here if needed
    ...
end
```
"""
function get_attribute(f, G::Any, attr::Symbol)
   D = _get_attributes(G)
   D isa Dict && return get(f, D, attr)
   return f()
end

"""
    get_attribute(G::Any, attr::Symbol, default::Any = nothing)

Return the value stored for the attribute `attr`, or if no value has been set,
return `default`.
"""
function get_attribute(G::Any, attr::Symbol, default::Any = nothing)
   D = _get_attributes(G)
   D isa Dict && return get(D, attr, default)
   return default
end

# the following method is necessary to disambiguate between the above
# methods when the default value is a Symbol
function get_attribute(G::Any, attr::Symbol, default::Symbol)
   D = _get_attributes(G)
   D isa Dict && return get(D, attr, default)
   return default
end

"""
    get_attribute!(f::Function, G::Any, attr::Symbol)

Return the value stored for the attribute `attr` of `G`, or if no value has been set,
store `key => f()` and return `f()`.

This is intended to be called using `do` block syntax.

```julia
get_attribute!(obj, attr) do
    # default value calculated here if needed
    ...
end
```
"""
function get_attribute!(f, G::Any, attr::Symbol)
   D = _get_attributes!(G)
   return Base.get!(f, D, attr)
end

"""
    get_attribute!(G::Any, attr::Symbol, default::Any)

Return the value stored for the attribute `attr` of `G`, or if no value has been set,
store `key => default`, and return `default`.
"""
function get_attribute!(G::Any, attr::Symbol, default::Any)
   D = _get_attributes!(G)
   return Base.get!(D, attr, default)
end

# the following method is necessary to disambiguate between the above
# methods when the default value is a Symbol
function get_attribute!(G::Any, attr::Symbol, default::Symbol)
   D = _get_attributes!(G)
   return Base.get!(D, attr, default)
end

"""
    set_attribute!(G::Any, data::Pair{Symbol, <:Any}...)

Attach the given sequence of `key=>value` pairs as attributes of `G`.
"""
function set_attribute!(G::Any, data::Pair{Symbol, <:Any}...)
   D = _get_attributes!(G)
   for d in data
      push!(D, d)
   end
   return nothing
end

"""
    set_attribute!(G::Any, attr::Symbol, value::Any)

Attach the given `value` as attribute `attr` of `G`.
"""
function set_attribute!(G::Any, attr::Symbol, value::Any)
   D = _get_attributes!(G)
   D[attr] = value
   return nothing
end


"""
    @attr RetType funcdef

This macro is applied to the definition of a unary function, and enables
caching ("memoization") of its return values based on the argument. This
assumes the argument supports attribute storing (see [`@attributes`](@ref))
via [`get_attribute!`](@ref).

The name of the function is used as name for the underlying attribute.

The macro works the same for unary functions with keyword arguments,
but ignores the keyword arguments when caching the result, i.e.
different calls with different keyword arguments will return
the identical (cached) result. In case that there is no result cached yet,
the function is called with the given keyword arguments.

Effectively, this turns code like this:
```julia
@attr RetType function myattr(obj::Foo)
   # ... expensive computation
   return result
end
```
into something essentially equivalent to this:
```julia
function myattr(obj::Foo)
  return get_attribute!(obj, :myattr) do
    # ... expensive computation
    return result
  end::RetType
end
```

# Examples
```jldoctest
julia> @attributes mutable struct Foo
           x::Int
           Foo(x::Int) = new(x)
       end;

julia> @attr Int function myattr(obj::Foo)
                println("Performing expensive computation")
                return factorial(obj.x)
             end;

julia> obj = Foo(5);

julia> myattr(obj)
Performing expensive computation
120

julia> myattr(obj) # second time uses the cached result
120

```
"""
macro attr(rettype, expr::Expr)
   d = MacroTools.splitdef(expr)
   length(d[:args]) == 1 || throw(ArgumentError("Only unary functions are supported"))

   # store the original function name
   name = d[:name]

   # take the original function and rename it; use a distinctive name to
   # minimize the risk of accidental conflicts. We deliberately don't
   # use gensym anymore because this interacts badly with Revise.
   compute_name = Symbol("__compute_$(name)__")
   compute_def = copy(d)
   compute_def[:name] = compute_name
   compute = MacroTools.combinedef(compute_def)

   argname = d[:args][1]
   wrapper_def = copy(d)
   wrapper_def[:name] = name
   wrapper_def[:body] = quote
      return get_attribute!(() -> $(compute_name)($argname; $(first.(MacroTools.splitarg.(d[:kwargs]))...)), $(argname), Symbol($(string(name))))::$(rettype)
   end
   # insert the correct line number, so that `functionloc(name)` works correctly
   wrapper_def[:body].args[1] = __source__
   wrapper = MacroTools.combinedef(wrapper_def)

   result = quote
      $(compute)
      Base.@__doc__ $(wrapper)
   end


   # TODO: perhaps also generate a tester, i.e.:
   #    has_attrname(obj::T) = has_attribute(obj, :attrname))
   #   auto-generated documentation for that??!?

    # we must prevent Julia from applying gensym to all locals, as these
    # substitutions do not get applied to the quoted part of the new body,
    # leading to trouble if the wrapped function has arguments (as the
    # argument names will be replaced, but not their uses in the quoted part
    return esc(result)
end
