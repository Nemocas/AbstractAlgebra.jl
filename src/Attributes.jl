export @attributes, has_attribute, get_attribute, get_attribute!, set_attribute!

if VERSION >= v"1.7"
   import Base: ismutabletype
else
   function ismutabletype(@nospecialize(t::Type))
       t = Base.unwrap_unionall(t)
       return isa(t, DataType) && t.mutable
   end
end


"""
    @attributes typedef

This is a helper macro that ensures that there is storage for attributes in
the type declared in the expression `typedef`, which must be either a `mutable
struct` definition expression, or the name of a `mutable struct` type.

The latter variant is useful to enable attribute storage for types defined in
other packages. Note that `@attributes` is idempotent: when applied to a type
for which attribute storage is already available, it does nothing.

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
```jldoctest; setup = :(using AbstractAlgebra)
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
```jldoctest; setup = :(using AbstractAlgebra)
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
   elseif expr isa Symbol || (expr isa Expr && expr.head === :. &&
                              length(expr.args) == 2 && expr.args[2] isa QuoteNode)
      # Handle the following usage:
      #    @attributes Type
      #    @attributes Module.Type
      #    @attributes Module.Submodule.Type
      #    etc.
      return esc(quote
         # do nothing if the type already has storage for attributes
         if !AbstractAlgebra._is_attribute_storing_type($expr)
           isstructtype($expr) && AbstractAlgebra.ismutabletype($expr) || error("attributes can only be attached to mutable structs")
           let attr = Base.WeakKeyDict{$expr, Dict{Symbol, Any}}()
             AbstractAlgebra._get_attributes(G::$expr) = Base.get(attr, G, nothing)
             AbstractAlgebra._get_attributes!(G::$expr) = Base.get!(() -> Dict{Symbol, Any}(), attr, G)
             AbstractAlgebra._get_attributes(::Type{$expr}) = attr
             AbstractAlgebra._is_attribute_storing_type(::Type{$expr}) = true
           end
         end
      end)
   end
   error("attributes can only be attached to mutable structs")
end

_is_attribute_storing_type(::Type{T}) where T = Base.issingletontype(T) || isstructtype(T) && ismutable(T) && hasfield(T, :__attrs)

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
   isstructtype(T) && ismutable(T) && hasfield(T, :__attrs) || error("attributes storage not supported for type $T")
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
   isstructtype(T) && ismutable(T) && hasfield(T, :__attrs) || error("attributes storage not supported for type $T")
   if !isdefined(G, :__attrs)
      G.__attrs = Dict{Symbol, Any}()
   end
   return G.__attrs
end

"""
    has_attribute(G::Any, attr::Symbol)

Return a boolean indicating whether `G` has a value stored for the attribute `attr`.
"""
function has_attribute(G::Any, attr::Symbol)
   D = _get_attributes(G)
   return D isa Dict && haskey(D, attr)
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

"""
    get_attribute!(f::Function, G::Any, attr::Symbol)

Return the value stored for the attribute `attr` of `G, or if no value has been set,
store `key => f()` and return `f()`.
"""
function get_attribute!(f::Function, G::Any, attr::Symbol)
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
