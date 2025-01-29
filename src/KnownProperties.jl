########################################################################
# This file introduces a common framework to ask whether or not a given
# property (which is potentially hard to compute) is known for some
# object `x`. 
########################################################################

@doc raw"""
    is_known(x::Any, f::Function)

Given a unary function `f` and a value `x`, return whether `f(x)` is "known"
in the sense that evaluating `f(x)` is fast and takes constant time.

For example this might return `true` if `f(x)` was computed before and has
been cached.

Note: The default implementation only checks whether an attribute 
(from the [`@attr`](@ref) macro) with the same name as the function `f` exists. 
Otherwise, it throws an error. In general for `is_known` to work correctly
for a given function `f` requires that everyone adding, modifying or removing
methods for `f` adjusts it as needed. That is, they are responsible for
ensuring `is_known(x,f)` returns an appropriate result whenever `f(x)` would
invoke the method they just modified. For example by adding a new method
`is_known(::MyTypes,::Type{f}) = true` where `MyTypes` is the union of types
accepted by the method for `f`.

Conversely this means that `is_known` does not work for arbitrary unary
function `f` but instead only a select supported list.

See `src/KnownProperties.jl` in Julia package `AbstractAlgebra.jl` for details.

# Examples
```jldoctest
julia> AbstractAlgebra.is_known(::Int, ::typeof(is_even)) = true; # manual installation of the method

julia> AbstractAlgebra.is_known(5, is_even) # sample call to the function
true

julia> AbstractAlgebra.is_known(::MPolyRing{<:FieldElem}, ::typeof(dim)) = true; # another implementation

julia> AbstractAlgebra.is_known(R::MPolyRing, ::typeof(dim)) = AbstractAlgebra.is_known(coefficient_ring(R), dim) # generic deflection to the `coefficient_ring`

julia> R, (x, y) = ZZ[:x, :y];

julia> try
         AbstractAlgebra.is_known(R, dim)
       catch e
         e
       end
ErrorException("no method implemented to check whether property dim is known for object Integers")
```
"""
function is_known(x::Any, f::Function)
  return _is_known(x, f)
end

function _is_known(x::Any, f::Function)
  # If the object `x` has attributes and an attribute with the name 
  # of `f` happens to be stored, return it in good faith.
  _is_attribute_storing_type(typeof(x)) && has_attribute(x, nameof(f)) && return true
  # Otherwise notify the user that no method is implemented rather than 
  # just saying "no". The explicit call to this method suggests that 
  # the callee expects the function to return something useful and 
  # not having a working method at hand is probably a bug. 
  error("no method implemented to check whether property $(nameof(f)) is known for object $x")
end

@doc raw"""
    is_known(x::Any, f::Function, args...; kwargs...)

Given a function `f` on `(x, args...; kwargs)` and a value `x`, 
return whether `f(x, args...; kwargs...)` is "known" in the sense that 
evaluating `f(x, args...; kwargs...)` is fast and takes constant time.

For example this might return `true` if the result was computed before and has
been cached.

Note: The default implementation only throws an error. In general for `is_known` 
to work correctly for a given function `f` requires that everyone adding,
modifying or removing methods for `f` adjusts it as needed.

See `src/KnownProperties.jl` in Julia package `AbstractAlgebra.jl` for details.
"""
function is_known(x::Any, f::Function, args...; kwargs...)
  return _is_known(x, f, args...; kwargs...)
end

function _is_known(x::Any, f::Function, args...; kwargs...)
  error("no method implemented to check whether property $(nameof(f)) with arguments $(args) and keyword arguments $(kwargs) is known for object $x")
end

