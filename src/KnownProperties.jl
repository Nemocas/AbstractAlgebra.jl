########################################################################
# This file introduces a common framework to ask whether or not a given
# property (which is potentially hard to compute) is known for some
# object `x`. 
########################################################################

@doc raw"""
    is_known(f::Function, args...; kwargs...)

Given a function `f`, arguments `args`, and potentially keyword arguments 
`kwargs`, return whether `f(args...; kwargs...)` is "known" in the sense that 
evaluating `f(args...; kwargs...)` is fast and takes constant time.

For example this might return `true` if `f(args...; kwargs...)` was computed 
before and has been cached.

Note: There is no default implementation. In general for `is_known` to work correctly
for a given function `f` requires that everyone adding, modifying or removing
methods for `f` adjusts the corresponding methods for `is_known(::typeof(f), ...)` 
as needed. That is, they are responsible for ensuring `is_known(f, args...; kwargs...)` 
returns an appropriate result whenever `f(args...; kwargs...)` would invoke the method 
they just added or modified, resp. would have invoked a method they removed.

Conversely this means that `is_known` does not work for an arbitrary
function `f` but instead only for a select supported list.

# Examples
```jldoctest
julia> AbstractAlgebra.is_known(::typeof(is_even), ::Int) = true; # manual installation of the method

julia> AbstractAlgebra.is_known(is_even, 5) # sample call to the function
true

julia> AbstractAlgebra.is_known(::typeof(dim), ::MPolyRing{<:FieldElem}) = true; # another implementation

julia> AbstractAlgebra.is_known(::typeof(dim), R::MPolyRing) = AbstractAlgebra.is_known(dim, coefficient_ring(R)) # generic deflection to the `coefficient_ring`

julia> R, (x, y) = ZZ[:x, :y];

julia> AbstractAlgebra.is_known(dim, R)
ERROR: MethodError: no method matching is_known(::typeof(dim), ::AbstractAlgebra.Integers{BigInt})
```
"""
function is_known end

# A generic implementation to look up attributes. This can be used to implement 
# `is_known` for a specific type of arguments of unary functions.
function _is_known_via_attributes(f::Function, x::Any)
  @req _is_attribute_storing_type(typeof(x)) "objects of type $(typeof(x)) do not support attribute storage"
  return has_attribute(x, nameof(f))
end

