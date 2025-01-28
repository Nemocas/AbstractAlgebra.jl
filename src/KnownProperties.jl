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
julia> AbstractAlgebra.is_known(5, is_even)
true
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

Return whether the property of `x` called for by `f(x, args...)` is known. 

Note: The default implementation throws an error. It is the programmer's 
responsibility to implement appropriate methods for their individual 
types and properties. See `src/KnownProperties.jl` for details.
"""
function is_known(x::Any, f::Function, args...; kwargs...)
  return _is_known(x, f, args...; kwargs...)
end

function _is_known(x::Any, f::Function, args...; kwargs...)
  error("no method implemented to check whether property $(nameof(f)) with arguments $(args) and keyword arguments $(kwargs) is known for object $x")
end


# In general, it is the programmer's responsibility to implement 
# a method for their types and properties. The above already implements 
# a generic deflection to the attributes where applicable. But this 
# will by far not catch all cases. 
# 
# We give a few samples. 

is_known(i::Int, ::typeof(is_even)) = true # Almost no computation cost, so known. 
is_known(i::Int, ::typeof(is_odd)) = true # Almost no computation cost, so known. 

# The dimension of a ring is difficult to compute in general. However, for 
# some rings it is easy.
is_known(R::MPolyRing{<:FieldElem}, ::typeof(dim)) = true

# The following needs data types which are only available in OSCAR. However, 
# we put the code here to illustrate what an implementation could look like.
#
# function is_known(A::MPolyQuoRing, ::typeof(dim))
#   return is_known(modulus(A), dim)
# end
#
# The `MPolyIdeal` caches information about its dimension in a field, not 
# the attributes via `@attr`. Hence, we check there. 
#
# function is_known(I::MPolyIdeal, ::typeof(dim))
#   return I.dim !== nothing
# end
