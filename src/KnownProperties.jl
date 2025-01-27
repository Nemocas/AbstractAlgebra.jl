########################################################################
# This is to introduce a common language to ask whether or not a given
# property (which is potentially hard to compute) is known for some
# object `x`. 
########################################################################

@doc raw"""
    is_known(x::Any, f::Function)
    is_known(x::Any, f::Function, v::Any)

Return whether a property `p` of an object `x` is known. 

Here `f` is function which would usually be called to compute that 
property. The second signature can be used to check whether `p` is 
known and its value equals the third argument `v`.

Note: The default implementation throws an error. It is the programmer's 
responsibility to implement appropriate methods for their individual 
types and properties. See `src/KnownProperties.jl` for details.
"""
function is_known(x::Any, f::Function)
  # If the object `x` has attributes and an attribute with the name 
  # of `f` happens to be stored, return it in good faith.
  hasfield(typeof(x), :__attrs) && has_attribute(x, nameof(f)) && return get_attribute(x, nameof(f))
  # Otherwise notify the user that no method is implemented rather than 
  # just saying "no". The explicit call to this method suggests that 
  # the callee expects the function to return something useful and 
  # not having a working method at hand is probably a bug. 
  error("no method implemented to check whether property $(nameof(f)) is known for object $x")
end

function is_known(x::Any, f::Function, p::Any)
  return is_known(x, f) && f(x) == p
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
