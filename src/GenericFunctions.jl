###############################################################################
#
#   GenericFunctions.jl : generic functions
#
###############################################################################

###############################################################################
#
#   Parents and elements
#
###############################################################################

@doc Markdown.doc"""
    parent(a)

Return parent object of given element $a$.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> G = SymmetricGroup(5); g = Perm([3,4,5,2,1])
(1,3,5)(2,4)

julia> parent(g) == G
true

julia> S, x = LaurentSeriesRing(ZZ, 3, "x")
(Laurent series ring in x over Integers, x + O(x^4))

julia> parent(x) == S
true
```
"""
function parent end

@doc Markdown.doc"""
    elem_type(parent)
    elem_type(parent_type)
Given a parent object (or its type), return the type of its elements.
"""
elem_type(x)  = elem_type(typeof(x))
elem_type(T::DataType) = throw(MethodError(elem_type, (T,)))

@doc Markdown.doc"""
    parent_type(element)
    parent_type(element_type)
Given an element (or its type), return the type of its parent object.
"""
parent_type(x) = parent_type(typeof(x))
parent_type(T::DataType) = throw(MethodError(parent_type, (T,)))

