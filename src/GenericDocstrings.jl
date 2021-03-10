###############################################################################
#
#   GenericDocstrings.jl : generic docstrings
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
