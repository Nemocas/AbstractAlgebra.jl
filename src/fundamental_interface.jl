###############################################################################
#
#   basic_interface.jl : basic interface for AbstractAlgebra
#
###############################################################################

# TODO: Move more generic functions to this file.

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

# TODO: Give example
@doc Markdown.doc"""
    elem_type(parent)
    elem_type(parent_type)

Given a parent object (or its type), return the type of its elements.

# Example
```jldoctest; setup = :(using AbstractAlgebra)
julia> S, x = PowerSeriesRing(QQ, 2, "x")
(Univariate power series ring in x over Rationals, x + O(x^3))

julia> elem_type(S) == typeof(x)
true
```
"""
elem_type(x)  = elem_type(typeof(x))
elem_type(T::DataType) = throw(MethodError(elem_type, (T,)))

@doc Markdown.doc"""
    parent_type(element)
    parent_type(element_type)

Given an element (or its type), return the type of its parent object.

# Example
```jldoctest; setup = :(using AbstractAlgebra)
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S = MatrixSpace(R, 2, 2)
Matrix Space of 2 rows and 2 columns over Univariate Polynomial Ring in x over Integers

julia> a = rand(S, 0:1, 0:1);

julia> parent_type(a) == typeof(S)
true
```
"""
parent_type(x) = parent_type(typeof(x))
parent_type(T::DataType) = throw(MethodError(parent_type, (T,)))

###############################################################################
#
#   Miscellaneous
#
###############################################################################

@doc Markdown.doc"""
    base_ring(a::T)

Return base ring $R$ of given element or parent $a$.

# Example
```jldoctest; setup = :(using AbstractAlgebra)
julia> S, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> base_ring(S) == QQ
true

julia> R = GF(7)
Finite field F_7

julia> base_ring(R)
Union{}
```
"""
function base_ring end

@doc Markdown.doc"""
    gen(a::T)

Return element generating parent $a$.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> S, x = LaurentPolynomialRing(QQ, "x")
(Univariate Laurent Polynomial Ring in x over Rationals, x)

julia> gen(S)
x

julia> G = GF(7)
Finite field F_7

julia> gen(G)
1
```
"""
function gen end

@doc Markdown.doc"""
    gens(a::T)

Return elements generating parent $a$ in an array.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> S, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Integers, AbstractAlgebra.Generic.MPoly{BigInt}[x, y, z])

julia> gens(S)
3-element Array{AbstractAlgebra.Generic.MPoly{BigInt},1}:
 x
 y
 z
```
"""
function gens end

# TODO: Move from CRT from Hecke to AbstractAlgebra?
# Currently no implementation, only example on how the arbitrary inputs `crt`
# should look like.
@doc Markdown.doc"""
    crt(r::Array{T, 1}, m::Array{T, 1}) where T
    crt(r::T, m::T...) where T

Return $x$ in the Euclidean domain $T$ such that $x \equiv r_i \mod m_i$
for all $i$.
"""
function crt end

@doc Markdown.doc"""
    factor(a::T, b::R)

Return factorization of element $a$. 
"""
function factor end

@doc Markdown.doc"""
    factor_squarefree(a::T)

Return square free factorization of element $a$.
"""
function factor_squarefree end
