###############################################################################
#
#   basic_interface.jl : basic interface for AbstractAlgebra
#
###############################################################################

# TODO: Move more generic functions to this file.

###############################################################################
#
#   Parents, elements and data type methods
#
###############################################################################

@doc raw"""
    parent(a)

Return parent object of given element $a$.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> G = SymmetricGroup(5); g = Perm([3,4,5,2,1])
(1,3,5)(2,4)

julia> parent(g) == G
true

julia> S, x = laurent_series_ring(ZZ, 3, "x")
(Laurent series ring in x over integers, x + O(x^4))

julia> parent(x) == S
true
```
"""
function parent end

# TODO: Give example
@doc raw"""
    elem_type(parent)
    elem_type(parent_type)

Given a parent object (or its type), return the type of its elements.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> S, x = power_series_ring(QQ, 2, "x")
(Univariate power series ring over rationals, x + O(x^3))

julia> elem_type(S) == typeof(x)
true
```
"""
elem_type(x)  = elem_type(typeof(x))
elem_type(T::DataType) = throw(MethodError(elem_type, (T,)))

@doc raw"""
    parent_type(element)
    parent_type(element_type)

Given an element (or its type), return the type of its parent object.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> R, x = polynomial_ring(ZZ, "x")
(Univariate polynomial ring in x over integers, x)

julia> S = matrix_space(R, 2, 2)
Matrix space of 2 rows and 2 columns
  over univariate polynomial ring in x over integers

julia> a = rand(S, 0:1, 0:1);

julia> parent_type(a) == typeof(S)
true
```
"""
parent_type(x) = parent_type(typeof(x))
parent_type(T::DataType) = throw(MethodError(parent_type, (T,)))

@doc raw"""
    base_ring(a)

Return base ring $R$ of given element or parent $a$.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> S, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> base_ring(S) == QQ
true

julia> R = GF(7)
Finite field F_7

julia> base_ring(R)
Union{}
```
"""
function base_ring end
base_ring(x::ModuleElem) = base_ring(parent(x))
base_ring(x::NCRingElement) = base_ring(parent(x))

# generic coefficient_ring method
coefficient_ring(x::NCRingElement) = coefficient_ring(parent(x))

###############################################################################
#
#   Special elements
#
###############################################################################

@doc raw"""
    one(a)

Return the multiplicative identity in the algebraic structure of $a$, which can
be either an element or parent.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> S = matrix_space(ZZ, 2, 2)
Matrix space of 2 rows and 2 columns
  over integers

julia> one(S)
[1   0]
[0   1]

julia> R, x = puiseux_series_field(QQ, 4, "x")
(Puiseux series field in x over rationals, x + O(x^5))

julia> one(x)
1 + O(x^4)

julia> G = GF(5)
Finite field F_5

julia> one(G)
1
```
"""
function one end

@doc raw"""
    zero(a)

Return the additive identity in the algebraic structure of $a$, which can be
either an element or parent.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> S = matrix_ring(QQ, 2)
Matrix ring of degree 2
  over rationals

julia> zero(S)
[0//1   0//1]
[0//1   0//1]

julia> R, x = polynomial_ring(ZZ, "x")
(Univariate polynomial ring in x over integers, x)

julia> zero(x^3 + 2)
0
```
"""
function zero end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

@doc raw"""
    isone(a)

Return true if $a$ is the multiplicative identity, else return false.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> S = matrix_space(ZZ, 2, 2); T = matrix_space(ZZ, 2, 3); U = matrix_space(ZZ, 3, 2);

julia> isone(S([1 0; 0 1]))
true

julia> isone(T([1 0 0; 0 1 0]))
false

julia> isone(U([1 0; 0 1; 0 0]))
false

julia> T, x = puiseux_series_field(QQ, 10, "x")
(Puiseux series field in x over rationals, x + O(x^11))

julia> isone(x), isone(T(1))
(false, true)
```
"""
function isone end

@doc raw"""
    iszero(a)

Return true if $a$ is the additative identity, else return false.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> T, x = puiseux_series_field(QQ, 10, "x")
(Puiseux series field in x over rationals, x + O(x^11))

julia> a = T(0)
O(x^10)

julia> iszero(a)
true
```
"""
function iszero end

###############################################################################
#
#   More generic functions
#
###############################################################################

# @doc raw"""
#     gen(a)

# Return element generating parent $a$.

# # Examples
# ```jldoctest; setup = :(using AbstractAlgebra)
# julia> S, x = laurent_polynomial_ring(QQ, "x")
# (Univariate Laurent Polynomial Ring in x over Rationals, x)

# julia> gen(S)
# x
# ```
# """
function gen end

# @doc raw"""
#     gens(a)

# Return elements generating parent $a$ in an array.

# # Examples
# ```jldoctest; setup = :(using AbstractAlgebra)
# ```
# """
function gens end

###############################################################################
#
#   Variable names
#
###############################################################################

"""
    const VarName = Union{Symbol, AbstractString, Char}

Types allowed when giving variable names.
"""
const VarName = Union{Symbol, AbstractString, Char}

###############################################################################
#
#   Unsafe functions
#
###############################################################################

@doc raw"""
    zero!(a)

Return the zero of `parent(a)`, possibly modifying the object `a` in the process.
"""
function zero! end

@doc raw"""
    add!(a, b, c)

Return `b + c`, possibly modifying the object `a` in the process.
"""
function add! end

@doc raw"""
    addeq!(a, b)

Return `a + b`, possibly modifying the object `a` in the process.
"""
function addeq! end

@doc raw"""
    sub!(a, b, c)

Return `b - c`, possibly modifying the object `a` in the process.
"""
function sub! end

@doc raw"""
    mul!(a, b, c)

Return `b*c`, possibly modifying the object `a` in the process.
"""
function mul! end
