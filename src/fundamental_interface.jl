###############################################################################
#
#   basic_interface.jl : basic interface for AbstractAlgebra
#
###############################################################################

###############################################################################
#
#   Parents, elements and data type methods
#
###############################################################################

@doc raw"""
    parent(a)

Return parent object of given element $a$.

# Examples
```jldoctest
julia> G = SymmetricGroup(5); g = Perm([3,4,5,2,1])
(1,3,5)(2,4)

julia> parent(g) == G
true

julia> S, x = laurent_series_ring(ZZ, 3, :x)
(Laurent series ring in x over integers, x + O(x^4))

julia> parent(x) == S
true
```
"""
function parent end

@doc raw"""
    elem_type(parent)
    elem_type(parent_type)

Given a parent object (or its type), return the type of its elements.

# Examples
```jldoctest
julia> S, x = power_series_ring(QQ, 2, :x)
(Univariate power series ring over rationals, x + O(x^3))

julia> elem_type(S) == typeof(x)
true
```
"""
elem_type(x) = elem_type(typeof(x))
elem_type(T::DataType) = throw(MethodError(elem_type, (T,)))
elem_type(T::Type{Union{}}) = throw(MethodError(elem_type, (T,)))

@doc raw"""
    parent_type(element)
    parent_type(element_type)

Given an element (or its type), return the type of its parent object.

# Examples
```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
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
parent_type(T::Type{Union{}}) = throw(MethodError(parent_type, (T,)))

@doc raw"""
    base_ring(a)

Return the internal base ring of the given element or parent $a$.

# Examples
```jldoctest
julia> S, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over rationals, x)

julia> base_ring(S) == QQ
true

julia> R = GF(7)
Finite field F_7

julia> base_ring(R)
ERROR: MethodError: no method matching base_ring(::AbstractAlgebra.GFField{Int64})
```
"""
function base_ring end
base_ring(x::ModuleElem) = base_ring(parent(x))
base_ring(x::NCRingElement) = base_ring(parent(x))

@doc raw"""
    base_ring_type(a)

Return the type of the internal base ring of the given element, element type, parent or parent type $a$.

# Examples
```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> base_ring_type(R) == typeof(base_ring(R))
true

julia> base_ring_type(zero(R)) == typeof(base_ring(zero(R)))
true

julia> base_ring_type(typeof(R)) == typeof(base_ring(R))
true

julia> base_ring_type(typeof(zero(R))) == typeof(base_ring(zero(R)))
true

julia> R = GF(7)
Finite field F_7

julia> base_ring_type(R)
Union{}
```
"""
base_ring_type(x) = base_ring_type(typeof(x))
base_ring_type(x::Type{<:NCRingElement}) = base_ring_type(parent_type(x))
base_ring_type(x::Type{<:ModuleElem}) = base_ring_type(parent_type(x))
base_ring_type(x::Type{<:Ideal}) = base_ring_type(parent_type(x))
base_ring_type(T::DataType) = throw(MethodError(base_ring_type, (T,)))
base_ring_type(T::Type{Union{}}) = throw(MethodError(base_ring_type, (T,)))

@doc raw"""
    coefficient_ring(a)

Return the coefficient ring of the given element or parent $a$.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over rationals, x)

julia> coefficient_ring(x^2+1) == QQ
true

julia> S, (z,w) = universal_polynomial_ring(QQ, [:z,:w])
(Universal polynomial ring over Rationals, UniversalRingElem{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}, Rational{BigInt}}[z, w])

julia> coefficient_ring(S) == QQ
true
```
"""
function coefficient_ring end
coefficient_ring(x::NCRingElement) = coefficient_ring(parent(x))

@doc raw"""
    coefficient_ring_type(a)

Return the type of the coefficient ring of the given element, element type, parent or parent type $a$.

# Examples
```jldoctest
julia> R, x = polynomial_ring(ZZ, :x)
(Univariate polynomial ring in x over integers, x)

julia> coefficient_ring_type(R) == typeof(coefficient_ring(R))
true

julia> coefficient_ring_type(zero(R)) == typeof(coefficient_ring(zero(R)))
true

julia> coefficient_ring_type(typeof(R)) == typeof(coefficient_ring(R))
true

julia> coefficient_ring_type(typeof(zero(R))) == typeof(coefficient_ring(zero(R)))
true
```
"""
coefficient_ring_type(x) = coefficient_ring_type(typeof(x))
coefficient_ring_type(x::Type{<:NCRingElement}) = coefficient_ring_type(parent_type(x))
coefficient_ring_type(x::Type{<:ModuleElem}) = coefficient_ring_type(parent_type(x))
coefficient_ring_type(x::Type{<:Ideal}) = coefficient_ring_type(parent_type(x))
coefficient_ring_type(T::DataType) = throw(MethodError(coefficient_ring_type, (T,)))
coefficient_ring_type(T::Type{Union{}}) = throw(MethodError(coefficient_ring_type, (T,)))

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
```jldoctest
julia> S = matrix_space(ZZ, 2, 2)
Matrix space of 2 rows and 2 columns
  over integers

julia> one(S)
[1   0]
[0   1]

julia> R, x = puiseux_series_field(QQ, 4, :x)
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
```jldoctest
julia> S = matrix_ring(QQ, 2)
Matrix ring of degree 2
  over rationals

julia> zero(S)
[0//1   0//1]
[0//1   0//1]

julia> R, x = polynomial_ring(ZZ, :x)
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
```jldoctest
julia> S = matrix_space(ZZ, 2, 2); T = matrix_space(ZZ, 2, 3); U = matrix_space(ZZ, 3, 2);

julia> isone(S([1 0; 0 1]))
true

julia> isone(T([1 0 0; 0 1 0]))
false

julia> isone(U([1 0; 0 1; 0 0]))
false

julia> T, x = puiseux_series_field(QQ, 10, :x)
(Puiseux series field in x over rationals, x + O(x^11))

julia> isone(x), isone(T(1))
(false, true)
```
"""
function isone end

@doc raw"""
    iszero(a)

Return true if $a$ is the additive identity, else return false.

# Examples
```jldoctest
julia> T, x = puiseux_series_field(QQ, 10, :x)
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
# ```jldoctest
# julia> S, x = laurent_polynomial_ring(QQ, :x)
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
# ```jldoctest
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

Return `zero(parent(a))`, possibly modifying the object `a` in the process.
"""
zero!(a) = zero(parent(a))

@doc raw"""
    one!(a)

Return `one(parent(a))`, possibly modifying the object `a` in the process.
"""
one!(a) = one(parent(a))

@doc raw"""
    neg!(z, a)
    neg!(a)

Return `-a`, possibly modifying the object `z` in the process.
Aliasing is permitted.
The unary version is a shorthand for `neg!(a, a)`.
"""
neg!(z, a) = -a
neg!(a) = neg!(a, a)

@doc raw"""
    inv!(z, a)
    inv!(a)

Return `AbstractAlgebra.inv(a)`, possibly modifying the object `z` in the process.
Aliasing is permitted.
The unary version is a shorthand for `inv!(a, a)`.

!!! note
    `AbstractAlgebra.inv` and `Base.inv` differ only in their behavior on julia
    types like `Integer` and `Rational{Int}`. The former makes it adhere to the
    Ring interface.
"""
inv!(z, a) = AbstractAlgebra.inv(a)
inv!(a) = inv!(a, a)

for (name, op) in ((:add!, :+), (:sub!, :-), (:mul!, :*), (:pow!, :^))
  @eval begin
    @doc """
        $($name)(z, a, b)
        $($name)(a, b)

    Return `a $($op) b`, possibly modifying the object `z` in the process.
    Aliasing is permitted.
    The two argument version is a shorthand for `$($name)(a, a, b)`.
    """
    $name(z, a, b) = $op(a, b)
    $name(a, b) = $name(a, a, b)
  end
end

@doc raw"""
    addmul!(z, a, b, t)
    addmul!(z, a, b)

Return `z + a * b`, possibly modifying the objects `z` and `t` in the process.

The second version is usually a shorthand for `addmul!(z, a, b, parent(z)())`,
but in some cases may be more efficient. For multiple operations in a row that
use temporary storage, it is still best to use the four argument version.
"""
addmul!(z, a, b, t) = add!(z, mul!(t, a, b))
addmul!(z, a, b) = addmul!(z, a, b, parent(z)())

@doc raw"""
    submul!(z, a, b, t)
    submul!(z, a, b)

Return `z - a * b`, possibly modifying the objects `z` and `t` in the process.

The second version is usually a shorthand for `submul!(z, a, b, parent(z)())`,
but in some cases may be more efficient. For multiple operations in a row that
use temporary storage, it is still best to use the four argument version.
"""
submul!(z, a, b, t) = sub!(z, mul!(t, a, b))
submul!(z, a, b) = submul!(z, a, b, parent(z)())

function divexact end

for name in (:divexact, :rem, :mod, :gcd, :lcm)
  name_bang = Symbol(name, "!")
  @eval begin
    @doc """
        $($name_bang)(z, a, b)
        $($name_bang)(a, b)

    Return `$($name)(a, b)`, possibly modifying the object `z` in the process.
    Aliasing is permitted.
    The two argument version is a shorthand for `$($name)(a, a, b)`.
    """
    $name_bang(z, a, b) = $name(a, b)
    $name_bang(a, b) = $name_bang(a, a, b)
  end
end

@doc """
    div!(z, a, b)
    div!(a, b)

Return `div(a, b)`, possibly modifying the object `z` in the process.
Aliasing is permitted.
The two argument version is a shorthand for `div(a, a, b)`.

!!! note
    `AbstractAlgebra.div` and `Base.div` differ only in their behavior on julia
    types like `Integer` and `Rational{Int}`. The former makes it adhere to the
    Ring interface.
"""
div!(z, a, b) = AbstractAlgebra.div(a, b)
div!(a, b) = div!(a, a, b)

@doc """
    divides!(z, a, b)

Return `divides(a, b)`, possibly modifying the object `z` in the process.
"""
divides!(z, a, b) = divides(a, b)

@doc raw"""
    canonical_injection(D, i)

Return the i-th canonical injection into the direct sum or product objects `D`.
"""
function canonical_injection end

@doc raw"""
    canonical_projection(D, i)

Return the i-th canonical projection into the direct sum or product objects `D`.
"""
function canonical_projection end

@doc raw"""
    _number_of_direct_product_factors(D)

Return the number of factors/ summands in the direct product/ sum object `D`
"""
function _number_of_direct_product_factors end

@doc raw"""
    hom(D, C, data)

Return the homomorphism from the domain `D` into the codomain `C` defined by the data.
"""
function hom end

function can_solve_with_solution end
