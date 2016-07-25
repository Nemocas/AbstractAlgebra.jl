


<a id='Introduction-1'></a>

## Introduction


Nemo allows the creation of fraction fields over any ring $R$. We don't require $R$ to be an integral domain, however no attempt is made to deal with the general case. Two fractions $a/b$ and $c/d$ are equal in Nemo iff $ad = bc$. Thus, in practice, a greatest common divisor function is currently required for the ring $R$.


In order to make the representation $a/b$ unique for printing, we have a notion of canonical unit for elements of a ring $R$. When canonicalising $a/b$, each of the elements $a$ and $b$ is first divided by the canonical unit of $b$.


The `canonical_unit` function is defined for elements of every Nemo ring. It must have the properties


```
canonical_unit(u) == u
canonical_unit(a*b) == canonical_unit(a)*canonical_unit(b)
```


for any unit $u$ of the ring in question, and $a$ and $b$ arbitrary elements of the ring.


For example, the canonical unit of an integer is its sign. Thus a fraction of integers always has positive denominator after canonicalisation.


The canonical unit of a polynomial is the canonical unit of its leading coefficient, etc.


There are two different kinds of implementation of fraction fields in Nemo: a generic one for the case where no specific implementation exists, and efficient implementations of fractions over specific rings, usually provided by C/C++ libraries.


The following table shows each of the fraction types available in Nemo, the base ring $R$, and the Julia/Nemo types for that kind of fraction (the type information is mainly of concern to developers).


|        Base ring | Library | Element type |          Parent type |
| ----------------:| -------:| ------------:| --------------------:|
| Generic ring $R$ |    Nemo | `GenFrac{T}` |    `GenFracField{T}` |
|     $\mathbb{Z}$ |   Flint |       `fmpq` | `FlintRationalField` |


All fraction element types belong to the abstract type `FracElem` and all of the fraction field types belong to the abstract type `FracField`. This enables one to write generic functions that can accept any Nemo fraction type.


<a id='Fraction-field-constructors-1'></a>

## Fraction field constructors


In order to construct fractions in Nemo, one must first construct the fraction field itself. This can be accomplished manually with the following constructor.

<a id='Nemo.FractionField-Tuple{Nemo.Ring,Bool}' href='#Nemo.FractionField-Tuple{Nemo.Ring,Bool}'>#</a>
**`Nemo.FractionField`** &mdash; *Method*.



Note: it is also possible to construct fractions directly in Nemo without manually constructing the fraction field. In such cases, Nemo creates the relevant fraction field internally.


For convenience, the rational fraction field is constructed automatically in Nemo. We have the definition


```
QQ = FractionField(ZZ)
```


Here are some examples of creating fraction fields and making use of the resulting parent objects to coerce various elements into those fields.


```
R = FractionField(ZZ)
S, x = PolynomialRing(ZZ, "x")
T = FractionField(S)

a = R(1)
b = T(fmpz(3))
c = T(x)
d = QQ(11)
```


<a id='Fraction-constructors-1'></a>

## Fraction constructors


Once a fraction field is constructed, there are various ways to construct fractions in that field.


Apart from coercing elements into the fraction field as above, we offer the following functions.

<a id='Base.zero-Tuple{Nemo.FracField}' href='#Base.zero-Tuple{Nemo.FracField}'>#</a>
**`Base.zero`** &mdash; *Method*.



```
zero(R::FracField)
```

> Return $0/1$ in the given fraction field.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/generic/Fraction.jl#L92' class='documenter-source'>source</a><br>

<a id='Base.one-Tuple{Nemo.FracField}' href='#Base.one-Tuple{Nemo.FracField}'>#</a>
**`Base.one`** &mdash; *Method*.



```
one(R::FracField)
```

> Return $1/1$ in the given fraction field.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/generic/Fraction.jl#L98' class='documenter-source'>source</a><br>


It is also possible to construct fractions directly in Nemo, without first manually constructing the relevant fraction field. For this purpose we overload Julia's fraction operator `//`.


Here are some examples of constructing fractions.


```
S, x = PolynomialRing(ZZ, "x")
T = FractionField(S)

a = zero(T)
b = one(T)
c = (x + 3)//(x^2 + 2)
```


<a id='Basic-functionality-1'></a>

## Basic functionality


All fraction field modules in Nemo must provide the functionality listed in this section. (Note that only some of these functions are useful to a user.)


Developers who are writing their own fraction field module, whether as an interface to a C library, or as some kind of generic module, must provide all of these functions for custom fraction field types in Nemo. 


We write `U` for the type of fraction elements in the fraction field and `T` for the type of elements of the base ring.


All of these functions are provided for all existing fraction types in Nemo.


```
parent_type{U <: FracElem}(::Type{U})
```


Given the type of fraction elements, should return the type of the corresponding parent object.


```
elem_type(R::FracField)
```


Given a parent object for the fraction field, return the type of elements of the fraction field.


```
Base.hash(a::FracElem, h::UInt)
```


Return a `UInt` hexadecimal hash of the fraction $a$. This should be xor'd with a fixed random hexadecimal specific to the fraction type. The hash of the numerator and denominator of a fraction should be xor'd with the supplied parameter `h` as part of computing the hash.


```
num(a::FracElem)
```


Return the numerator of the given fraction element, i.e. for $a = b/c$ return $b$. The returned numerator will be divided by the canonical unit of the denominator.


```
den(a::FracElem)
```


Return the denominator of the given fraction element, i.e. for $a = b/c$ return $c$. The returned denominator will be divided by the canonical unit of the denominator.


```
deepcopy(a::FracElem)
```


Construct a copy of the given fraction and return it. This function must recursively construct copies of all of the internal data in the given fraction. Nemo fractions are mutable and so returning shallow copies is not sufficient.


```
mul!(c::FracElem, a::FracElem, b::FracElem)
```


Multiply $a$ by $b$ and set the existing fraction $c$ to the result. This function is provided for performance reasons as it saves allocating a new object for the result and eliminates associated garbage collection.


```
addeq!(c::FracElem, a::FracElem)
```


In-place addition. Adds $a$ to $c$ and sets $c$ to the result. This function is provided for performance reasons as it saves allocating a new object for the result and eliminates associated garbage collection.


Given a parent object `S` for a fraction field, the following coercion functions are provided to coerce various elements into the fraction field.  Developers provide these by overloading the `call` operator for the fraction field parent objects.


```
S()
```


Coerce zero into the field $S$.


```
S(n::Integer)
S(n::fmpz)
```


Coerce an integer value or Flint integer into the fraction field $S$.


```
S(n::T)
```


Coerces an element of the base ring, of type `T` into $S$.


```
S(f::FracElem)
```


Take a fraction that is already in the field $S$ and simply return it. A copy of the original is not made.


```
S(c::RingElem)
```


Try to coerce the given ring element into the fraction field. This only succeeds if $c$ can be coerced into the base ring.


There are also the followin constructors for creatinf fracions from a numerator and denominator.


```
S(n::T, d::T)
S(n::Integer, d::T)
S(n::T, d::Integer)
```


Create the fraction $n/d$ in the fraction field $S$.


In addition to the above, developers of custom fractions must ensure the parent object of a fraction type constains a field `base_ring` specifying the base ring. They must also ensure that each fraction element contains a field `parent` specifying the parent object of the fraction.


Typically a developer will also overload the `FractionField` generic function to create fractions of the custom type they are implementing.


<a id='Basic-manipulation-1'></a>

## Basic manipulation


Numerous functions are provided to manipulate fractions. Also see the section on basic functionality above.

<a id='Nemo.base_ring-Tuple{Nemo.FracField}' href='#Nemo.base_ring-Tuple{Nemo.FracField}'>#</a>
**`Nemo.base_ring`** &mdash; *Method*.



```
base_ring{T}(S::FracField{T})
```

> Return the base ring $R$ of the given fraction field.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/generic/Fraction.jl#L19' class='documenter-source'>source</a><br>

<a id='Nemo.base_ring-Tuple{Nemo.FracElem}' href='#Nemo.base_ring-Tuple{Nemo.FracElem}'>#</a>
**`Nemo.base_ring`** &mdash; *Method*.



```
base_ring{T}(r::FracElem)
```

> Return the base ring $R$ of the fraction field that the supplied element $a$ belongs to.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/generic/Fraction.jl#L25' class='documenter-source'>source</a><br>

<a id='Base.parent-Tuple{Nemo.FracElem}' href='#Base.parent-Tuple{Nemo.FracElem}'>#</a>
**`Base.parent`** &mdash; *Method*.



```
parent(a::FracElem)
```

> Return the parent object of the given fraction element.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/generic/Fraction.jl#L32' class='documenter-source'>source</a><br>

<a id='Nemo.iszero-Tuple{Nemo.FracElem}' href='#Nemo.iszero-Tuple{Nemo.FracElem}'>#</a>
**`Nemo.iszero`** &mdash; *Method*.



```
iszero(a::FracElem)
```

> Return `true` if the supplied element $a$ is zero in the fraction field it belongs to, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/generic/Fraction.jl#L104' class='documenter-source'>source</a><br>

<a id='Nemo.isone-Tuple{Nemo.FracElem}' href='#Nemo.isone-Tuple{Nemo.FracElem}'>#</a>
**`Nemo.isone`** &mdash; *Method*.



```
isone(a::FracElem)
```

> Return `true` if the supplied element $a$ is one in the fraction field it belongs to, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/generic/Fraction.jl#L111' class='documenter-source'>source</a><br>

<a id='Nemo.isunit-Tuple{Nemo.FracElem}' href='#Nemo.isunit-Tuple{Nemo.FracElem}'>#</a>
**`Nemo.isunit`** &mdash; *Method*.



```
isunit(a::FracElem)
```

> Return `true` if the supplied element $a$ is invertible in the fraction field it belongs to, i.e. the numerator is nonzero, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/generic/Fraction.jl#L118' class='documenter-source'>source</a><br>


Some functions are only available for certain rings.

<a id='Base.abs-Tuple{Nemo.fmpq}' href='#Base.abs-Tuple{Nemo.fmpq}'>#</a>
**`Base.abs`** &mdash; *Method*.



```
abs(a::fmpq)
```

> Return the absolute value of $a$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L82' class='documenter-source'>source</a><br>

<a id='Nemo.height-Tuple{Nemo.fmpq}' href='#Nemo.height-Tuple{Nemo.fmpq}'>#</a>
**`Nemo.height`** &mdash; *Method*.



```
height(a::fmpq)
```

> Return the height of the fraction $a$, namely the largest of the absolute values of the numerator and denominator.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L102' class='documenter-source'>source</a><br>

<a id='Nemo.height_bits-Tuple{Nemo.fmpq}' href='#Nemo.height_bits-Tuple{Nemo.fmpq}'>#</a>
**`Nemo.height_bits`** &mdash; *Method*.



```
height_bits(a::fmpq)
```

> Return the number of bits of the height of the fraction $a$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L113' class='documenter-source'>source</a><br>

<a id='Base.<<-Tuple{Nemo.fmpq,Int64}' href='#Base.<<-Tuple{Nemo.fmpq,Int64}'>#</a>
**`Base.<<`** &mdash; *Method*.



```
<<(x, n)
```

Left bit shift operator, `x << n`. For `n >= 0`, the result is `x` shifted left by `n` bits, filling with `0`s. This is equivalent to `x * 2^n`. For `n < 0`, this is equivalent to `x >> -n`.

```jldoctest
julia> Int8(3) << 2
12

julia> bits(Int8(3))
"00000011"

julia> bits(Int8(12))
"00001100"
```

See also [`>>`](:func:`>>`), [`>>>`](:func:`>>>`).


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/operators.jl#L134-152' class='documenter-source'>source</a><br>


```
<<(a::fmpq, b::Int)
```

> Return $2^b\times a$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L359' class='documenter-source'>source</a><br>

<a id='Base.>>-Tuple{Nemo.fmpq,Int64}' href='#Base.>>-Tuple{Nemo.fmpq,Int64}'>#</a>
**`Base.>>`** &mdash; *Method*.



```
>>(x, n)
```

Right bit shift operator, `x >> n`. For `n >= 0`, the result is `x` shifted right by `n` bits, where `n >= 0`, filling with `0`s if `x >= 0`, `1`s if `x < 0`, preserving the sign of `x`. This is equivalent to `fld(x, 2^n)`. For `n < 0`, this is equivalent to `x << -n`.

```jldoctest
julia> Int8(13) >> 2
3

julia> bits(Int8(13))
"00001101"

julia> bits(Int8(3))
"00000011"

julia> Int8(-14) >> 2
-4

julia> bits(Int8(-14))
"11110010"

julia> bits(Int8(-4))
"11111100"
```

See also [`>>>`](:func:`>>>`), [`<<`](:func:`<<`).


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/operators.jl#L161-190' class='documenter-source'>source</a><br>


```
<<(a::fmpq, b::Int)
```

> Return $2^b/a$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L348' class='documenter-source'>source</a><br>


Rational fractions can be compared with each other and with integers. Julia provides the full range of operators $<, >, \leq, \geq$ which depend on the following functions.

<a id='Base.isless-Tuple{Nemo.fmpq,Nemo.fmpq}' href='#Base.isless-Tuple{Nemo.fmpq,Nemo.fmpq}'>#</a>
**`Base.isless`** &mdash; *Method*.



```
isless(a::fmpq, b::fmpq)
```

> Return `true` if $a < b$, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L261' class='documenter-source'>source</a><br>

<a id='Base.isless-Tuple{Integer,Nemo.fmpq}' href='#Base.isless-Tuple{Integer,Nemo.fmpq}'>#</a>
**`Base.isless`** &mdash; *Method*.



```
isless(a::Integer, b::fmpq)
```

> Return `true` if $a < b$, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L299' class='documenter-source'>source</a><br>

<a id='Base.isless-Tuple{Nemo.fmpq,Integer}' href='#Base.isless-Tuple{Nemo.fmpq,Integer}'>#</a>
**`Base.isless`** &mdash; *Method*.



```
isless(a::fmpq, b::Integer)
```

> Return `true` if $a < b$, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L289' class='documenter-source'>source</a><br>

<a id='Base.isless-Tuple{Nemo.fmpq,Nemo.fmpz}' href='#Base.isless-Tuple{Nemo.fmpq,Nemo.fmpz}'>#</a>
**`Base.isless`** &mdash; *Method*.



```
isless(a::fmpq, b::fmpz)
```

> Return `true` if $a < b$, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L309' class='documenter-source'>source</a><br>

<a id='Base.isless-Tuple{Nemo.fmpz,Nemo.fmpq}' href='#Base.isless-Tuple{Nemo.fmpz,Nemo.fmpq}'>#</a>
**`Base.isless`** &mdash; *Method*.



```
isless(a::fmpz, b::fmpq)
```

> Return `true` if $a < b$, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L319' class='documenter-source'>source</a><br>


Here are some examples of basic manipulation of fractions.


```
S, x = PolynomialRing(ZZ, "x")

f = zero(S)
g = one(S)

a = isunit((x + 1)//(-x^2 + 1))
b = iszero(f)
c = isone(g)
d = abs(ZZ(11)//3)
U = base_ring(S)
V = base_ring((x + 1)//(-x^2 + 1))
W = parent(S(x))
4 <= ZZ(7)//ZZ(3)
```


<a id='Arithmetic-operators-1'></a>

## Arithmetic operators


All the usual arithmetic operators are overloaded for Nemo fractions. Note that Julia uses the single slash for floating point division. Therefore to perform exact division in a ring we use `divexact`. To construct an element of a fraction field one can use the double slash operator `//`.


The following operators and functions are provided.


|                                                  Function |      Operation |
| ---------------------------------------------------------:| --------------:|
|                                          `-(a::FracElem)` |    unary minus |
|        `+{T <: RingElem}(a::FracElem{T}, b::FracElem{T})` |       addition |
|        `-{T <: RingElem}(a::FracElem{T}, b::FracElem{T})` |    subtraction |
|        `*{T <: RingElem}(a::FracElem{T}, b::FracElem{T})` | multiplication |
| `divexact{T <: RingElem}(a::FracElem{T}, b::FracElem{T})` | exact division |


The following ad hoc operators are also provided.


|                                        Function |      Operation |
| -----------------------------------------------:| --------------:|
|                    `+(a::Integer, b::FracElem)` |       addition |
|                    `+(a::FracElem, b::Integer)` |       addition |
|                       `+(a::fmpz, b::FracElem)` |       addition |
|                       `+(a::FracElem, b::fmpz)` |       addition |
|        `+{T <: RingElem}(a::T, b::FracElem{T})` |       addition |
|        `+{T <: RingElem}(a::FracElem{T}, b::T)` |       addition |
|                    `-(a::Integer, b::FracElem)` |    subtraction |
|                    `-(a::FracElem, b::Integer)` |    subtraction |
|                       `-(a::fmpz, b::FracElem)` |    subtraction |
|                       `-(a::FracElem, b::fmpz)` |    subtraction |
|        `-{T <: RingElem}(a::T, b::FracElem{T})` |    subtraction |
|        `-{T <: RingElem}(a::FracElem{T}, b::T)` |    subtraction |
|                    `*(a::Integer, b::FracElem)` | multiplication |
|                    `*(a::FracElem, b::Integer)` | multiplication |
|                       `*(a::fmpz, b::FracElem)` | multiplication |
|                       `*(a::FracElem, b::fmpz)` | multiplication |
|        `*{T <: RingElem}(a::T, b::FracElem{T})` | multiplication |
|        `*{T <: RingElem}(a::FracElem{T}, b::T)` | multiplication |
|             `divexact(a::Integer, b::FracElem)` | exact division |
|             `divexact(a::FracElem, b::Integer)` | exact division |
|                `divexact(a::fmpz, b::FracElem)` | exact division |
|                `divexact(a::FracElem, b::fmpz)` | exact division |
| `divexact{T <: RingElem}(a::T, b::FracElem{T})` | exact division |
| `divexact{T <: RingElem}(a::FracElem{T}, b::T)` | exact division |
|                        `^(a::FracElem, n::Int)` |       powering |


If the appropriate `promote_rule` and coercion exists, these operators can also be used with elements of other rings. Nemo will try to coerce the operands to the dominating type and then apply the operator.


Here are some examples of arithmetic operations on fractions.


```
S, x = PolynomialRing(ZZ, "x")

a = -((x + 1)//(-x^2 + 1))
b = -(x + 3)//(x + 1) + (2x + 3)//(x^2 + 4)
c = (x + 1)//(-x^2 + 1) - x//(2x + 1)
d = ((x^2 + 3x)//(5x))*((x + 1)//(2x^2 + 2))
f = a + 2
g = 3 - a
h = b*(x + 1)
k = divexact(a, b)
m = divexact(a, x + 1)
n = divexact(b, 23)
p = divexact(ZZ(2), b)
q = a^3
```


<a id='Comparison-operators-1'></a>

## Comparison operators


The following comparison operators are implemented for fractions in Nemo. Julia provides the corresponding `!=` operator automatically.


<a id='Function-1'></a>

## Function


`isequal{T <: RingElem}(a::FracElem{T}, b::FracElem{T})` `=={T <: RingElem}(a::FracElem{T}, b::FracElem{T})`


The `isequal` operation returns `true` if and only if numerator and denominator of the fraction are precisely equal as compared by `isequal`. This is a stronger form of equality, used for comparing inexact ring elements, such as elements of a power series ring, the $p$-adics, or the reals or complex numbers. Two elements are precisely equal only if they have the same precision or bounds in addition to being arithmetically equal. 


In addition we have the following ad hoc comparison operators.


<a id='Function-2'></a>

## Function


`=={T <: RingElem}(a::FracElem{T}, b::T)` `=={T <: RingElem}(a::T, b::FracElem{T})` `==(a::FracElem, b::Integer)` `==(a::Integer, b::FracElem)` `==(a::FracElem, b::fmpz)` `==(a::fmpz, b::FracElem)`


Here are some examples of comparisons.


```
S, x = PolynomialRing(ZZ, "x")

a = -(ZZ(4)//ZZ(6))
b = -((x + 1)//(-x^2 + 1))

a == -ZZ(2)//ZZ(3)
b == 1//(x - 1)
a == 4
ZZ(2) == b
b == x + 1
```


<a id='Inversion-1'></a>

## Inversion

<a id='Base.inv-Tuple{Nemo.FracElem}' href='#Base.inv-Tuple{Nemo.FracElem}'>#</a>
**`Base.inv`** &mdash; *Method*.



```
inv(a::FracElem)
```

> Return the inverse of the fraction $a$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/generic/Fraction.jl#L513' class='documenter-source'>source</a><br>


Here are some examples of computing inverses.


```
S, x = PolynomialRing(ZZ, "x")

a = (x + 1)//(-x^2 + 1)

b = inv(a)
```


<a id='Greatest-common-divisor-1'></a>

## Greatest common divisor

<a id='Base.gcd-Tuple{Nemo.FracElem{T<:Nemo.RingElem},Nemo.FracElem{T<:Nemo.RingElem}}' href='#Base.gcd-Tuple{Nemo.FracElem{T<:Nemo.RingElem},Nemo.FracElem{T<:Nemo.RingElem}}'>#</a>
**`Base.gcd`** &mdash; *Method*.



```
gcd{T <: RingElem}(a::FracElem{T}, b::FracElem{T})
```

> Return a greatest common divisor of $a$ and $b$ if one exists. N.B: we define the GCD of $a/b$ and $c/d$ to be gcd$(ad, bc)/bd$, reduced to lowest terms. This requires the existence of a greatest common divisor function for the base ring.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/generic/Fraction.jl#L653' class='documenter-source'>source</a><br>


Here are some examples of computing a greatest common divisor.


```
S, x = PolynomialRing(ZZ, "x")

a = -x//(2x + 1)

f = gcd(a, (x + 1)//(x - 1))
```


<a id='Modular-arithmetic-1'></a>

## Modular arithmetic


The following functions are available for rationals.

<a id='Base.mod-Tuple{Nemo.fmpq,Nemo.fmpz}' href='#Base.mod-Tuple{Nemo.fmpq,Nemo.fmpz}'>#</a>
**`Base.mod`** &mdash; *Method*.



```
mod(a::fmpq, b::fmpz)
```

> Return $a \pmod{b}$ where $b$ is an integer coprime to the denominator of $a$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L420' class='documenter-source'>source</a><br>

<a id='Base.mod-Tuple{Nemo.fmpq,Integer}' href='#Base.mod-Tuple{Nemo.fmpq,Integer}'>#</a>
**`Base.mod`** &mdash; *Method*.



```
mod(a::fmpq, b::Integer)
```

> Return $a \pmod{b}$ where $b$ is an integer coprime to the denominator of $a$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L432' class='documenter-source'>source</a><br>


Here are some examples of modular arithmetic.


```
a = -fmpz(2)//3
b = fmpz(1)//2

c = mod(a, 7)
d = mod(b, fmpz(5))
```


<a id='Rational-Reconstruction-1'></a>

### Rational Reconstruction


Rational reconstruction is available for rational numbers.

<a id='Nemo.reconstruct-Tuple{Nemo.fmpz,Nemo.fmpz}' href='#Nemo.reconstruct-Tuple{Nemo.fmpz,Nemo.fmpz}'>#</a>
**`Nemo.reconstruct`** &mdash; *Method*.



```
reconstruct(a::fmpz, b::fmpz)
```

> Attempt to find a rational number $n/d$ such that  $0 \leq |n| \leq \lfloor\sqrt{m/2}\rfloor$ and  $0 < d \leq \lfloor\sqrt{m/2}\rfloor$ such that gcd$(n, d) = 1$ and $a \equiv nd^{-1} \pmod{m}$. If no solution exists, an exception is thrown.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L458' class='documenter-source'>source</a><br>

<a id='Nemo.reconstruct-Tuple{Nemo.fmpz,Integer}' href='#Nemo.reconstruct-Tuple{Nemo.fmpz,Integer}'>#</a>
**`Nemo.reconstruct`** &mdash; *Method*.



```
reconstruct(a::fmpz, b::Integer)
```

> Attempt to find a rational number $n/d$ such that  $0 \leq |n| \leq \lfloor\sqrt{m/2}\rfloor$ and  $0 < d \leq \lfloor\sqrt{m/2}\rfloor$ such that gcd$(n, d) = 1$ and $a \equiv nd^{-1} \pmod{m}$. If no solution exists, an exception is thrown.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L474' class='documenter-source'>source</a><br>

<a id='Nemo.reconstruct-Tuple{Integer,Nemo.fmpz}' href='#Nemo.reconstruct-Tuple{Integer,Nemo.fmpz}'>#</a>
**`Nemo.reconstruct`** &mdash; *Method*.



```
reconstruct(a::Integer, b::fmpz)
```

> Attempt to find a rational number $n/d$ such that  $0 \leq |n| \leq \lfloor\sqrt{m/2}\rfloor$ and  $0 < d \leq \lfloor\sqrt{m/2}\rfloor$ such that gcd$(n, d) = 1$ and $a \equiv nd^{-1} \pmod{m}$. If no solution exists, an exception is thrown.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L483' class='documenter-source'>source</a><br>

<a id='Nemo.reconstruct-Tuple{Integer,Integer}' href='#Nemo.reconstruct-Tuple{Integer,Integer}'>#</a>
**`Nemo.reconstruct`** &mdash; *Method*.



```
reconstruct(a::Integer, b::Integer)
```

> Attempt to find a rational number $n/d$ such that  $0 \leq |n| \leq \lfloor\sqrt{m/2}\rfloor$ and  $0 < d \leq \lfloor\sqrt{m/2}\rfloor$ such that gcd$(n, d) = 1$ and $a \equiv nd^{-1} \pmod{m}$. If no solution exists, an exception is thrown.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L492' class='documenter-source'>source</a><br>


Here are some examples of rational reconstruction.


```
a = reconstruct(7, 13)
b = reconstruct(fmpz(15), 31)
c = reconstruct(fmpz(123), fmpz(237))
```


<a id='Rational-enumeration-1'></a>

## Rational enumeration


Various methods exist to enumerator rationals.

<a id='Nemo.next_minimal-Tuple{Nemo.fmpq}' href='#Nemo.next_minimal-Tuple{Nemo.fmpq}'>#</a>
**`Nemo.next_minimal`** &mdash; *Method*.



```
next_minimal(a::fmpq)
```

> Given $x$, returns the next rational number in the sequence obtained by enumerating all positive denominators $q$, and for each $q$ enumerating the numerators $1 \le p < q$ in order and generating both $p/q$ and $q/p$, but skipping all gcd$(p,q) \neq 1$. Starting with zero, this generates every nonnegative rational number once and only once, with the first few entries being $0, 1, 1/2, 2, 1/3, 3, 2/3, 3/2, 1/4, 4, 3/4, 4/3, \ldots$. This enumeration produces the rational numbers in order of minimal height.  It has the disadvantage of being somewhat slower to compute than the Calkin-Wilf enumeration. If $x < 0$ we throw a `DomainError()`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L507' class='documenter-source'>source</a><br>

<a id='Nemo.next_signed_minimal-Tuple{Nemo.fmpq}' href='#Nemo.next_signed_minimal-Tuple{Nemo.fmpq}'>#</a>
**`Nemo.next_signed_minimal`** &mdash; *Method*.



```
next_signed_minimal(a::fmpq)
```

> Given a signed rational number $x$ assumed to be in canonical form,  returns the next element in the minimal-height sequence generated by  `next_minimal` but with negative numbers interleaved. The sequence begins $0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, \ldots$. Starting with zero, this generates every rational number once and only once, in order of minimal height.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L526' class='documenter-source'>source</a><br>

<a id='Nemo.next_calkin_wilf-Tuple{Nemo.fmpq}' href='#Nemo.next_calkin_wilf-Tuple{Nemo.fmpq}'>#</a>
**`Nemo.next_calkin_wilf`** &mdash; *Method*.



```
next_calkin_wilf(a::fmpq)
```

> Given $x$ return the next number in the breadth-first traversal of the Calkin-Wilf tree. Starting with zero, this generates every nonnegative rational number once and only once, with the first few entries being $0, 1, 1/2, 2, 1/3, 3/2, 2/3, 3, 1/4, 4/3, 3/5, 5/2, 2/5, \ldots$. Despite the appearance of the initial entries, the Calkin-Wilf enumeration  does not produce the rational numbers in order of height: some small fractions will appear late in the sequence. This order has the advantage of being faster to produce than the minimal-height order.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L542' class='documenter-source'>source</a><br>

<a id='Nemo.next_signed_calkin_wilf-Tuple{Nemo.fmpq}' href='#Nemo.next_signed_calkin_wilf-Tuple{Nemo.fmpq}'>#</a>
**`Nemo.next_signed_calkin_wilf`** &mdash; *Method*.



```
next_signed_calkin_wilf(a::fmpq)
```

> Given a signed rational number $x$ returns the next element in the Calkin-Wilf sequence with negative numbers interleaved. The sequence begins $0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, \ldots$. Starting with zero, this generates every rational number once and only once, but not in order of minimal height.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L561' class='documenter-source'>source</a><br>


Here are some examples of rational enumeration.


```
next_minimal(fmpz(2)//3)
next_signed_minimal(-fmpz(21)//31)
next_calkin_wilf(fmpz(321)//113)
next_signed_calkin_wilf(-fmpz(51)//(17))
```


<a id='Special-functions-1'></a>

## Special functions


The following special functions are available for specific rings in Nemo.

<a id='Nemo.harmonic-Tuple{Int64}' href='#Nemo.harmonic-Tuple{Int64}'>#</a>
**`Nemo.harmonic`** &mdash; *Method*.



```
harmonic(n::Int)
```

> Computes the harmonic number $H_n = 1 + 1/2 + 1/3 + \cdots + 1/n$. Table lookup is used for $H_n$ whose numerator and denominator  fit in a single limb. For larger $n$, a divide and conquer strategy is used.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L582' class='documenter-source'>source</a><br>

<a id='Nemo.bernoulli-Tuple{Int64}' href='#Nemo.bernoulli-Tuple{Int64}'>#</a>
**`Nemo.bernoulli`** &mdash; *Method*.



```
bernoulli(n::Int)
```

> Computes the Bernoulli number $B_n$ for nonnegative $n$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L595' class='documenter-source'>source</a><br>

<a id='Nemo.bernoulli_cache-Tuple{Int64}' href='#Nemo.bernoulli_cache-Tuple{Int64}'>#</a>
**`Nemo.bernoulli_cache`** &mdash; *Method*.



```
bernoulli_cache(n::Int)
```

> Precomputes and caches all the Bernoulli numbers up to $B_n$. This is much faster than repeatedly calling `bernoulli(k)`. Once cached, subsequent calls to `bernoulli(k)` for any $k \le n$ will read from the cache, making them virtually free.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L606' class='documenter-source'>source</a><br>

<a id='Nemo.dedekind_sum-Tuple{Nemo.fmpz,Nemo.fmpz}' href='#Nemo.dedekind_sum-Tuple{Nemo.fmpz,Nemo.fmpz}'>#</a>
**`Nemo.dedekind_sum`** &mdash; *Method*.



```
dedekind_sum(h::fmpz, k::fmpz)
```


<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L619' class='documenter-source'>source</a><br>

<a id='Nemo.dedekind_sum-Tuple{Nemo.fmpz,Integer}' href='#Nemo.dedekind_sum-Tuple{Nemo.fmpz,Integer}'>#</a>
**`Nemo.dedekind_sum`** &mdash; *Method*.



```
dedekind_sum(h::fmpz, k::Integer)
```

> Computes the Dedekind sum $s(h,k)$ for arbitrary $h$ and $k$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L629' class='documenter-source'>source</a><br>

<a id='Nemo.dedekind_sum-Tuple{Integer,Nemo.fmpz}' href='#Nemo.dedekind_sum-Tuple{Integer,Nemo.fmpz}'>#</a>
**`Nemo.dedekind_sum`** &mdash; *Method*.



```
dedekind_sum(h::Integer, k::fmpz)
```

> Computes the Dedekind sum $s(h,k)$ for arbitrary $h$ and $k$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L635' class='documenter-source'>source</a><br>

<a id='Nemo.dedekind_sum-Tuple{Integer,Integer}' href='#Nemo.dedekind_sum-Tuple{Integer,Integer}'>#</a>
**`Nemo.dedekind_sum`** &mdash; *Method*.



```
dedekind_sum(h::Integer, k::Integer)
```

> Computes the Dedekind sum $s(h,k)$ for arbitrary $h$ and $k$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/00727ca77a4ddfb3293c0b6590c674f002191822/src/flint/fmpq.jl#L641' class='documenter-source'>source</a><br>


Here are some examples of special functions.


```
a = harmonic(12)

b = dedekind_sum(12, 13)
c = dedekind_sum(-120, fmpz(1305))

d = bernoulli(12)

bernoulli_cache(100)
e = bernoulli(100)
```

