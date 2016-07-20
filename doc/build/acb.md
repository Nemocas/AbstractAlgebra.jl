


<a id='Introduction-1'></a>

## Introduction


Arbitrary precision complex ball arithmetic is supplied by Arb which provides a ball representation which tracks error bounds rigorously. Complex numbers are  represented in in rectangular form $a+bi$ where $a,b$ are `arb` balls.


The Arb complex field is constructed using the `AcbField` constructor. This constructs the parent object for the Arb complex field.


We define


```
ComplexField = AcbField
```


so that one can construct the Arb complex field parent using `ComplexField` instead of `AcbField`.


The types of complex boxes in Nemo are given in the following table, along with the libraries that provide them and the associated types of the parent objects.


| Library |                Field | Element type | Parent type |
| -------:| --------------------:| ------------:| -----------:|
|     Arb | $\mathbb{C}$ (boxes) |        `acb` |  `AcbField` |


All the complex field types belong to the `Field` abstract type and the types of elements in this field, i.e. complex boxes in this case, belong to the `FieldElem` abstract type.


<a id='Complex-field-constructors-1'></a>

## Complex field constructors


In order to construct complex boxes in Nemo, one must first construct the Arb complex field itself. This is accomplished with the following constructor.


```
AcbField(prec::Int)
```


Return the Arb complex field with precision in bits `prec` used for operations on interval midpoints. The precision used for interval radii is a fixed implementation-defined constant (30 bits).


We define


```
ComplexField = AcbField
```


so that one can use `ComplexField` instead of `AcbField`.


Here is an example of creating an Arb complex field and using the resulting parent object to coerce values into the resulting field.


```
CC = ComplexField(64)

a = CC("0.25")
b = CC("0.1")
c = CC(0.5)
d = CC(12)
```


Note that whilst one can coerce double precision floating point values into an Arb complex field, unless those values can be represented exactly in double precision, the resulting ball can't be any more precise than the double precision supplied.


If instead, values can be represented precisely using decimal arithmetic then one can supply them to Arb using a string. In this case, Arb will store them to the precision specified when creating the Arb complex field.


If the values can be stored precisely as a binary floating point number, Arb will store the values exactly. See the function `isexact` below for more information.


<a id='Complex-ball-constructors-1'></a>

## Complex ball constructors


Once an Arb complex field is constructed, there are various ways to construct boxes in that field.


Apart from coercing elements into the Arb complex field as above, we offer the following functions.

<a id='Base.zero-Tuple{Nemo.AcbField}' href='#Base.zero-Tuple{Nemo.AcbField}'>#</a>
**`Base.zero`** &mdash; *Method*.



```
zero(R::AcbField)
```

> Return exact zero in the given Arb complex field.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L58' class='documenter-source'>source</a><br>

<a id='Base.one-Tuple{Nemo.AcbField}' href='#Base.one-Tuple{Nemo.AcbField}'>#</a>
**`Base.one`** &mdash; *Method*.



```
one(R::AcbField)
```

> Return exact one in the given Arb complex field.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L68' class='documenter-source'>source</a><br>

<a id='Nemo.onei-Tuple{Nemo.AcbField}' href='#Nemo.onei-Tuple{Nemo.AcbField}'>#</a>
**`Nemo.onei`** &mdash; *Method*.



```
onei(R::AcbField)
```

> Return exact one times $i$ in the given Arb complex field.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L79' class='documenter-source'>source</a><br>


Here are some examples of constructing complex boxes.


```
CC = ComplexField(64)

a = zero(CC)
b = one(CC)
c = onei(CC)
```


<a id='Basic-functionality-1'></a>

## Basic functionality


The following basic functionality is provided by the default Arb complex field implementation in Nemo, to support construction of generic rings over complex fields. Any custom complex field implementation in Nemo should provide analogues of these functions along with the usual arithmetic operations.


```
parent_type(::Type{acb})
```


Gives the type of the parent object of an Arb complex field element.


```
elem_type(R::AcbField)
```


Given the parent object for an Arb complex field, return the type of elements of the field.


```
mul!(c::acb, a::acb, b::acb)
```


Multiply $a$ by $b$ and set the existing Arb complex field element $c$ to the result. This function is provided for performance reasons as it saves allocating a new object for the result and eliminates associated garbage collection.


```
addeq!(c::acb, a::acb)
```


In-place addition. Adds $a$ to $c$ and sets $c$ to the result. This function is provided for performance reasons as it saves allocating a new object for the result and eliminates associated garbage collection.


```
deepcopy(a::acb)
```


Return a copy of the Arb complex field element $a$, recursively copying the internal data. Arb complex field elements are mutable in Nemo so a shallow copy is not sufficient.


Given the parent object `R` for an Arb complex field, the following coercion functions are provided to coerce various elements into the Arb complex field. Developers provide these by overloading the `call` operator for the complex field parent objects.


```
R()
```


Coerce zero into the Arb complex field.


```
R(n::Integer)
R(f::fmpz)
R(q::fmpq)
```


Coerce an integer or rational value into the Arb complex field.


```
R(f::Float64)
R(f::BigFloat)
```


Coerce the given floating point number into the Arb complex field.


```
R(f::AbstractString)
R(f::AbstractString, g::AbstractString)
```


Coerce the decimal number, given as a string, into the Arb complex field. In each case $f$ is the real part and $g$ is the imaginary part.


```
R(f::arb)
```


Coerce the given Arb real ball into the Arb complex field.


```
R(f::acb)
```


Take an Arb complex field element that is already in an Arb field and simply return it. A copy of the original is not made.


Here are some examples of coercing elements into the Arb complex field.


```
RR = RealField(64)
CC = ComplexField(64)

a = CC(3)
b = CC(QQ(2,3))
c = CC("3 +/- 0.0001")
d = CC("-1.24e+12345")
f = CC("nan +/- inf")
g = CC(RR(3))
```


In addition to the above, developers of custom complex field types must ensure that they provide the equivalent of the function `base_ring(R::AcbField)` which should return `Union{}`. In addition to this they should ensure that each complex field element contains a field `parent` specifying the parent object of the complex field element, or at least supply the equivalent of the function `parent(a::acb)` to return the parent object of a complex field element.


<a id='Basic-manipulation-1'></a>

## Basic manipulation


Numerous functions are provided to manipulate Arb complex field elements. Also see the section on basic functionality above.

<a id='Nemo.base_ring-Tuple{Nemo.AcbField}' href='#Nemo.base_ring-Tuple{Nemo.AcbField}'>#</a>
**`Nemo.base_ring`** &mdash; *Method*.



```
base_ring(R::AcbField)
```

> Returns `Union{}` since an Arb complex field does not depend on any other ring.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L38' class='documenter-source'>source</a><br>

<a id='Nemo.base_ring-Tuple{Nemo.acb}' href='#Nemo.base_ring-Tuple{Nemo.acb}'>#</a>
**`Nemo.base_ring`** &mdash; *Method*.



```
base_ring(a::acb)
```

> Returns `Union{}` since an Arb complex field does not depend on any other ring.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L45' class='documenter-source'>source</a><br>

<a id='Base.parent-Tuple{Nemo.acb}' href='#Base.parent-Tuple{Nemo.acb}'>#</a>
**`Base.parent`** &mdash; *Method*.



```
parent(x::acb)
```

> Return the parent of the given Arb complex field element.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L52' class='documenter-source'>source</a><br>

<a id='Nemo.iszero-Tuple{Nemo.acb}' href='#Nemo.iszero-Tuple{Nemo.acb}'>#</a>
**`Nemo.iszero`** &mdash; *Method*.



```
iszero(x::acb)
```

> Return `true` if $x$ is certainly zero, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L442' class='documenter-source'>source</a><br>

<a id='Nemo.isone-Tuple{Nemo.acb}' href='#Nemo.isone-Tuple{Nemo.acb}'>#</a>
**`Nemo.isone`** &mdash; *Method*.



```
isone(x::acb)
```

> Return `true` if $x$ is certainly zero, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L450' class='documenter-source'>source</a><br>

<a id='Base.isfinite-Tuple{Nemo.acb}' href='#Base.isfinite-Tuple{Nemo.acb}'>#</a>
**`Base.isfinite`** &mdash; *Method*.



```
isfinite(x::acb)
```

> Return `true` if $x$ is finite, i.e. its real and imaginary parts have finite midpoint and radius, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L458' class='documenter-source'>source</a><br>

<a id='Nemo.isexact-Tuple{Nemo.acb}' href='#Nemo.isexact-Tuple{Nemo.acb}'>#</a>
**`Nemo.isexact`** &mdash; *Method*.



```
isexact(x::acb)
```

> Return `true` if $x$ is exact, i.e. has its real and imaginary parts have zero radius, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L467' class='documenter-source'>source</a><br>

<a id='Nemo.isint-Tuple{Nemo.acb}' href='#Nemo.isint-Tuple{Nemo.acb}'>#</a>
**`Nemo.isint`** &mdash; *Method*.



```
isint(x::acb)
```

> Return `true` if $x$ is an exact integer, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L476' class='documenter-source'>source</a><br>

<a id='Base.isreal-Tuple{Nemo.acb}' href='#Base.isreal-Tuple{Nemo.acb}'>#</a>
**`Base.isreal`** &mdash; *Method*.



```
isint(x::acb)
```

> Return `true` if $x$ is purely real, i.e. having zero imaginary part, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L484' class='documenter-source'>source</a><br>

<a id='Base.real-Tuple{Nemo.acb}' href='#Base.real-Tuple{Nemo.acb}'>#</a>
**`Base.real`** &mdash; *Method*.



```
real(x::acb)
```

> Return the real part of $x$ as an `arb`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L131' class='documenter-source'>source</a><br>

<a id='Base.imag-Tuple{Nemo.acb}' href='#Base.imag-Tuple{Nemo.acb}'>#</a>
**`Base.imag`** &mdash; *Method*.



```
imag(x::acb)
```

> Return the imaginary part of $x$ as an `arb`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L142' class='documenter-source'>source</a><br>

<a id='Nemo.accuracy_bits-Tuple{Nemo.acb}' href='#Nemo.accuracy_bits-Tuple{Nemo.acb}'>#</a>
**`Nemo.accuracy_bits`** &mdash; *Method*.



```
accuracy_bits(x::acb)
```

> Return the relative accuracy of $x$ measured in bits, capped between `typemax(Int)` and `-typemax(Int)`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L90' class='documenter-source'>source</a><br>


Here are some examples of basic manipulation of Arb complex boxes.


```
CC = ComplexField(64)

a = CC("1.2 +/- 0.001")
b = CC(3)

iszero(a)
isone(b)
isreal(a)
isfinite(b)
isint(b)
c = real(a)
d = imag(b)
f = accuracy_bits(a)
S = parent(a)
T = base_ring(CC)
```


<a id='Arithmetic-operations-1'></a>

## Arithmetic operations


Nemo provides all the standard field operations for Arb complex field elements, as follows. Note that division is represented by `//` since a field is its own fraction field and since exact division is not generally possible in an inexact field.


|           Function |      Operation |
| ------------------:| --------------:|
|          -(a::acb) |    unary minus |
|  +(a::acb, b::acb) |       addition |
|  -(a::acb, b::acb) |    subtraction |
|  *(a::acb, b::acb) | multiplication |
| //(a::acb, b::acb) |       division |
|  ^(a::acb, b::acb) |       powering |


In addition, the following ad hoc field operations are defined.


|               Function |      Operation |
| ----------------------:| --------------:|
|  +(a::acb, b::Integer) |       addition |
|  +(a::Integer, b::acb) |       addition |
|     +(a::acb, b::fmpz) |       addition |
|     +(a::fmpz, b::acb) |       addition |
|     +(a::acb, b::fmpq) |       addition |
|     +(a::fmpq, b::acb) |       addition |
|      +(a::acb, b::arb) |       addition |
|      +(a::arb, b::acb) |       addition |
|  -(a::acb, b::Integer) |    subtraction |
|  -(a::Integer, b::acb) |    subtraction |
|     -(a::acb, b::fmpz) |    subtraction |
|     -(a::fmpz, b::acb) |    subtraction |
|     -(a::acb, b::fmpq) |    subtraction |
|     -(a::fmpq, b::acb) |    subtraction |
|      -(a::acb, b::arb) |    subtraction |
|      -(a::arb, b::acb) |    subtraction |
|  *(a::acb, b::Integer) | multiplication |
|  *(a::Integer, b::acb) | multiplication |
|     *(a::acb, b::fmpz) | multiplication |
|     *(a::fmpz, b::acb) | multiplication |
|     *(a::acb, b::fmpq) | multiplication |
|     *(a::fmpq, b::acb) | multiplication |
|      *(a::acb, b::arb) | multiplication |
|      *(a::arb, b::acb) | multiplication |
| //(a::acb, b::Integer) |       division |
|    //(a::acb, b::fmpz) |       division |
|    //(a::acb, b::fmpq) |       division |
| //(a::Integer, b::acb) |       division |
|    //(a::fmpz, b::acb) |       division |
|    //(a::fmpq, b::acb) |       division |
|     //(a::arb, b::acb) |       division |
|     ^(a::acb, b::fmpq) |       powering |
|  ^(a::Integer, b::acb) |       powering |
|     ^(a::fmpz, b::acb) |       powering |
|     ^(a::fmpq, b::acb) |       powering |
|      ^(a::arb, b::acb) |       powering |


Here are some examples of arithmetic operations on Arb complex boxes.


```
CC = ComplexField(64)

x = CC(3)
y = CC(QQ(2,3))

a = x + y
b = x*y
d = 1//y - x
d = 3x + ZZ(100)//y
f = (x^2 + y^2) ^ QQ(1,2)
```


<a id='Containment-1'></a>

## Containment


It is often necessary to determine whether a given exact value or box is contained in a given complex box or whether two boxes overlap. The following functions are provided for this purpose.

<a id='Nemo.overlaps-Tuple{Nemo.acb,Nemo.acb}' href='#Nemo.overlaps-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Nemo.overlaps`** &mdash; *Method*.



```
overlaps(x::acb, y::acb)
```

> Returns `true` if any part of the box $x$ overlaps any part of the box $y$, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L375' class='documenter-source'>source</a><br>

<a id='Base.contains-Tuple{Nemo.acb,Nemo.acb}' href='#Base.contains-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Base.contains`** &mdash; *Method*.



```
contains(x::acb, y::acb)
```

> Returns `true` if the box $x$ contains the box $y$, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L385' class='documenter-source'>source</a><br>

<a id='Base.contains-Tuple{Nemo.acb,Integer}' href='#Base.contains-Tuple{Nemo.acb,Integer}'>#</a>
**`Base.contains`** &mdash; *Method*.



```
contains(x::acb, y::Integer)
```

> Returns `true` if the box $x$ contains the given integer value, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L421' class='documenter-source'>source</a><br>

<a id='Base.contains-Tuple{Nemo.acb,Nemo.fmpz}' href='#Base.contains-Tuple{Nemo.acb,Nemo.fmpz}'>#</a>
**`Base.contains`** &mdash; *Method*.



```
contains(x::acb, y::fmpz)
```

> Returns `true` if the box $x$ contains the given integer value, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L405' class='documenter-source'>source</a><br>

<a id='Base.contains-Tuple{Nemo.acb,Nemo.fmpq}' href='#Base.contains-Tuple{Nemo.acb,Nemo.fmpq}'>#</a>
**`Base.contains`** &mdash; *Method*.



```
contains(x::acb, y::fmpq)
```

> Returns `true` if the box $x$ contains the given rational value, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L395' class='documenter-source'>source</a><br>


The following functions are also provided for determining if a box intersects a certain part of the complex number plane.

<a id='Nemo.contains_zero-Tuple{Nemo.acb}' href='#Nemo.contains_zero-Tuple{Nemo.acb}'>#</a>
**`Nemo.contains_zero`** &mdash; *Method*.



```
contains_zero(x::acb)
```

> Returns `true` if the box $x$ contains zero, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L428' class='documenter-source'>source</a><br>


Here are some examples of testing containment.


```
CC = ComplexField(64)
x = CC("1 +/- 0.001")
y = CC("3")

overlaps(x, y)
contains(x, y)
contains(y, 3)
contains(x, ZZ(1)//2)
contains_zero(x)
```


<a id='Comparison-1'></a>

## Comparison


Nemo provides a full range of comparison operations for Arb complex boxes. 


Firstly, we introduce an exact equality which is distinct from arithmetic equality. This is distinct from arithmetic equality implemented by `==`, which merely compares up to the minimum of the precisions of its operands.

<a id='Base.isequal-Tuple{Nemo.acb,Nemo.acb}' href='#Base.isequal-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Base.isequal`** &mdash; *Method*.



```
isequal(x::acb, y::acb)
```

> Return `true` if the boxes $x$ and $y$ are precisely equal, i.e. their real and imaginary parts have the same midpoints and radii.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L322' class='documenter-source'>source</a><br>


A full range of functions is available for comparing boxes, i.e. `==` and `!=`. In fact, these are implemented directly in C. In the table below we document these as though only `==` had been provided to Julia.


<a id='Function-1'></a>

## Function


`==(x::acb, y::acb)`


As well as these, we provide a full range of ad hoc comparison operators. Again, these are implemented directly in Julia, but we document them as though only `==` were provided.


<a id='Function-2'></a>

## Function


`==(x::acb, y::Integer)` `==(x::Integer, y::acb)` `==(x::acb, y::fmpz)` `==(x::fmpz, y::acb)` `==(x::arb, y::fmpz)` `==(x::fmpz, y::arb)` `==(x::acb, y::Float64)` `==(x::Float64, y::acb)`


Here are some examples of comparison.


```
CC = ComplexField(64)
x = CC("1 +/- 0.001")
y = CC("3")
z = CC("4")

isequal(y, z)
x != z
x == 3
ZZ(3) == z
x != 1.23
```


<a id='Absolute-value-1'></a>

## Absolute value

<a id='Base.abs-Tuple{Nemo.acb}' href='#Base.abs-Tuple{Nemo.acb}'>#</a>
**`Base.abs`** &mdash; *Method*.



```
abs(x::acb)
```

> Return the complex absolute value of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L499' class='documenter-source'>source</a><br>


Here are some examples of taking the absolute value.


```
CC = ComplexField(64)
x = CC("-1 +/- 0.001")

a = abs(x)
```


<a id='Inverse-1'></a>

## Inverse

<a id='Base.inv-Tuple{Nemo.acb}' href='#Base.inv-Tuple{Nemo.acb}'>#</a>
**`Base.inv`** &mdash; *Method*.



```
inv(x::acb)
```

> Return the multiplicative inverse of $x$, i.e. $1/x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L517' class='documenter-source'>source</a><br>


Here are some examples of taking the inverse.


```
CC = ComplexField(64)
x = CC("-3 +/- 0.001")

a = inv(x)
```


<a id='Shifting-1'></a>

## Shifting

<a id='Base.Math.ldexp-Tuple{Nemo.acb,Int64}' href='#Base.Math.ldexp-Tuple{Nemo.acb,Int64}'>#</a>
**`Base.Math.ldexp`** &mdash; *Method*.



```
ldexp(x::acb, y::Int)
```

> Return $2^yx$. Note that $y$ can be positive, zero or negative.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L533' class='documenter-source'>source</a><br>

<a id='Base.Math.ldexp-Tuple{Nemo.acb,Nemo.fmpz}' href='#Base.Math.ldexp-Tuple{Nemo.acb,Nemo.fmpz}'>#</a>
**`Base.Math.ldexp`** &mdash; *Method*.



```
ldexp(x::acb, y::fmpz)
```

> Return $2^yx$. Note that $y$ can be positive, zero or negative.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L544' class='documenter-source'>source</a><br>


Here are some examples of shifting.


```
CC = ComplexField(64)
x = CC("-3 +/- 0.001")

a = ldexp(x, 23)
b = ldexp(x, -ZZ(15))
```


<a id='Miscellaneous-operations-1'></a>

## Miscellaneous operations

<a id='Nemo.trim-Tuple{Nemo.acb}' href='#Nemo.trim-Tuple{Nemo.acb}'>#</a>
**`Nemo.trim`** &mdash; *Method*.



```
trim(x::acb)
```

> Return an `acb` box containing $x$ but which may be more economical, by rounding off insignificant bits from midpoints.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L561' class='documenter-source'>source</a><br>

<a id='Nemo.unique_integer-Tuple{Nemo.acb}' href='#Nemo.unique_integer-Tuple{Nemo.acb}'>#</a>
**`Nemo.unique_integer`** &mdash; *Method*.



```
unique_integer(x::acb)
```

> Return a pair where the first value is a boolean and the second is an `fmpz` integer. The boolean indicates whether the box $x$ contains a unique integer. If this is the case, the second return value is set to this unique integer.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L572' class='documenter-source'>source</a><br>

<a id='Base.conj-Tuple{Nemo.acb}' href='#Base.conj-Tuple{Nemo.acb}'>#</a>
**`Base.conj`** &mdash; *Method*.



```
conj(x::acb)
```

> Return the complex conjugate of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L586' class='documenter-source'>source</a><br>

<a id='Base.angle-Tuple{Nemo.acb}' href='#Base.angle-Tuple{Nemo.acb}'>#</a>
**`Base.angle`** &mdash; *Method*.



```
angle(x::acb)
```

> Return the angle in radians that the complex vector $x$ makes with the positive real axis in a counterclockwise direction.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L596' class='documenter-source'>source</a><br>


Here are some examples of miscellaneous operations.


```
CC = ComplexField(64)
x = CC("-3 +/- 0.001", "0.1")

a = trim(x)
b, c = unique_integer(x)
d = conj(x)
f = angle(x)
```


<a id='Constants-1'></a>

## Constants

<a id='Nemo.const_pi-Tuple{Nemo.AcbField}' href='#Nemo.const_pi-Tuple{Nemo.AcbField}'>#</a>
**`Nemo.const_pi`** &mdash; *Method*.



Here are some examples of computing complex constants.


```
CC = ComplexField(200)

a = const_pi(CC)
```


<a id='Mathematical-functions-1'></a>

## Mathematical functions

<a id='Base.sqrt-Tuple{Nemo.acb}' href='#Base.sqrt-Tuple{Nemo.acb}'>#</a>
**`Base.sqrt`** &mdash; *Method*.



```
sqrt(x::acb)
```

> Return the square root of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L629' class='documenter-source'>source</a><br>

<a id='Nemo.rsqrt-Tuple{Nemo.acb}' href='#Nemo.rsqrt-Tuple{Nemo.acb}'>#</a>
**`Nemo.rsqrt`** &mdash; *Method*.



```
rsqrt(x::acb)
```

> Return the reciprocal of the square root of $x$, i.e. $1/\sqrt{x}$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L639' class='documenter-source'>source</a><br>

<a id='Base.log-Tuple{Nemo.acb}' href='#Base.log-Tuple{Nemo.acb}'>#</a>
**`Base.log`** &mdash; *Method*.



```
log(x)
```

Compute the natural logarithm of `x`. Throws `DomainError` for negative `Real` arguments. Use complex negative arguments to obtain complex results.

There is an experimental variant in the `Base.Math.JuliaLibm` module, which is typically faster and more accurate.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L3536-3544' class='documenter-source'>source</a><br>


```
log(x::acb)
```

> Return the principal branch of the logarithm of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L649' class='documenter-source'>source</a><br>

<a id='Base.log1p-Tuple{Nemo.acb}' href='#Base.log1p-Tuple{Nemo.acb}'>#</a>
**`Base.log1p`** &mdash; *Method*.



```
log1p(x::acb)
```

> Return $\log(1+x)$, evaluated accurately for small $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L659' class='documenter-source'>source</a><br>

<a id='Base.exp-Tuple{Nemo.acb}' href='#Base.exp-Tuple{Nemo.acb}'>#</a>
**`Base.exp`** &mdash; *Method*.



```
exp(x::acb)
```

> Return the exponential of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L669' class='documenter-source'>source</a><br>

<a id='Nemo.exppii-Tuple{Nemo.acb}' href='#Nemo.exppii-Tuple{Nemo.acb}'>#</a>
**`Nemo.exppii`** &mdash; *Method*.



```
exppii(x::acb)
```

> Return the exponential of $\pi i x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L679' class='documenter-source'>source</a><br>

<a id='Base.sin-Tuple{Nemo.acb}' href='#Base.sin-Tuple{Nemo.acb}'>#</a>
**`Base.sin`** &mdash; *Method*.



```
sin(x::acb)
```

> Return the sine of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L689' class='documenter-source'>source</a><br>

<a id='Base.cos-Tuple{Nemo.acb}' href='#Base.cos-Tuple{Nemo.acb}'>#</a>
**`Base.cos`** &mdash; *Method*.



```
cos(x::acb)
```

> Return the cosine of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L699' class='documenter-source'>source</a><br>

<a id='Base.Math.sinpi-Tuple{Nemo.acb}' href='#Base.Math.sinpi-Tuple{Nemo.acb}'>#</a>
**`Base.Math.sinpi`** &mdash; *Method*.



```
sinpi(x::acb)
```

> Return the sine of $\pi x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L729' class='documenter-source'>source</a><br>

<a id='Base.Math.cospi-Tuple{Nemo.acb}' href='#Base.Math.cospi-Tuple{Nemo.acb}'>#</a>
**`Base.Math.cospi`** &mdash; *Method*.



```
cospi(x::acb)
```

> Return the cosine of $\pi x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L739' class='documenter-source'>source</a><br>

<a id='Base.tan-Tuple{Nemo.acb}' href='#Base.tan-Tuple{Nemo.acb}'>#</a>
**`Base.tan`** &mdash; *Method*.



```
tan(x::acb)
```

> Return the tangent of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L709' class='documenter-source'>source</a><br>

<a id='Base.Math.cot-Tuple{Nemo.acb}' href='#Base.Math.cot-Tuple{Nemo.acb}'>#</a>
**`Base.Math.cot`** &mdash; *Method*.



```
cot(x::acb)
```

> Return the cotangent of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L719' class='documenter-source'>source</a><br>

<a id='Nemo.tanpi-Tuple{Nemo.acb}' href='#Nemo.tanpi-Tuple{Nemo.acb}'>#</a>
**`Nemo.tanpi`** &mdash; *Method*.



```
tanpi(x::acb)
```

> Return the tangent of $\pi x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L749' class='documenter-source'>source</a><br>

<a id='Nemo.cotpi-Tuple{Nemo.acb}' href='#Nemo.cotpi-Tuple{Nemo.acb}'>#</a>
**`Nemo.cotpi`** &mdash; *Method*.



```
cotpi(x::acb)
```

> Return the cotangent of $\pi x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L759' class='documenter-source'>source</a><br>

<a id='Base.sinh-Tuple{Nemo.acb}' href='#Base.sinh-Tuple{Nemo.acb}'>#</a>
**`Base.sinh`** &mdash; *Method*.



```
sinh(x::acb)
```

> Return the hyperbolic sine of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L769' class='documenter-source'>source</a><br>

<a id='Base.cosh-Tuple{Nemo.acb}' href='#Base.cosh-Tuple{Nemo.acb}'>#</a>
**`Base.cosh`** &mdash; *Method*.



```
cosh(x::acb)
```

> Return the hyperbolic cosine of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L779' class='documenter-source'>source</a><br>

<a id='Base.tanh-Tuple{Nemo.acb}' href='#Base.tanh-Tuple{Nemo.acb}'>#</a>
**`Base.tanh`** &mdash; *Method*.



```
tanh(x::acb)
```

> Return the hyperbolic tangent of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L789' class='documenter-source'>source</a><br>

<a id='Base.Math.coth-Tuple{Nemo.acb}' href='#Base.Math.coth-Tuple{Nemo.acb}'>#</a>
**`Base.Math.coth`** &mdash; *Method*.



```
coth(x::acb)
```

> Return the hyperbolic cotangent of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L799' class='documenter-source'>source</a><br>

<a id='Base.atan-Tuple{Nemo.acb}' href='#Base.atan-Tuple{Nemo.acb}'>#</a>
**`Base.atan`** &mdash; *Method*.



```
atan(x::acb)
```

> Return the arctangent of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L809' class='documenter-source'>source</a><br>

<a id='Nemo.logsinpi-Tuple{Nemo.acb}' href='#Nemo.logsinpi-Tuple{Nemo.acb}'>#</a>
**`Nemo.logsinpi`** &mdash; *Method*.



```
logsinpi(x::acb)
```

> Return $\log\sin(\pi x)$, constructed without branch cuts off the real line.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L819' class='documenter-source'>source</a><br>

<a id='Base.Math.gamma-Tuple{Nemo.acb}' href='#Base.Math.gamma-Tuple{Nemo.acb}'>#</a>
**`Base.Math.gamma`** &mdash; *Method*.



```
gamma(x::acb)
```

> Return the Gamma function evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L829' class='documenter-source'>source</a><br>

<a id='Base.Math.lgamma-Tuple{Nemo.acb}' href='#Base.Math.lgamma-Tuple{Nemo.acb}'>#</a>
**`Base.Math.lgamma`** &mdash; *Method*.



```
lgamma(x::acb)
```

> Return the logarithm of the Gamma function evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L849' class='documenter-source'>source</a><br>

<a id='Nemo.rgamma-Tuple{Nemo.acb}' href='#Nemo.rgamma-Tuple{Nemo.acb}'>#</a>
**`Nemo.rgamma`** &mdash; *Method*.



```
rgamma(x::acb)
```

> Return the reciprocal of the Gamma function evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L839' class='documenter-source'>source</a><br>

<a id='Base.Math.digamma-Tuple{Nemo.acb}' href='#Base.Math.digamma-Tuple{Nemo.acb}'>#</a>
**`Base.Math.digamma`** &mdash; *Method*.



```
digamma(x::acb)
```

> Return the  logarithmic derivative of the gamma function evaluated at $x$, i.e. $\psi(x)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L859' class='documenter-source'>source</a><br>

<a id='Base.Math.zeta-Tuple{Nemo.acb}' href='#Base.Math.zeta-Tuple{Nemo.acb}'>#</a>
**`Base.Math.zeta`** &mdash; *Method*.



```
zeta(s)
```

Riemann zeta function $\zeta(s)$.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L9507-9511' class='documenter-source'>source</a><br>


```
zeta(x::acb)
```

> Return the Riemann zeta function evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L870' class='documenter-source'>source</a><br>

<a id='Nemo.barnesg-Tuple{Nemo.acb}' href='#Nemo.barnesg-Tuple{Nemo.acb}'>#</a>
**`Nemo.barnesg`** &mdash; *Method*.



```
barnesg(x::acb)
```

> Return the Barnes $G$-function, evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L880' class='documenter-source'>source</a><br>

<a id='Nemo.logbarnesg-Tuple{Nemo.acb}' href='#Nemo.logbarnesg-Tuple{Nemo.acb}'>#</a>
**`Nemo.logbarnesg`** &mdash; *Method*.



```
logbarnesg(x::acb)
```

> Return the logarithm of the Barnes $G$-function, evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L890' class='documenter-source'>source</a><br>

<a id='Base.Math.erf-Tuple{Nemo.acb}' href='#Base.Math.erf-Tuple{Nemo.acb}'>#</a>
**`Base.Math.erf`** &mdash; *Method*.



```
erf(x::acb)
```

> Return the error function evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L910' class='documenter-source'>source</a><br>

<a id='Base.Math.erfi-Tuple{Nemo.acb}' href='#Base.Math.erfi-Tuple{Nemo.acb}'>#</a>
**`Base.Math.erfi`** &mdash; *Method*.



```
erfi(x::acb)
```

> Return the imaginary error function evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L920' class='documenter-source'>source</a><br>

<a id='Nemo.ei-Tuple{Nemo.acb}' href='#Nemo.ei-Tuple{Nemo.acb}'>#</a>
**`Nemo.ei`** &mdash; *Method*.



```
ei(x::acb)
```

> Return the exponential integral evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L940' class='documenter-source'>source</a><br>

<a id='Nemo.si-Tuple{Nemo.acb}' href='#Nemo.si-Tuple{Nemo.acb}'>#</a>
**`Nemo.si`** &mdash; *Method*.



```
si(x::acb)
```

> Return the sine integral evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L950' class='documenter-source'>source</a><br>

<a id='Nemo.ci-Tuple{Nemo.acb}' href='#Nemo.ci-Tuple{Nemo.acb}'>#</a>
**`Nemo.ci`** &mdash; *Method*.



```
ci(x::acb)
```

> Return the exponential cosine integral evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L960' class='documenter-source'>source</a><br>

<a id='Nemo.shi-Tuple{Nemo.acb}' href='#Nemo.shi-Tuple{Nemo.acb}'>#</a>
**`Nemo.shi`** &mdash; *Method*.



```
shi(x::acb)
```

> Return the hyperbolic sine integral evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L970' class='documenter-source'>source</a><br>

<a id='Nemo.chi-Tuple{Nemo.acb}' href='#Nemo.chi-Tuple{Nemo.acb}'>#</a>
**`Nemo.chi`** &mdash; *Method*.



```
chi(x::acb)
```

> Return the hyperbolic cosine integral evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L980' class='documenter-source'>source</a><br>

<a id='Nemo.modeta-Tuple{Nemo.acb}' href='#Nemo.modeta-Tuple{Nemo.acb}'>#</a>
**`Nemo.modeta`** &mdash; *Method*.



```
modeta(x::acb)
```

> Return the Dedekind eta function $\eta(\tau)$ at $\tau = x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L990' class='documenter-source'>source</a><br>

<a id='Nemo.modj-Tuple{Nemo.acb}' href='#Nemo.modj-Tuple{Nemo.acb}'>#</a>
**`Nemo.modj`** &mdash; *Method*.



```
modj(x::acb)
```

> Return the $j$-invariant $j(\tau)$ at $\tau = x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1000' class='documenter-source'>source</a><br>

<a id='Nemo.modlambda-Tuple{Nemo.acb}' href='#Nemo.modlambda-Tuple{Nemo.acb}'>#</a>
**`Nemo.modlambda`** &mdash; *Method*.



```
modlambda(x::acb)
```

> Return the modular lambda function $\lambda(\tau)$ at $\tau = x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1010' class='documenter-source'>source</a><br>

<a id='Nemo.moddelta-Tuple{Nemo.acb}' href='#Nemo.moddelta-Tuple{Nemo.acb}'>#</a>
**`Nemo.moddelta`** &mdash; *Method*.



```
moddelta(x::acb)
```

> Return the modular delta function $\Delta(\tau)$ at $\tau = x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1020' class='documenter-source'>source</a><br>

<a id='Nemo.ellipk-Tuple{Nemo.acb}' href='#Nemo.ellipk-Tuple{Nemo.acb}'>#</a>
**`Nemo.ellipk`** &mdash; *Method*.



```
ellipk(x::acb)
```

> Return the complete elliptic integral $K(x)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1030' class='documenter-source'>source</a><br>

<a id='Nemo.ellipe-Tuple{Nemo.acb}' href='#Nemo.ellipe-Tuple{Nemo.acb}'>#</a>
**`Nemo.ellipe`** &mdash; *Method*.



```
ellipe(x::acb)
```

> Return the complete elliptic integral $E(x)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1040' class='documenter-source'>source</a><br>

<a id='Nemo.sincos-Tuple{Nemo.acb}' href='#Nemo.sincos-Tuple{Nemo.acb}'>#</a>
**`Nemo.sincos`** &mdash; *Method*.



```
sincos(x::acb)
```

> Return a tuple $s, c$ consisting of the sine $s$ and cosine $c$ of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1050' class='documenter-source'>source</a><br>

<a id='Nemo.sincospi-Tuple{Nemo.acb}' href='#Nemo.sincospi-Tuple{Nemo.acb}'>#</a>
**`Nemo.sincospi`** &mdash; *Method*.



```
sincospi(x::acb)
```

> Return a tuple $s, c$ consisting of the sine $s$ and cosine $c$ of $\pi x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1062' class='documenter-source'>source</a><br>

<a id='Nemo.sinhcosh-Tuple{Nemo.acb}' href='#Nemo.sinhcosh-Tuple{Nemo.acb}'>#</a>
**`Nemo.sinhcosh`** &mdash; *Method*.



```
sinhcosh(x::acb)
```

> Return a tuple $s, c$ consisting of the hyperbolic sine and cosine of $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1074' class='documenter-source'>source</a><br>

<a id='Nemo.agm-Tuple{Nemo.acb}' href='#Nemo.agm-Tuple{Nemo.acb}'>#</a>
**`Nemo.agm`** &mdash; *Method*.



```
agm(x::acb)
```

> Return the arithmetic-geometric mean of $1$ and $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L900' class='documenter-source'>source</a><br>

<a id='Nemo.agm-Tuple{Nemo.acb,Nemo.acb}' href='#Nemo.agm-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Nemo.agm`** &mdash; *Method*.



```
agm(x::acb, y::acb)
```

> Return the arithmetic-geometric mean of $x$ and $y$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1307' class='documenter-source'>source</a><br>

<a id='Base.Math.polygamma-Tuple{Nemo.acb,Nemo.acb}' href='#Base.Math.polygamma-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Base.Math.polygamma`** &mdash; *Method*.



```
polygamma(s::acb, a::acb)
```

> Return the generalised polygamma function $\psi(s,z)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1097' class='documenter-source'>source</a><br>

<a id='Base.Math.zeta-Tuple{Nemo.acb,Nemo.acb}' href='#Base.Math.zeta-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Base.Math.zeta`** &mdash; *Method*.



```
zeta(s, z)
```

Generalized zeta function $\zeta(s, z)$, defined by the sum $\sum_{k=0}^\infty ((k+z)^2)^{-s/2}$, where any term with $k+z=0$ is excluded.  For $\Re z > 0$, this definition is equivalent to the Hurwitz zeta function $\sum_{k=0}^\infty (k+z)^{-s}$.   For $z=1$, it yields the Riemann zeta function $\zeta(s)$.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/special/gamma.jl#L251-260' class='documenter-source'>source</a><br>


```
zeta(s::acb, a::acb)
```

> Return the Hurwitz zeta function $\zeta(s,a)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1086' class='documenter-source'>source</a><br>

<a id='Nemo.risingfac-Tuple{Nemo.acb,Int64}' href='#Nemo.risingfac-Tuple{Nemo.acb,Int64}'>#</a>
**`Nemo.risingfac`** &mdash; *Method*.



```
risingfac(x::acb, n::Int)
```

> Return the rising factorial $x(x + 1)\ldots (x + n - 1)$ as an Acb.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1115' class='documenter-source'>source</a><br>

<a id='Nemo.risingfac2-Tuple{Nemo.acb,Int64}' href='#Nemo.risingfac2-Tuple{Nemo.acb,Int64}'>#</a>
**`Nemo.risingfac2`** &mdash; *Method*.



```
risingfac2(x::acb, n::Int)
```

> Return a tuple containing the rising factorial $x(x + 1)\ldots (x + n - 1)$ and its derivative.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1129' class='documenter-source'>source</a><br>

<a id='Nemo.polylog-Tuple{Nemo.acb,Nemo.acb}' href='#Nemo.polylog-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Nemo.polylog`** &mdash; *Method*.



```
polylog(s::acb, a::acb)
```

>



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1136' class='documenter-source'>source</a><br>

<a id='Nemo.polylog-Tuple{Int64,Nemo.acb}' href='#Nemo.polylog-Tuple{Int64,Nemo.acb}'>#</a>
**`Nemo.polylog`** &mdash; *Method*.



```
polylog(s::Int, a::acb)
```

> Return the polylogarithm Li$_s(a)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1147' class='documenter-source'>source</a><br>

<a id='Nemo.li-Tuple{Nemo.acb}' href='#Nemo.li-Tuple{Nemo.acb}'>#</a>
**`Nemo.li`** &mdash; *Method*.



```
li(x::acb)
```

> Return the logarithmic integral, evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1158' class='documenter-source'>source</a><br>

<a id='Nemo.lioffset-Tuple{Nemo.acb}' href='#Nemo.lioffset-Tuple{Nemo.acb}'>#</a>
**`Nemo.lioffset`** &mdash; *Method*.



```
lioffset(x::acb)
```

> Return the offset logarithmic integral, evaluated at $x$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1169' class='documenter-source'>source</a><br>

<a id='Nemo.expint-Tuple{Nemo.acb}' href='#Nemo.expint-Tuple{Nemo.acb}'>#</a>
**`Nemo.expint`** &mdash; *Method*.


<a id='Base.Math.gamma-Tuple{Nemo.acb,Nemo.acb}' href='#Base.Math.gamma-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Base.Math.gamma`** &mdash; *Method*.



```
gamma(s::acb, x::acb)
```

> Return the upper incomplete gamma function $\Gamma(s,x)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1191' class='documenter-source'>source</a><br>

<a id='Base.Math.besselj-Tuple{Nemo.acb,Nemo.acb}' href='#Base.Math.besselj-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Base.Math.besselj`** &mdash; *Method*.



```
besselj(nu::acb, x::acb)
```

> Return the Bessel function $J_{\nu}(x)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1202' class='documenter-source'>source</a><br>

<a id='Base.Math.bessely-Tuple{Nemo.acb,Nemo.acb}' href='#Base.Math.bessely-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Base.Math.bessely`** &mdash; *Method*.



```
bessely(nu::acb, x::acb)
```

> Return the Bessel function $Y_{\nu}(x)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1213' class='documenter-source'>source</a><br>

<a id='Base.Math.besseli-Tuple{Nemo.acb,Nemo.acb}' href='#Base.Math.besseli-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Base.Math.besseli`** &mdash; *Method*.



```
besseli(nu::acb, x::acb)
```

> Return the Bessel function $I_{\nu}(x)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1224' class='documenter-source'>source</a><br>

<a id='Base.Math.besselk-Tuple{Nemo.acb,Nemo.acb}' href='#Base.Math.besselk-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Base.Math.besselk`** &mdash; *Method*.



```
besselk(nu::acb, x::acb)
```

> Return the Bessel function $K_{\nu}(x)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1235' class='documenter-source'>source</a><br>

<a id='Nemo.hyp1f1-Tuple{Nemo.acb,Nemo.acb,Nemo.acb}' href='#Nemo.hyp1f1-Tuple{Nemo.acb,Nemo.acb,Nemo.acb}'>#</a>
**`Nemo.hyp1f1`** &mdash; *Method*.



```
hyp1f1(a::acb, b::acb, x::acb)
```

> Return the confluent hypergeometric function ${}_1F1(a,b,x)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1246' class='documenter-source'>source</a><br>

<a id='Nemo.hyp1f1r-Tuple{Nemo.acb,Nemo.acb,Nemo.acb}' href='#Nemo.hyp1f1r-Tuple{Nemo.acb,Nemo.acb,Nemo.acb}'>#</a>
**`Nemo.hyp1f1r`** &mdash; *Method*.



```
hyp1f1r(a::acb, b::acb, x::acb)
```

> Return the regularized confluent hypergeometric function ${}_1F1(a,b,x) / \Gamma(b)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1257' class='documenter-source'>source</a><br>

<a id='Nemo.hyperu-Tuple{Nemo.acb,Nemo.acb,Nemo.acb}' href='#Nemo.hyperu-Tuple{Nemo.acb,Nemo.acb,Nemo.acb}'>#</a>
**`Nemo.hyperu`** &mdash; *Method*.



```
hyperu(a::acb, b::acb, x::acb)
```

> Return the confluent hypergeometric function $U(a,b,x)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1269' class='documenter-source'>source</a><br>

<a id='Nemo.jtheta-Tuple{Nemo.acb,Nemo.acb}' href='#Nemo.jtheta-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Nemo.jtheta`** &mdash; *Method*.



```
jtheta(z::acb, tau::acb)
```

> Return a tuple of four elements containing the Jacobi theta function values $\theta_1, \theta_2, \theta_3, \theta_4$ evaluated at $z, \tau$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1280' class='documenter-source'>source</a><br>

<a id='Nemo.ellipwp-Tuple{Nemo.acb,Nemo.acb}' href='#Nemo.ellipwp-Tuple{Nemo.acb,Nemo.acb}'>#</a>
**`Nemo.ellipwp`** &mdash; *Method*.



```
ellipwp(z::acb, tau::acb)
```

> Return the Weierstrass elliptic function $\wp(z,\tau)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/bd26164b61cbd0978f038227c3d10d550157c1c2/src/arb/acb.jl#L1296' class='documenter-source'>source</a><br>


Here are some examples of complex valued mathematical functions.


```
CC = ComplexField(64)

s = CC(1, 2)
z = CC("1.23", "3.45")

a = sin(z)^2 + cos(z)^2
b = zeta(z)
c = besselj(s, z)
d = hyp1f1(s, s+1, z)
```

