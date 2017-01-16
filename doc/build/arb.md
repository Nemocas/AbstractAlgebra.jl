


<a id='Introduction-1'></a>

## Introduction


Arbitrary precision real ball arithmetic is supplied by Arb which provides a ball representation which tracks error bounds rigorously. Real numbers are  represented in mid-rad interval form $[m \pm r] = [m-r, m+r]$.


The Arb real field is constructed using the `ArbField` constructor. This constructs the parent object for the Arb real field.


However, we define


```
RealField = ArbField
```


so that one can construct the Arb real field parent object using `RealField` instead of `ArbField`.


The types of real balls in Nemo are given in the following table, along with the libraries that provide them and the associated types of the parent objects.


| Library |                Field | Element type | Parent type |
| -------:| --------------------:| ------------:| -----------:|
|     Arb | $\mathbb{R}$ (balls) |        `arb` |  `ArbField` |


All the real field types belong to the `Field` abstract type and the types of elements in this field, i.e. balls in this case, belong to the `FieldElem` abstract type.


<a id='Real-field-constructors-1'></a>

## Real field constructors


In order to construct real balls in Nemo, one must first construct the Arb real field itself. This is accomplished with the following constructor.


```
ArbField(prec::Int)
```


Return the Arb field with precision in bits `prec` used for operations on interval midpoints. The precision used for interval radii is a fixed implementation-defined constant (30 bits).


We define


```
RealField = ArbField
```


so that one can use `RealField` in place of `ArbField`.


Here is an example of creating an Arb real field and using the resulting parent object to coerce values into the resulting field.


```
RR = RealField(64)

a = RR("0.25")
b = RR("0.1")
c = RR(0.5)
d = RR(12)
```


Note that whilst one can coerce double precision floating point values into an Arb real field, unless those values can be represented exactly in double precision, the resulting ball can't be any more precise than the double precision supplied.


If instead, values can be represented precisely using decimal arithmetic then one can supply them to Arb using a string. In this case, Arb will store them to the precision specified when creating the Arb field.


If the values can be stored precisely as a binary floating point number, Arb will store the values exactly. See the function `isexact` below for more information.


<a id='Real-ball-constructors-1'></a>

## Real ball constructors


Once an Arb real field is constructed, there are various ways to construct balls in that field.


Apart from coercing elements into the Arb real field as above, we offer the following functions.


```
zero(R::ArbField)
```

> Return exact zero in the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L57' class='documenter-source'>source</a><br>


```
one(R::ArbField)
```

> Return exact one in the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L63' class='documenter-source'>source</a><br>


```
ball(mid::arb, rad::arb)
```

> Constructs an `arb` enclosing the range $[m-|r|, m+|r|]$, given the pair $(m, r)$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L445' class='documenter-source'>source</a><br>


Here are some examples of constructing balls.


```
RR = RealField(64)

a = zero(RR)
b = one(RR)
c = ball(RR(3), RR("0.0001"))
```


<a id='Basic-functionality-1'></a>

## Basic functionality


The following basic functionality is provided by the default Arb real field implementation in Nemo, to support construction of generic rings over real fields. Any custom real field implementation in Nemo should provide analogues of these functions along with the usual arithmetic operations.


```
parent_type(::Type{arb})
```


Gives the type of the parent object of an Arb real field element.


```
elem_type(R::ArbField)
```


Given the parent object for a Arb field, return the type of elements of the field.


```
mul!(c::arb, a::arb, b::arb)
```


Multiply $a$ by $b$ and set the existing Arb field element $c$ to the result. This function is provided for performance reasons as it saves allocating a new object for the result and eliminates associated garbage collection.


```
addeq!(c::arb, a::arb)
```


In-place addition. Adds $a$ to $c$ and sets $c$ to the result. This function is provided for performance reasons as it saves allocating a new object for the result and eliminates associated garbage collection.


```
deepcopy(a::arb)
```


Return a copy of the Arb field element $a$, recursively copying the internal data. Arb field elements are mutable in Nemo so a shallow copy is not sufficient.


Given the parent object `R` for an Arb real field, the following coercion functions are provided to coerce various elements into the Arb field. Developers provide these by overloading the `call` operator for the real field parent objects.


```
R()
```


Coerce zero into the Arb field.


```
R(n::Integer)
R(f::fmpz)
R(q::fmpq)
```


Coerce an integer or rational value into the Arb field.


```
R(f::Float64)
R(f::BigFloat)
```


Coerce the given floating point number into the Arb field.


```
R(f::AbstractString)
```


Coerce the decimal number, given as a string, into the Arb field.


```
R(f::arb)
```


Take an Arb field element that is already in the Arb field and simply return it. A copy of the original is not made.


Here are some examples of coercing elements into the Arb field.


```
RR = RealField(64)

a = RR(3)
b = RR(QQ(2,3))
c = ball(RR(3), RR("0.0001"))
d = RR("3 +/- 0.0001")
f = RR("-1.24e+12345")
g = RR("nan +/- inf")
```


In addition to the above, developers of custom real field types must ensure that they provide the equivalent of the function `base_ring(R::ArbField)` which should return `Union{}`. In addition to this they should ensure that each real field element contains a field `parent` specifying the parent object of the real field element, or at least supply the equivalent of the function `parent(a::arb)` to return the parent object of a real field element.


<a id='Conversions-1'></a>

## Conversions


```
convert(::Type{Float64}, x::arb)
```

> Return the midpoint of $x$ rounded down to a machine double.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L92' class='documenter-source'>source</a><br>


<a id='Basic-manipulation-1'></a>

## Basic manipulation


Numerous functions are provided to manipulate Arb field elements. Also see the section on basic functionality above.


```
base_ring(R::ArbField)
```

> Returns `Union{}` since an Arb field does not depend on any other ring.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L39' class='documenter-source'>source</a><br>


```
base_ring(x::arb)
```

> Returns `Union{}` since an Arb field does not depend on any other ring.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L45' class='documenter-source'>source</a><br>


```
parent(x::arb)
```

> Return the parent of the given Arb field element.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L51' class='documenter-source'>source</a><br>


```
iszero(x::arb)
```

> Return `true` if $x$ is certainly zero, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L353' class='documenter-source'>source</a><br>


```
isnonzero(x::arb)
```

> Return `true` if $x$ is certainly not equal to zero, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L361' class='documenter-source'>source</a><br>


```
isone(x::arb)
```

> Return `true` if $x$ is certainly not equal to oneo, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L370' class='documenter-source'>source</a><br>


```
isfinite(x::arb)
```

> Return `true` if $x$ is finite, i.e. having finite midpoint and radius, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L379' class='documenter-source'>source</a><br>


```
isexact(x::arb)
```

> Return `true` if $x$ is exact, i.e. has zero radius, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L388' class='documenter-source'>source</a><br>


```
isint(x::arb)
```

> Return `true` if $x$ is an exact integer, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L397' class='documenter-source'>source</a><br>


```
ispositive(x::arb)
```

> Return `true` if $x$ is certainly positive, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L405' class='documenter-source'>source</a><br>


```
isnonnegative(x::arb)
```

> Return `true` if $x$ is certainly nonnegative, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L413' class='documenter-source'>source</a><br>


```
isnegative(x::arb)
```

> Return `true` if $x$ is certainly negative, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L421' class='documenter-source'>source</a><br>


```
isnonpositive(x::arb)
```

> Return `true` if $x$ is certainly nonpositive, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L429' class='documenter-source'>source</a><br>


```
midpoint(x::arb)
```

> Return the midpoint of the ball $x$ as an Arb ball.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L466' class='documenter-source'>source</a><br>


```
radius(x::arb)
```

> Return the radius of the ball $x$ as an Arb ball.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L456' class='documenter-source'>source</a><br>


```
accuracy_bits(x::arb)
```

> Return the relative accuracy of $x$ measured in bits, capped between `typemax(Int)` and `-typemax(Int)`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L71' class='documenter-source'>source</a><br>


Here are some examples of basic manipulation of Arb balls.


```
RR = RealField(64)

a = RR("1.2 +/- 0.001")
b = RR(3)

iszero(a)
isone(b)
ispositive(a)
isfinite(b)
isint(b)
isnegative(a)
c = radius(a)
d = midpoint(b)
f = accuracy_bits(a)
S = parent(a)
T = base_ring(RR)
```


<a id='Arithmetic-operations-1'></a>

## Arithmetic operations


Nemo provides all the standard field operations for Arb field elements, as follows. Note that division is represented by `//` since a field is its own fraction field and since exact division is not generally possible in an inexact field.


|           Function |      Operation |
| ------------------:| --------------:|
|          -(a::arb) |    unary minus |
|  +(a::arb, b::arb) |       addition |
|  -(a::arb, b::arb) |    subtraction |
|  *(a::arb, b::arb) | multiplication |
| //(a::arb, b::arb) |       division |
|  ^(a::arb, b::arb) |       powering |


In addition, the following ad hoc field operations are defined.


|               Function |      Operation |
| ----------------------:| --------------:|
|  +(a::arb, b::Integer) |       addition |
|  +(a::Integer, b::arb) |       addition |
|     +(a::arb, b::fmpz) |       addition |
|     +(a::fmpz, b::arb) |       addition |
|  -(a::arb, b::Integer) |    subtraction |
|  -(a::Integer, b::arb) |    subtraction |
|     -(a::arb, b::fmpz) |    subtraction |
|     -(a::fmpz, b::arb) |    subtraction |
|  *(a::arb, b::Integer) | multiplication |
|  *(a::Integer, b::arb) | multiplication |
|     *(a::arb, b::fmpz) | multiplication |
|     *(a::fmpz, b::arb) | multiplication |
| //(a::arb, b::Integer) |       division |
|    //(a::arb, b::fmpz) |       division |
| //(a::Integer, b::arb) |       division |
|    //(a::fmpz, b::arb) |       division |
|  ^(a::arb, b::Integer) |       powering |
|     ^(a::arb, b::fmpz) |       powering |
|     ^(a::arb, b::fmpq) |       powering |


Here are some examples of arithmetic operations on Arb balls.


```
RR = RealField(64)

x = RR(3)
y = RR(QQ(2,3))

a = x + y
b = x*y
d = 1//y - x
d = 3x + ZZ(100)//y
f = (x^2 + y^2) ^ QQ(1,2)
```


<a id='Containment-1'></a>

## Containment


It is often necessary to determine whether a given exact value or ball is contained in a given real ball or whether two balls overlap. The following functions are provided for this purpose.


```
overlaps(x::arb, y::arb)
```

> Returns `true` if any part of the ball $x$ overlaps any part of the ball $y$, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L130' class='documenter-source'>source</a><br>


```
contains(x::arb, y::arb)
```

> Returns `true` if the ball $x$ contains the ball $y$, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L188' class='documenter-source'>source</a><br>


```
contains(x::arb, y::Integer)
```

> Returns `true` if the ball $x$ contains the given integer value, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L170' class='documenter-source'>source</a><br>


```
contains(x::arb, y::fmpz)
```

> Returns `true` if the ball $x$ contains the given integer value, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L155' class='documenter-source'>source</a><br>


```
contains(x::arb, y::fmpq)
```

> Returns `true` if the ball $x$ contains the given rational value, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L145' class='documenter-source'>source</a><br>


```
contains(x::arb, y::BigFloat)
```

> Returns `true` if the ball $x$ contains the given floating point value,  otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L177' class='documenter-source'>source</a><br>


The following functions are also provided for determining if a ball intersects a certain part of the real number line.


```
contains_zero(x::arb)
```

> Returns `true` if the ball $x$ contains zero, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L198' class='documenter-source'>source</a><br>


```
contains_negative(x::arb)
```

> Returns `true` if the ball $x$ contains any negative value, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L207' class='documenter-source'>source</a><br>


```
contains_positive(x::arb)
```

> Returns `true` if the ball $x$ contains any positive value, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L217' class='documenter-source'>source</a><br>


```
contains_nonnegative(x::arb)
```

> Returns `true` if the ball $x$ contains any nonnegative value, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L227' class='documenter-source'>source</a><br>


```
contains_nonpositive(x::arb)
```

> Returns `true` if the ball $x$ contains any nonpositive value, otherwise return `false`.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L237' class='documenter-source'>source</a><br>


Here are some examples of testing containment.


```
RR = RealField(64)
x = RR("1 +/- 0.001")
y = RR("3")

overlaps(x, y)
contains(x, y)
contains(y, 3)
contains(x, ZZ(1)//2)
contains_zero(x)
contains_positive(y)
```


<a id='Comparison-1'></a>

## Comparison


Nemo provides a full range of comparison operations for Arb balls. Note that a ball is considered less than another ball if every value in the first ball is less than every value in the second ball, etc.


Firstly, we introduce an exact equality which is distinct from arithmetic equality. This is distinct from arithmetic equality implemented by `==`, which merely compares up to the minimum of the precisions of its operands.


```
isequal(x::arb, y::arb)
```

> Return `true` if the balls $x$ and $y$ are precisely equal, i.e. have the same midpoints and radii.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L253' class='documenter-source'>source</a><br>


A full range of functions is available for comparing balls, i.e. `==`, `!=`, `<`, `<=`, `>=`, `>`. In fact, all these are implemented directly in C. In the table below we document these as though only `==` and `isless` had been provided to Julia.


<a id='Function-1'></a>

## Function


`isless(x::arb, y::arb)`  `==(x::arb, y::arb)`


As well as these, we provide a full range of ad hoc comparison operators. Again, these are implemented directly in Julia, but we document them as though `isless` and `==` were provided.


<a id='Function-2'></a>

## Function


`==(x::arb, y::Integer)` `==(x::Integer, y::arb)` `==(x::arb, y::fmpz)` `==(x::fmpz, y::arb)` `==(x::arb, y::Float64)` `==(x::Float64, y::arb)` `isless(x::arb, y::Integer)` `isless(x::Integer, y::arb)` `isless(x::arb, y::fmpz)` `isless(x::fmpz, y::arb)` `isless(x::arb, y::Float64)` `isless(x::Float64, y::arb)`


Here are some examples of comparison.


```
RR = RealField(64)
x = RR("1 +/- 0.001")
y = RR("3")
z = RR("4")

isequal(y, z)
x <= z
x == 3
ZZ(3) < z
x != 1.23
```


<a id='Absolute-value-1'></a>

## Absolute value


```
abs(x::arb)
```

> Return the absolute value of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L719' class='documenter-source'>source</a><br>


Here are some examples of taking the absolute value.


```
RR = RealField(64)
x = RR("-1 +/- 0.001")

a = abs(x)
```


<a id='Inverse-1'></a>

## Inverse


```
inv(x::arb)
```

> Return the multiplicative inverse of $x$, i.e. $1/x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L735' class='documenter-source'>source</a><br>


Here are some examples of taking the inverse.


```
RR = RealField(64)
x = RR("-3 +/- 0.001")

a = inv(x)
```


<a id='Shifting-1'></a>

## Shifting


```
ldexp(x::arb, y::Int)
```

> Return $2^yx$. Note that $y$ can be positive, zero or negative.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L752' class='documenter-source'>source</a><br>


```
ldexp(x::arb, y::fmpz)
```

> Return $2^yx$. Note that $y$ can be positive, zero or negative.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L763' class='documenter-source'>source</a><br>


Here are some examples of shifting.


```
RR = RealField(64)
x = RR("-3 +/- 0.001")

a = ldexp(x, 23)
b = ldexp(x, -ZZ(15))
```


<a id='Miscellaneous-operations-1'></a>

## Miscellaneous operations


```
trim(x::arb)
```

> Return an `arb` interval containing $x$ but which may be more economical, by rounding off insignificant bits from the midpoint.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L780' class='documenter-source'>source</a><br>


```
unique_integer(x::arb)
```

> Return a pair where the first value is a boolean and the second is an `fmpz` integer. The boolean indicates whether the interval $x$ contains a unique integer. If this is the case, the second return value is set to this unique integer.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L791' class='documenter-source'>source</a><br>


```
setunion(x::arb, y::arb)
```

> Return an `arb` containing the union of the intervals represented by $x$ and $y$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L805' class='documenter-source'>source</a><br>


Here are some examples of miscellaneous operations.


```
RR = RealField(64)
x = RR("-3 +/- 0.001")
y = RR("2 +/- 0.5")

a = trim(x)
b, c = unique_integer(x)
d = setunion(x, y)
```


<a id='Constants-1'></a>

## Constants


```
const_pi(r::ArbField)
```

> Return $\pi = 3.14159\ldots$ as an element of $r$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L823' class='documenter-source'>source</a><br>


```
const_e(r::ArbField)
```

> Return $e = 2.71828\ldots$ as an element of $r$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L833' class='documenter-source'>source</a><br>


```
const_log2(r::ArbField)
```

> Return $\log(2) = 0.69314\ldots$ as an element of $r$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L843' class='documenter-source'>source</a><br>


```
const_log10(r::ArbField)
```

> Return $\log(10) = 2.302585\ldots$ as an element of $r$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L853' class='documenter-source'>source</a><br>


```
const_euler(r::ArbField)
```

> Return Euler's constant $\gamma = 0.577215\ldots$ as an element of $r$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L863' class='documenter-source'>source</a><br>


```
const_catalan(r::ArbField)
```

> Return Catalan's constant $C = 0.915965\ldots$ as an element of $r$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L873' class='documenter-source'>source</a><br>


```
const_khinchin(r::ArbField)
```

> Return Khinchin's constant $K = 2.685452\ldots$ as an element of $r$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L883' class='documenter-source'>source</a><br>


```
const_glaisher(r::ArbField)
```

> Return Glaisher's constant $A = 1.282427\ldots$ as an element of $r$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L893' class='documenter-source'>source</a><br>


Here are some examples of computing real constants.


```
RR = RealField(200)

a = const_pi(RR)
b = const_e(RR)
c = const_euler(RR)
d = const_glaisher(RR)
```


<a id='Mathematical-functions-1'></a>

## Mathematical functions


```
floor(x::arb)
```

> Compute the floor of $x$, i.e. the greatest integer not exceeding $x$, as an Arb.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L911' class='documenter-source'>source</a><br>


```
ceil(x::arb)
```

> Return the ceiling of $x$, i.e. the least integer not less than $x$, as an Arb.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L922' class='documenter-source'>source</a><br>


```
sqrt(x::arb)
```

> Return the square root of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L933' class='documenter-source'>source</a><br>


```
rsqrt(x::arb)
```

> Return the reciprocal of the square root of $x$, i.e. $1/\sqrt{x}$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L943' class='documenter-source'>source</a><br>


```
sqrt1pm1(x::arb)
```

> Return $\sqrt{1+x}-1$, evaluated accurately for small $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L953' class='documenter-source'>source</a><br>


```
log(x)
```

Compute the natural logarithm of `x`. Throws `DomainError` for negative `Real` arguments. Use complex negative arguments to obtain complex results.

There is an experimental variant in the `Base.Math.JuliaLibm` module, which is typically faster and more accurate.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/3c9d75391c72d7c32eea75ff187ce77b2d5effc8/base/docs/helpdb/Base.jl#L2563-2571' class='documenter-source'>source</a><br>


```
log(x::arb)
```

> Return the principal branch of the logarithm of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L963' class='documenter-source'>source</a><br>


```
log1p(x::arb)
```

> Return $\log(1+x)$, evaluated accurately for small $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L973' class='documenter-source'>source</a><br>


```
exp(x::arb)
```

> Return the exponential of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L983' class='documenter-source'>source</a><br>


```
expm1(x::arb)
```

> Return $\exp(x)-1$, evaluated accurately for small $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L993' class='documenter-source'>source</a><br>


```
sin(x::arb)
```

> Return the sine of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1003' class='documenter-source'>source</a><br>


```
cos(x::arb)
```

> Return the cosine of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1013' class='documenter-source'>source</a><br>


```
sinpi(x::arb)
```

> Return the sine of $\pi x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1023' class='documenter-source'>source</a><br>


```
cospi(x::arb)
```

> Return the cosine of $\pi x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1033' class='documenter-source'>source</a><br>


```
tan(x::arb)
```

> Return the tangent of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1043' class='documenter-source'>source</a><br>


```
cot(x::arb)
```

> Return the cotangent of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1053' class='documenter-source'>source</a><br>


```
tanpi(x::arb)
```

> Return the tangent of $\pi x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1063' class='documenter-source'>source</a><br>


```
cotpi(x::arb)
```

> Return the cotangent of $\pi x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1073' class='documenter-source'>source</a><br>


```
sinh(x::arb)
```

> Return the hyperbolic sine of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1083' class='documenter-source'>source</a><br>


```
cosh(x::arb)
```

> Return the hyperbolic cosine of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1093' class='documenter-source'>source</a><br>


```
tanh(x::arb)
```

> Return the hyperbolic tangent of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1103' class='documenter-source'>source</a><br>


```
coth(x::arb)
```

> Return the hyperbolic cotangent of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1113' class='documenter-source'>source</a><br>


```
atan(x::arb)
```

> Return the arctangent of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1123' class='documenter-source'>source</a><br>


```
asin(x::arb)
```

> Return the arcsine of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1133' class='documenter-source'>source</a><br>


```
acos(x::arb)
```

> Return the arccosine of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1143' class='documenter-source'>source</a><br>


```
atanh(x::arb)
```

> Return the hyperbolic arctangent of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1153' class='documenter-source'>source</a><br>


```
asinh(x::arb)
```

> Return the hyperbolic arcsine of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1163' class='documenter-source'>source</a><br>


```
acosh(x::arb)
```

> Return the hyperbolic arccosine of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1173' class='documenter-source'>source</a><br>


```
gamma(x::arb)
```

> Return the Gamma function evaluated at $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1183' class='documenter-source'>source</a><br>


```
lgamma(x::arb)
```

> Return the logarithm of the Gamma function evaluated at $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1193' class='documenter-source'>source</a><br>


```
rgamma(x::arb)
```

> Return the reciprocal of the Gamma function evaluated at $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1203' class='documenter-source'>source</a><br>


```
digamma(x::arb)
```

> Return the  logarithmic derivative of the gamma function evaluated at $x$, i.e. $\psi(x)$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1213' class='documenter-source'>source</a><br>


```
zeta(s)
```

Riemann zeta function $\zeta(s)$.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/3c9d75391c72d7c32eea75ff187ce77b2d5effc8/base/docs/helpdb/Base.jl#L6858-6862' class='documenter-source'>source</a><br>


```
zeta(x::arb)
```

> Return the Riemann zeta function evaluated at $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1224' class='documenter-source'>source</a><br>


```
sincos(x::arb)
```

> Return a tuple $s, c$ consisting of the sine $s$ and cosine $c$ of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1234' class='documenter-source'>source</a><br>


```
sincospi(x::arb)
```

> Return a tuple $s, c$ consisting of the sine $s$ and cosine $c$ of $\pi x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1246' class='documenter-source'>source</a><br>


```
sinpi(x::fmpq, r::ArbField)
```

> Return the sine of $\pi x$ in the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1258' class='documenter-source'>source</a><br>


```
cospi(x::fmpq, r::ArbField)
```

> Return the cosine of $\pi x$ in the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1269' class='documenter-source'>source</a><br>


```
sincospi(x::fmpq, r::ArbField)
```

> Return a tuple $s, c$ consisting of the sine and cosine of $\pi x$ in the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1280' class='documenter-source'>source</a><br>


```
sinhcosh(x::arb)
```

> Return a tuple $s, c$ consisting of the hyperbolic sine and cosine of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1293' class='documenter-source'>source</a><br>


```
atan2(x::arb, y::arb)
```

> Return atan2$(b,a) = \arg(a+bi)$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1305' class='documenter-source'>source</a><br>


```
agm(x::arb, y::arb)
```

> Return the arithmetic-geometric mean of $x$ and $y$



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1316' class='documenter-source'>source</a><br>


```
zeta(s, z)
```

Generalized zeta function $\zeta(s, z)$, defined by the sum $\sum_{k=0}^\infty ((k+z)^2)^{-s/2}$, where any term with $k+z=0$ is excluded.  For $\Re z > 0$, this definition is equivalent to the Hurwitz zeta function $\sum_{k=0}^\infty (k+z)^{-s}$.   For $z=1$, it yields the Riemann zeta function $\zeta(s)$.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/3c9d75391c72d7c32eea75ff187ce77b2d5effc8/base/special/gamma.jl#L279-288' class='documenter-source'>source</a><br>


```
zeta(s::arb, a::arb)
```

> Return the Hurwitz zeta function $\zeta(s,a)$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1327' class='documenter-source'>source</a><br>


```
hypot(x::arb, y::arb)
```

> Return $\sqrt{x^2 + y^2}$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1338' class='documenter-source'>source</a><br>


```
root(x::arb, n::Int)
```

> Return the $n$-th root of $x$. We require $x \geq 0$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1356' class='documenter-source'>source</a><br>


```
fac(x::arb)
```

> Return the factorial of $x$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1362' class='documenter-source'>source</a><br>


```
fac(n::Int, r::ArbField)
```

> Return the factorial of $n$ in the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1374' class='documenter-source'>source</a><br>


```
binom(x::arb, n::UInt)
```

> Return the binomial coefficient ${x \choose n}$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1380' class='documenter-source'>source</a><br>


```
binom(n::UInt, k::UInt, r::ArbField)
```

> Return the binomial coefficient ${n \choose k}$ in the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1391' class='documenter-source'>source</a><br>


```
fib(n::fmpz, r::ArbField)
```

> Return the $n$-th Fibonacci number in the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1402' class='documenter-source'>source</a><br>


```
fib(n::Int, r::ArbField)
```

> Return the $n$-th Fibonacci number in the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1420' class='documenter-source'>source</a><br>


```
gamma(x::fmpz, r::ArbField)
```

> Return the Gamma function evaluated at $x$ in the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1426' class='documenter-source'>source</a><br>


```
gamma(x::fmpq, r::ArbField)
```

> Return the Gamma function evaluated at $x$ in the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1437' class='documenter-source'>source</a><br>


```
zeta(s, z)
```

Generalized zeta function $\zeta(s, z)$, defined by the sum $\sum_{k=0}^\infty ((k+z)^2)^{-s/2}$, where any term with $k+z=0$ is excluded.  For $\Re z > 0$, this definition is equivalent to the Hurwitz zeta function $\sum_{k=0}^\infty (k+z)^{-s}$.   For $z=1$, it yields the Riemann zeta function $\zeta(s)$.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/3c9d75391c72d7c32eea75ff187ce77b2d5effc8/base/special/gamma.jl#L279-288' class='documenter-source'>source</a><br>


```
zeta(n::Int, r::ArbField)
```

> Return the Riemann zeta function $\zeta(n)$ as an element of the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1456' class='documenter-source'>source</a><br>


```
bernoulli(n::Int, r::ArbField)
```

> Return the $n$-th Bernoulli number as an element of the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1470' class='documenter-source'>source</a><br>


```
risingfac(x::arb, n::Int)
```

> Return the rising factorial $x(x + 1)\ldots (x + n - 1)$ as an Arb.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1483' class='documenter-source'>source</a><br>


```
risingfac(x::fmpq, n::Int, r::ArbField)
```

> Return the rising factorial $x(x + 1)\ldots (x + n - 1)$ as an element of the given Arb field.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1496' class='documenter-source'>source</a><br>


```
risingfac2(x::arb, n::Int)
```

> Return a tuple containing the rising factorial $x(x + 1)\ldots (x + n - 1)$ and its derivative.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1511' class='documenter-source'>source</a><br>


```
polylog(s::arb, a::arb)
```

> Return the polylogarithm Li$_s(a)$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1518' class='documenter-source'>source</a><br>


```
polylog(s::Int, a::arb)
```

> Return the polylogarithm Li$_s(a)$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1529' class='documenter-source'>source</a><br>


```
chebyshev_t(n::Int, x::arb)
```

> Return the value of the Chebyshev polynomial $T_n(x)$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1570' class='documenter-source'>source</a><br>


```
chebyshev_u(n::Int, x::arb)
```

> Return the value of the Chebyshev polynomial $U_n(x)$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1576' class='documenter-source'>source</a><br>


```
chebyshev_t2(n::Int, x::arb)
```

> Return the tuple $(T_{n}(x), T_{n-1}(x))$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1582' class='documenter-source'>source</a><br>


```
chebyshev_u2(n::Int, x::arb)
```

> Return the tuple $(U_{n}(x), U_{n-1}(x))$



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1588' class='documenter-source'>source</a><br>


```
bell(n::fmpz, r::ArbField)
```

> Return the Bell number $B_n$ as an element of $r$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1594' class='documenter-source'>source</a><br>


```
bell(n::Int, r::ArbField)
```

> Return the Bell number $B_n$ as an element of $r$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1605' class='documenter-source'>source</a><br>


```
numpart(n::fmpz, r::ArbField)
```

> Return the number of partitions $p(n)$ as an element of $r$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1611' class='documenter-source'>source</a><br>


```
numpart(n::fmpz, r::ArbField)
```

> Return the number of partitions $p(n)$ as an element of $r$.



<a target='_blank' href='https://github.com/Nemocas/Nemo.jl/tree/4ead26d440e2c61322542a1a216fa6e66c8617df/src/arb/arb.jl#L1622' class='documenter-source'>source</a><br>


Here are some examples of real valued mathematical functions.


```
RR = RealField(64)

a = floor(exp(RR(1)))
b = sinpi(QQ(5,6), RR)
c = gamma(QQ(1,3), RealField(256))
d = bernoulli(1000, RealField(53))
f = polylog(3, RR(-10))
```

