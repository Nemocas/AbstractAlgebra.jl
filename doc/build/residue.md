


<a id='Introduction-1'></a>

## Introduction


Nemo allows the creation of residue rings of the form $R/(a)$ for an element $a$ of a ring $R$.


We don't require $(a)$ to be a prime or maximal ideal. Instead, we allow the creation of the residue ring $R/(a)$ for any nonzero $a$ and simply raise an exception if an impossible inverse is encountered during computations  involving elements of $R/(a)$.


There is a generic implementation of residue rings of this form in Nemo, which accepts any ring $R$ as base ring. The associated types of parent object and elements in such residue rings are given in the following table.


|        Base ring | Library | Element type |    Parent type |
| ----------------:| -------:| ------------:| --------------:|
| Generic ring $R$ |    Nemo |  `GenRes{T}` | `GeResRing{T}` |


The modulus $a$ of a residue ring is stored in its parent object.


All residue element types belong to the abstract type `ResElem` and all the residue ring parent object types belong to the abstract type `ResRing`. This enables one to write generic functions that accept any Nemo residue type.


<a id='Residue-ring-constructors-1'></a>

## Residue ring constructors


In order to construct residues in Nemo, one must first construct the residue ring itself. This is accomplished with the following constructor.

<a id='Nemo.ResidueRing-Tuple{Nemo.Ring,T<:Nemo.RingElem}' href='#Nemo.ResidueRing-Tuple{Nemo.Ring,T<:Nemo.RingElem}'>#</a>
**`Nemo.ResidueRing`** &mdash; *Method*.



```
ResidueRing{T <: RingElem}(R::Ring, a::T; cached=true)
```

> Create the residue ring $R/(a)$ where $a$ is an element of the ring $R$. We require $a \neq 0$. If `cached == true` (the default) then the resulting residue ring parent object is cached and returned for any subsequent calls to the constructor with the same base ring $R$ and element $a$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L549' class='documenter-source'>source</a><br>


Here are some examples of creating residue rings and making use of the resulting parent objects to coerce various elements into those rings.


```
R = ResidueRing(ZZ, 17)
S, x = PolynomialRing(QQ, "x")
T = ResidueRing(S, x^2 + 1)

a = R(11)
b = T(x)
c = T(ZZ(3))
```


<a id='Residue-element-constructors-1'></a>

## Residue element constructors


Once a residue ring is constructed, there are various ways to construct residues in that ring.


Apart from coercing elements into the residue ring as above, we offer the following functions.

<a id='Base.zero-Tuple{Nemo.ResRing}' href='#Base.zero-Tuple{Nemo.ResRing}'>#</a>
**`Base.zero`** &mdash; *Method*.



```
zero(R::ResRing)
```

> Return the zero element of the given residue ring, i.e. $0 \pmod{a}$ where $a$ is the modulus of the residue ring.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L79' class='documenter-source'>source</a><br>

<a id='Base.one-Tuple{Nemo.ResRing}' href='#Base.one-Tuple{Nemo.ResRing}'>#</a>
**`Base.one`** &mdash; *Method*.



```
zero(R::ResRing)
```

> Return $1 \pmod{a}$ where $a$ is the modulus of the residue ring.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L86' class='documenter-source'>source</a><br>


Here are some examples of constructing residues.


```
R = ResidueRing(ZZ, 16453889)
S, x = PolynomialRing(R, "x")
T = ResidueRing(S, x^3 + 3x + 1)

g = zero(R)
h = one(T)
```


<a id='Basic-functionality-1'></a>

## Basic functionality


All residue ring modules in Nemo must provide the functionality listed in this section. (Note that only some of these functions are useful to a user.)


Developers who are writing their own residue ring module, whether as an interface to a C library, or as some kind of generic module, must provide all of these functions for custom residue ring types in Nemo. 


We write `U` for the type of residue elements in the residue ring and `T` for the type of elements of the base ring.


All of these functions are provided for all existing residue types in Nemo.


```
parent_type{U <: ResElem}(::Type{U})
```


Given the type of residue elements, should return the type of the corresponding parent object.


```
elem_type(R::ResRing)
```


Given a parent object for the residue ring, return the type of elements of the residue ring.


```
Base.hash(a::ResElem, h::UInt)
```


Return a `UInt` hexadecimal hash of the residue $a$. This should be xor'd with a fixed random hexadecimal specific to the residue type. The hash of the data associated with a residue should be xor'd with the supplied parameter `h` as part of computing the hash.


```
data(a::ResElem)
```


Return the data of the given residue element, i.e. for $c \pmod{a}$ in the ring $R/(a)$ return the value $c$ as an element of the ring $R$.


```
deepcopy(a::ResElem)
```


Construct a copy of the given residue and return it. This function must recursively construct copies of all of the internal data in the given residue. Nemo residues are mutable and so returning shallow copies is not sufficient.


```
mul!(c::ResElem, a::ResElem, b::ResElem)
```


Multiply $a$ by $b$ and set the existing residue $c$ to the result. This function is provided for performance reasons as it saves allocating a new object for the result and eliminates associated garbage collection.


```
addeq!(c::ResElem, a::ResElem)
```


In-place addition. Adds $a$ to $c$ and sets $c$ to the result. This function is provided for performance reasons as it saves allocating a new object for the result and eliminates associated garbage collection.


Given a parent object `S` for a residue ring, the following coercion functions are provided to coerce various elements into the residue ring. Developers provide these by overloading the `call` operator for the residue ring parent objects.


```
S()
```


Coerce zero into the ring $S$.


```
S(n::Integer)
S(n::fmpz)
```


Coerce an integer value or Flint integer into the residue ring $S$.


```
S(n::T)
```


Coerces an element of the base ring, of type `T` into $S$.


```
S(f::ResElem)
```


Take a residue that is already in the ring $S$ and simply return it. A copy of the original is not made.


```
S(c::RingElem)
```


Try to coerce the given ring element into the residue ring. This only succeeds if $c$ can be coerced into the base ring.


In addition to the above, developers of custom residues must ensure the parent object of a residue type constains a field `base_ring` specifying the base ring, a field `modulus` containing the modulus $a$ of the residue ring $R/(a)$. They must also ensure that each residue element contains a field `parent` specifying the parent object of the residue.


Typically a developer will also overload the `ResidueRing` generic function to create residues of the custom type they are implementing.


<a id='Basic-manipulation-1'></a>

## Basic manipulation


Numerous functions are provided to manipulate residues. Also see the section on basic functionality above.

<a id='Nemo.base_ring-Tuple{Nemo.ResRing}' href='#Nemo.base_ring-Tuple{Nemo.ResRing}'>#</a>
**`Nemo.base_ring`** &mdash; *Method*.


<a id='Nemo.base_ring-Tuple{Nemo.ResElem}' href='#Nemo.base_ring-Tuple{Nemo.ResElem}'>#</a>
**`Nemo.base_ring`** &mdash; *Method*.



```
base_ring(r::ResElem)
```

> Return the base ring $R$ of the residue ring $R/(a)$ that the supplied element $r$ belongs to.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L25' class='documenter-source'>source</a><br>

<a id='Base.parent-Tuple{Nemo.ResElem}' href='#Base.parent-Tuple{Nemo.ResElem}'>#</a>
**`Base.parent`** &mdash; *Method*.



```
parent(a::ResElem)
```

> Return the parent object of the given residue element.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L32' class='documenter-source'>source</a><br>

<a id='Nemo.modulus-Tuple{Nemo.ResRing}' href='#Nemo.modulus-Tuple{Nemo.ResRing}'>#</a>
**`Nemo.modulus`** &mdash; *Method*.



```
modulus(R::ResRing)
```

> Return the modulus $a$ of the given residue ring $S = R/(a)$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L60' class='documenter-source'>source</a><br>

<a id='Nemo.modulus-Tuple{Nemo.ResElem}' href='#Nemo.modulus-Tuple{Nemo.ResElem}'>#</a>
**`Nemo.modulus`** &mdash; *Method*.



```
modulus(R::ResRing)
```

> Return the modulus $a$ of the residue ring $S = R/(a)$ that the supplied residue $r$ belongs to.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L68' class='documenter-source'>source</a><br>

<a id='Nemo.iszero-Tuple{Nemo.ResElem}' href='#Nemo.iszero-Tuple{Nemo.ResElem}'>#</a>
**`Nemo.iszero`** &mdash; *Method*.



```
iszero(a::ResElem)
```

> Return `true` if the supplied element $a$ is zero in the residue ring it belongs to, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L92' class='documenter-source'>source</a><br>

<a id='Nemo.isone-Tuple{Nemo.ResElem}' href='#Nemo.isone-Tuple{Nemo.ResElem}'>#</a>
**`Nemo.isone`** &mdash; *Method*.



```
isone(a::ResElem)
```

> Return `true` if the supplied element $a$ is one in the residue ring it belongs to, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L99' class='documenter-source'>source</a><br>

<a id='Nemo.isunit-Tuple{Nemo.ResElem}' href='#Nemo.isunit-Tuple{Nemo.ResElem}'>#</a>
**`Nemo.isunit`** &mdash; *Method*.



```
iszero(a::ResElem)
```

> Return `true` if the supplied element $a$ is invertible in the residue ring it belongs to, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L106' class='documenter-source'>source</a><br>


Here are some examples of basic manipulation of residues.


```
R = ResidueRing(ZZ, 16453889)
S, x = PolynomialRing(R, "x")
T = ResidueRing(S, x^3 + 3x + 1)

g = zero(R)
h = one(T)

a = modulus(R)
b = modulus(T(x + 1))
c = isunit(h)
d = isone(h)
f = iszero(g)
U = parent(g)
V = base_ring(R)
W = base_ring(T(x + 1))
```


<a id='Arithmetic-operators-1'></a>

## Arithmetic operators


All the usual arithmetic operators are overloaded for Nemo residues. Note that Julia uses the single slash for floating point division. Therefore to perform exact division in a ring we use `divexact`. To construct an element of a fraction field one can use the double slash operator `//`.

<a id='Base.--Tuple{Nemo.ResElem}' href='#Base.--Tuple{Nemo.ResElem}'>#</a>
**`Base.-`** &mdash; *Method*.



```
-(x)
```

Unary minus operator.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L636-640' class='documenter-source'>source</a><br>

<a id='Base.+-Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}' href='#Base.+-Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}'>#</a>
**`Base.+`** &mdash; *Method*.



```
+{T <: RingElem}(a::ResElem{T}, b::ResElem{T})
```

> Return $a + b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L166' class='documenter-source'>source</a><br>


```
+{T <: RingElem}(a::ResElem{T}, b::T)
```

> Return $a + b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L247' class='documenter-source'>source</a><br>


```
+{T <: RingElem}(a::T, b::ResElem{T})
```

> Return $a + b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L265' class='documenter-source'>source</a><br>

<a id='Base.--Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}' href='#Base.--Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}'>#</a>
**`Base.-`** &mdash; *Method*.



```
-(x, y)
```

Subtraction operator.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L643-647' class='documenter-source'>source</a><br>


```
-{T <: RingElem}(a::ResElem{T}, b::ResElem{T})
```

> Return $a - b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L175' class='documenter-source'>source</a><br>


```
-{T <: RingElem}(a::ResElem{T}, b::T)
```

> Return $a - b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L283' class='documenter-source'>source</a><br>


```
-{T <: RingElem}(a::T, b::ResElem{T})
```

> Return $a - b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L301' class='documenter-source'>source</a><br>

<a id='Base.*-Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}' href='#Base.*-Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}'>#</a>
**`Base.*`** &mdash; *Method*.



```
*(x, y...)
```

Multiplication operator. `x*y*z*...` calls this function with all arguments, i.e. `*(x, y, z, ...)`.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L7701-7707' class='documenter-source'>source</a><br>


```
*{T <: RingElem}(a::ResElem{T}, b::ResElem{T})
```

> Return $a\times b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L184' class='documenter-source'>source</a><br>


```
*{T <: RingElem}(a::ResElem{T}, b::T)
```

> Return $a\times b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L211' class='documenter-source'>source</a><br>


```
*{T <: RingElem}(a::T, b::ResElem{T})
```

> Return $a\times b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L229' class='documenter-source'>source</a><br>

<a id='Nemo.divexact-Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}' href='#Nemo.divexact-Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}'>#</a>
**`Nemo.divexact`** &mdash; *Method*.



```
divexact{T <: RingElem}(a::ResElem{T}, b::ResElem{T})
```

> Return $a/b$ where the quotient is expected to be exact.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L441' class='documenter-source'>source</a><br>


The following ad hoc operators are also provided.

<a id='Base.+-Tuple{Integer,Nemo.ResElem}' href='#Base.+-Tuple{Integer,Nemo.ResElem}'>#</a>
**`Base.+`** &mdash; *Method*.



```
+(a::Integer, b::ResElem)
```

> Return $a + b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L253' class='documenter-source'>source</a><br>

<a id='Base.+-Tuple{Nemo.ResElem,Integer}' href='#Base.+-Tuple{Nemo.ResElem,Integer}'>#</a>
**`Base.+`** &mdash; *Method*.



```
+(a::ResElem, b::Integer)
```

> Return $a + b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L235' class='documenter-source'>source</a><br>

<a id='Base.+-Tuple{Nemo.fmpz,Nemo.ResElem}' href='#Base.+-Tuple{Nemo.fmpz,Nemo.ResElem}'>#</a>
**`Base.+`** &mdash; *Method*.



```
+(a::fmpz, b::ResElem)
```

> Return $a + b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L259' class='documenter-source'>source</a><br>

<a id='Base.+-Tuple{Nemo.ResElem,Nemo.fmpz}' href='#Base.+-Tuple{Nemo.ResElem,Nemo.fmpz}'>#</a>
**`Base.+`** &mdash; *Method*.



```
+(a::ResElem, b::Integer)
```

> Return $a + b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L241' class='documenter-source'>source</a><br>

<a id='Base.+-Tuple{T<:Nemo.RingElem,Nemo.ResElem{T<:Nemo.RingElem}}' href='#Base.+-Tuple{T<:Nemo.RingElem,Nemo.ResElem{T<:Nemo.RingElem}}'>#</a>
**`Base.+`** &mdash; *Method*.



```
+{T <: RingElem}(a::T, b::ResElem{T})
```

> Return $a + b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L265' class='documenter-source'>source</a><br>

<a id='Base.+-Tuple{Nemo.ResElem{T<:Nemo.RingElem},T<:Nemo.RingElem}' href='#Base.+-Tuple{Nemo.ResElem{T<:Nemo.RingElem},T<:Nemo.RingElem}'>#</a>
**`Base.+`** &mdash; *Method*.



```
+{T <: RingElem}(a::ResElem{T}, b::T)
```

> Return $a + b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L247' class='documenter-source'>source</a><br>

<a id='Base.--Tuple{Integer,Nemo.ResElem}' href='#Base.--Tuple{Integer,Nemo.ResElem}'>#</a>
**`Base.-`** &mdash; *Method*.



```
-(x, y)
```

Subtraction operator.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L643-647' class='documenter-source'>source</a><br>


```
-(a::Integer, b::ResElem)
```

> Return $a - b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L289' class='documenter-source'>source</a><br>

<a id='Base.--Tuple{Nemo.ResElem,Integer}' href='#Base.--Tuple{Nemo.ResElem,Integer}'>#</a>
**`Base.-`** &mdash; *Method*.



```
-(x, y)
```

Subtraction operator.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L643-647' class='documenter-source'>source</a><br>


```
-(a::ResElem, b::Integer)
```

> Return $a - b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L271' class='documenter-source'>source</a><br>

<a id='Base.--Tuple{Nemo.fmpz,Nemo.ResElem}' href='#Base.--Tuple{Nemo.fmpz,Nemo.ResElem}'>#</a>
**`Base.-`** &mdash; *Method*.



```
-(x, y)
```

Subtraction operator.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L643-647' class='documenter-source'>source</a><br>


```
-(a::fmpz, b::ResElem)
```

> Return $a - b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L295' class='documenter-source'>source</a><br>

<a id='Base.--Tuple{Nemo.ResElem,Nemo.fmpz}' href='#Base.--Tuple{Nemo.ResElem,Nemo.fmpz}'>#</a>
**`Base.-`** &mdash; *Method*.



```
-(x, y)
```

Subtraction operator.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L643-647' class='documenter-source'>source</a><br>


```
-(a::ResElem, b::fmpz)
```

> Return $a - b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L277' class='documenter-source'>source</a><br>

<a id='Base.--Tuple{T<:Nemo.RingElem,Nemo.ResElem{T<:Nemo.RingElem}}' href='#Base.--Tuple{T<:Nemo.RingElem,Nemo.ResElem{T<:Nemo.RingElem}}'>#</a>
**`Base.-`** &mdash; *Method*.



```
-(x, y)
```

Subtraction operator.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L643-647' class='documenter-source'>source</a><br>


```
-{T <: RingElem}(a::T, b::ResElem{T})
```

> Return $a - b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L301' class='documenter-source'>source</a><br>

<a id='Base.--Tuple{Nemo.ResElem{T<:Nemo.RingElem},T<:Nemo.RingElem}' href='#Base.--Tuple{Nemo.ResElem{T<:Nemo.RingElem},T<:Nemo.RingElem}'>#</a>
**`Base.-`** &mdash; *Method*.



```
-(x, y)
```

Subtraction operator.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L643-647' class='documenter-source'>source</a><br>


```
-{T <: RingElem}(a::ResElem{T}, b::T)
```

> Return $a - b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L283' class='documenter-source'>source</a><br>

<a id='Base.*-Tuple{Integer,Nemo.ResElem}' href='#Base.*-Tuple{Integer,Nemo.ResElem}'>#</a>
**`Base.*`** &mdash; *Method*.



```
*(x, y...)
```

Multiplication operator. `x*y*z*...` calls this function with all arguments, i.e. `*(x, y, z, ...)`.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L7701-7707' class='documenter-source'>source</a><br>


```
*(a::Integer, b::ResElem)
```

> Return $a\times b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L217' class='documenter-source'>source</a><br>

<a id='Base.*-Tuple{Nemo.ResElem,Integer}' href='#Base.*-Tuple{Nemo.ResElem,Integer}'>#</a>
**`Base.*`** &mdash; *Method*.



```
*(x, y...)
```

Multiplication operator. `x*y*z*...` calls this function with all arguments, i.e. `*(x, y, z, ...)`.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L7701-7707' class='documenter-source'>source</a><br>


```
*(a::ResElem, b::Integer)
```

> Return $a\times b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L199' class='documenter-source'>source</a><br>

<a id='Base.*-Tuple{Nemo.fmpz,Nemo.ResElem}' href='#Base.*-Tuple{Nemo.fmpz,Nemo.ResElem}'>#</a>
**`Base.*`** &mdash; *Method*.



```
*(x, y...)
```

Multiplication operator. `x*y*z*...` calls this function with all arguments, i.e. `*(x, y, z, ...)`.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L7701-7707' class='documenter-source'>source</a><br>


```
*(a::fmpz, b::ResElem)
```

> Return $a\times b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L223' class='documenter-source'>source</a><br>

<a id='Base.*-Tuple{Nemo.ResElem,Nemo.fmpz}' href='#Base.*-Tuple{Nemo.ResElem,Nemo.fmpz}'>#</a>
**`Base.*`** &mdash; *Method*.



```
*(x, y...)
```

Multiplication operator. `x*y*z*...` calls this function with all arguments, i.e. `*(x, y, z, ...)`.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L7701-7707' class='documenter-source'>source</a><br>


```
*(a::ResElem, b::fmpz)
```

> Return $a\times b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L205' class='documenter-source'>source</a><br>

<a id='Base.*-Tuple{T<:Nemo.RingElem,Nemo.ResElem{T<:Nemo.RingElem}}' href='#Base.*-Tuple{T<:Nemo.RingElem,Nemo.ResElem{T<:Nemo.RingElem}}'>#</a>
**`Base.*`** &mdash; *Method*.



```
*(x, y...)
```

Multiplication operator. `x*y*z*...` calls this function with all arguments, i.e. `*(x, y, z, ...)`.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L7701-7707' class='documenter-source'>source</a><br>


```
*{T <: RingElem}(a::T, b::ResElem{T})
```

> Return $a\times b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L229' class='documenter-source'>source</a><br>

<a id='Base.*-Tuple{Nemo.ResElem{T<:Nemo.RingElem},T<:Nemo.RingElem}' href='#Base.*-Tuple{Nemo.ResElem{T<:Nemo.RingElem},T<:Nemo.RingElem}'>#</a>
**`Base.*`** &mdash; *Method*.



```
*(x, y...)
```

Multiplication operator. `x*y*z*...` calls this function with all arguments, i.e. `*(x, y, z, ...)`.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L7701-7707' class='documenter-source'>source</a><br>


```
*{T <: RingElem}(a::ResElem{T}, b::T)
```

> Return $a\times b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L211' class='documenter-source'>source</a><br>

<a id='Base.^-Tuple{Nemo.ResElem,Int64}' href='#Base.^-Tuple{Nemo.ResElem,Int64}'>#</a>
**`Base.^`** &mdash; *Method*.



```
^(x, y)
```

Exponentiation operator.


<a target='_blank' href='https://github.com/JuliaLang/julia/tree/55e3a39579696345027d0d8ae489825c9d9201ab/base/docs/helpdb/Base.jl#L4495-4499' class='documenter-source'>source</a><br>


```
^(a::ResElem, b::Int)
```

> Return $a^b$.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L313' class='documenter-source'>source</a><br>


If the appropriate `promote_rule` and coercion exists, these operators can also be used with elements of other rings. Nemo will try to coerce the operands to the dominating type and then apply the operator.


Here are some examples of arithmetic operations on residues.


```
R = ResidueRing(ZZ, 16453889)

S, x = PolynomialRing(R, "x")
T = ResidueRing(S, x^3 + 3x + 1)

f = R(4)
g = R(6)
n = T(x^5 + 1)
p = T(x^2 + 2x + 1)

a = -R(12345)
b = -T(x^5 + 1)
c = f + g
d = f - g
h = f*g
k = n + p
l = n - p
m = n*p
q = divexact(a*b, a)
r = divexact(n*p, p)
s = f + 3
t = 3 - g
u = 5g
v = n + 4
w = 4 - p
y = n*5
z = p^12
```


<a id='Comparison-operators-1'></a>

## Comparison operators


The following comparison operators are implemented for residues in Nemo.

<a id='Base.==-Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}' href='#Base.==-Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}'>#</a>
**`Base.==`** &mdash; *Method*.



```
=={T <: RingElem}(x::ResElem{T}, y::ResElem{T})
```

> Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall that power series to different precisions may still be arithmetically equal to the minimum of the two precisions.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L327' class='documenter-source'>source</a><br>


```
=={T <: RingElem}(x::ResElem{T}, y::T)
```

> Return `true` if $x == y$ arithmetically, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L392' class='documenter-source'>source</a><br>


```
=={T <: RingElem}(x::T, y::ResElem{T})
```

> Return `true` if $x == y$ arithmetically, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L401' class='documenter-source'>source</a><br>

<a id='Base.isequal-Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}' href='#Base.isequal-Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}'>#</a>
**`Base.isequal`** &mdash; *Method*.



```
isequal{T <: RingElem}(x::ResElem{T}, y::ResElem{T})
```

> Return `true` if $x == y$ exactly, otherwise return `false`. This function is useful in cases where the data of the residues are inexact, e.g. power series Only if the power series are precisely the same, to the same precision, are they declared equal by this function.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L338' class='documenter-source'>source</a><br>


In addition we have the following ad hoc comparison operators.

<a id='Base.==-Tuple{Nemo.ResElem{T<:Nemo.RingElem},T<:Nemo.RingElem}' href='#Base.==-Tuple{Nemo.ResElem{T<:Nemo.RingElem},T<:Nemo.RingElem}'>#</a>
**`Base.==`** &mdash; *Method*.



```
=={T <: RingElem}(x::ResElem{T}, y::T)
```

> Return `true` if $x == y$ arithmetically, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L392' class='documenter-source'>source</a><br>

<a id='Base.==-Tuple{T<:Nemo.RingElem,Nemo.ResElem{T<:Nemo.RingElem}}' href='#Base.==-Tuple{T<:Nemo.RingElem,Nemo.ResElem{T<:Nemo.RingElem}}'>#</a>
**`Base.==`** &mdash; *Method*.



```
=={T <: RingElem}(x::T, y::ResElem{T})
```

> Return `true` if $x == y$ arithmetically, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L401' class='documenter-source'>source</a><br>

<a id='Base.==-Tuple{Nemo.ResElem,Integer}' href='#Base.==-Tuple{Nemo.ResElem,Integer}'>#</a>
**`Base.==`** &mdash; *Method*.



```
==(x::ResElem, y::Integer)
```

> Return `true` if $x == y$ arithmetically, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L356' class='documenter-source'>source</a><br>

<a id='Base.==-Tuple{Integer,Nemo.ResElem}' href='#Base.==-Tuple{Integer,Nemo.ResElem}'>#</a>
**`Base.==`** &mdash; *Method*.



```
==(x::Integer, y::ResElem)
```

> Return `true` if $x == y$ arithmetically, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L365' class='documenter-source'>source</a><br>

<a id='Base.==-Tuple{Nemo.ResElem,Nemo.fmpz}' href='#Base.==-Tuple{Nemo.ResElem,Nemo.fmpz}'>#</a>
**`Base.==`** &mdash; *Method*.



```
==(x::ResElem, y::fmpz)
```

> Return `true` if $x == y$ arithmetically, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L374' class='documenter-source'>source</a><br>

<a id='Base.==-Tuple{Nemo.fmpz,Nemo.ResElem}' href='#Base.==-Tuple{Nemo.fmpz,Nemo.ResElem}'>#</a>
**`Base.==`** &mdash; *Method*.



```
==(x::fmpz, y::ResElem)
```

> Return `true` if $x == y$ arithmetically, otherwise return `false`.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L383' class='documenter-source'>source</a><br>


Here are some examples of comparisons.


```
R = ResidueRing(ZZ, 7)
S, x = PolynomialRing(R, "x")
T = ResidueRing(S, x^3 + 3x + 1)

a = R(3)
b = a
c = R(2)
f = T(x^5 + 1)
g = 8f
h = f + g

b == a
c != a
f == g
h != g
isequal(f, g)
a == 3
4 != a
f != 5
```


<a id='Inversion-1'></a>

## Inversion

<a id='Base.inv-Tuple{Nemo.ResElem}' href='#Base.inv-Tuple{Nemo.ResElem}'>#</a>
**`Base.inv`** &mdash; *Method*.



```
inv(a::ResElem)
```

> Return the inverse of the element $a$ in the residue ring. If an impossible inverse is encountered, an exception is raised.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L422' class='documenter-source'>source</a><br>


Here are some examples of computing inverses.


```
R = ResidueRing(ZZ, 49)
S, x = PolynomialRing(R, "x")
T = ResidueRing(S, x^3 + 3x + 1)

a = R(5)
f = T(x^5 + 1)

b = inv(a)
g = inv(f)
```


<a id='Greatest-common-divisor-1'></a>

## Greatest common divisor

<a id='Base.gcd-Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}' href='#Base.gcd-Tuple{Nemo.ResElem{T<:Nemo.RingElem},Nemo.ResElem{T<:Nemo.RingElem}}'>#</a>
**`Base.gcd`** &mdash; *Method*.



```
gcd{T <: RingElem}(a::ResElem{T}, b::ResElem{T})
```

> Return a greatest common divisor of $a$ and $b$ if one exists. This is done by taking the greatest common divisor of the data associated with the supplied residues and taking its greatest common divisor with the modulus.



<a target='_blank' href='https://github.com/wbhart/Nemo.jl/tree/73562614f04fbf543aacd73feb832aff7b4fe899/src/generic/Residue.jl#L460' class='documenter-source'>source</a><br>


Here are some examples of computing a greatest common divisor.


```
R = ResidueRing(ZZ, 43)
S, x = PolynomialRing(R, "x")
T = ResidueRing(S, x^3 - 2x + 1)

a = R(5)
b = R(7)
c = T(x^5 + 1)
d = T(x^2 + x)

f = gcd(a, b)
g = gcd(c, d)
```

