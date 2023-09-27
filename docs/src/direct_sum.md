```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Direct Sums

AbstractAlgebra allows the construction of the external direct sum of any
nonempty vector of finitely presented modules.

Note that external direct sums are considered equal iff they are the same
object.

## Generic direct sum type

AbstractAlgebra provides a generic direct sum type
`Generic.DirectSumModule{T}` where `T` is the element type of the base ring.
The implementation is in `src/generic/DirectSum.jl`

Elements of direct sum modules have type `Generic.DirectSumModuleElem{T}`.

## Abstract types

Direct sum module types belong to the abstract type `FPModule{T}` and their
elements to `FPModuleElem{T}`.

## Constructors

```@docs
direct_sum
```

**Examples**

```jldoctest
julia> F = FreeModule(ZZ, 5)
Free module of rank 5 over integers

julia> m1 = F(BigInt[4, 7, 8, 2, 6])
(4, 7, 8, 2, 6)

julia> m2 = F(BigInt[9, 7, -2, 2, -4])
(9, 7, -2, 2, -4)

julia> S1, f1 = sub(F, [m1, m2])
(Submodule over Integers with 2 generators and no relations, Hom: submodule over Integers with 2 generators and no relations -> free module of rank 5 over integers)

julia> m1 = F(BigInt[3, 1, 7, 7, -7])
(3, 1, 7, 7, -7)

julia> m2 = F(BigInt[-8, 6, 10, -1, 1])
(-8, 6, 10, -1, 1)

julia> S2, f2 = sub(F, [m1, m2])
(Submodule over Integers with 2 generators and no relations, Hom: submodule over Integers with 2 generators and no relations -> free module of rank 5 over integers)

julia> m1 = F(BigInt[2, 4, 2, -3, -10])
(2, 4, 2, -3, -10)

julia> m2 = F(BigInt[5, 7, -6, 9, -5])
(5, 7, -6, 9, -5)

julia> S3, f3 = sub(F, [m1, m2])
(Submodule over Integers with 2 generators and no relations, Hom: submodule over Integers with 2 generators and no relations -> free module of rank 5 over integers)

julia> D, f = direct_sum(S1, S2, S3)
(DirectSumModule over integers, AbstractAlgebra.Generic.ModuleHomomorphism{BigInt}[Hom: submodule over Integers with 2 generators and no relations -> DirectSumModule over integers, Hom: submodule over Integers with 2 generators and no relations -> DirectSumModule over integers, Hom: submodule over Integers with 2 generators and no relations -> DirectSumModule over integers], AbstractAlgebra.Generic.ModuleHomomorphism{BigInt}[Hom: DirectSumModule over integers -> submodule over Integers with 2 generators and no relations, Hom: DirectSumModule over integers -> submodule over Integers with 2 generators and no relations, Hom: DirectSumModule over integers -> submodule over Integers with 2 generators and no relations])
```

## Functionality for direct sums

In addition to the Module interface, AbstractAlgebra direct sums implement the
following functionality.

### Basic manipulation

```@docs
summands(::Generic.DirectSumModule{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> F = FreeModule(ZZ, 5)
Free module of rank 5 over integers

julia> m1 = F(BigInt[4, 7, 8, 2, 6])
(4, 7, 8, 2, 6)

julia> m2 = F(BigInt[9, 7, -2, 2, -4])
(9, 7, -2, 2, -4)

julia> S1, f1 = sub(F, [m1, m2])
(Submodule over Integers with 2 generators and no relations, Hom: submodule over Integers with 2 generators and no relations -> free module of rank 5 over integers)

julia> m1 = F(BigInt[3, 1, 7, 7, -7])
(3, 1, 7, 7, -7)

julia> m2 = F(BigInt[-8, 6, 10, -1, 1])
(-8, 6, 10, -1, 1)

julia> S2, f2 = sub(F, [m1, m2])
(Submodule over Integers with 2 generators and no relations, Hom: submodule over Integers with 2 generators and no relations -> free module of rank 5 over integers)

julia> m1 = F(BigInt[2, 4, 2, -3, -10])
(2, 4, 2, -3, -10)

julia> m2 = F(BigInt[5, 7, -6, 9, -5])
(5, 7, -6, 9, -5)

julia> S3, f3 = sub(F, [m1, m2])
(Submodule over Integers with 2 generators and no relations, Hom: submodule over Integers with 2 generators and no relations -> free module of rank 5 over integers)

julia> D, f = direct_sum(S1, S2, S3)
(DirectSumModule over integers, AbstractAlgebra.Generic.ModuleHomomorphism{BigInt}[Hom: submodule over Integers with 2 generators and no relations -> DirectSumModule over integers, Hom: submodule over Integers with 2 generators and no relations -> DirectSumModule over integers, Hom: submodule over Integers with 2 generators and no relations -> DirectSumModule over integers], AbstractAlgebra.Generic.ModuleHomomorphism{BigInt}[Hom: DirectSumModule over integers -> submodule over Integers with 2 generators and no relations, Hom: DirectSumModule over integers -> submodule over Integers with 2 generators and no relations, Hom: DirectSumModule over integers -> submodule over Integers with 2 generators and no relations])

julia> summands(D)
3-element Vector{AbstractAlgebra.Generic.Submodule{BigInt}}:
 Submodule over Integers with 2 generators and no relations
 Submodule over Integers with 2 generators and no relations
 Submodule over Integers with 2 generators and no relations
```


```
    (D::DirectSumModule{T}(::Vector{<:FPModuleElem{T}}) where T <: RingElement
```

Given a vector (or $1$-dim array) of module elements, where the $i$-th entry
has to be an element of the $i$-summand of $D$, create the corresponding
element in $D$.

**Examples**

```jldoctest
julia> N = FreeModule(QQ, 1);

julia> M = FreeModule(QQ, 2);

julia> D, _ = direct_sum(M, N, M);

julia> D([gen(M, 1), gen(N, 1), gen(M, 2)])
(1//1, 0//1, 1//1, 0//1, 1//1)
```
### Special Homomorphisms

Due to the special structure as direct sums, homomorphisms can be created by specifying
homomorphisms for all summands. In case of the codmain being a direct sum as well,
any homomorphism may be thought of as a matrix containing maps from the $i$-th
source summand to the $j$-th target module:

```
ModuleHomomorphism(D::DirectSumModule{T}, S::DirectSumModule{T}, m::Matrix{Any}) where T <: RingElement
```

Given a matrix $m$ such that the $(i,j)$-th entry is either $0$ (`Int(0)`)
or a `ModuleHomomorphism` from the $i$-th summand of $D$ to the $j$-th summand of
$S$, construct the corresponding homomorphism.

```
ModuleHomomorphism(D::DirectSumModule{T}, S::FPModuleElem{T}, m::Vector{ModuleHomomorphism})
```
Given an array $a$ of `ModuleHomomorphism` such that $a_i$, the $i$-th entry
of $a$ is a `ModuleHomomorphism` from the $i$-th summand of `D` into `S`,
construct the direct sum of the components.


Given a matrix $m$ such that the $(i,j)$-th entry is either $0$ (`Int(0)`)
or a `ModuleHomomorphism` from the $i$-th summand of $D$ to the $j$-th summand of
$S$, construct the corresponding homomorphism.


**Examples**

```jldoctest
julia> N = FreeModule(QQ, 2);

julia> D, _ = direct_sum(N, N);

julia> p = ModuleHomomorphism(N, N, [3,4] .* basis(N));

julia> q = ModuleHomomorphism(N, N, [5,7] .* basis(N));

julia> phi = ModuleHomomorphism(D, D, [p 0; 0 q])
Module homomorphism
  from DirectSumModule over rationals
  to DirectSumModule over rationals

julia> r = ModuleHomomorphism(N, D, [2,3] .* gens(D)[1:2])
Module homomorphism
  from vector space of dimension 2 over rationals
  to DirectSumModule over rationals

julia> psi = ModuleHomomorphism(D, D, [r, r])
Module homomorphism
  from DirectSumModule over rationals
  to DirectSumModule over rationals
```
