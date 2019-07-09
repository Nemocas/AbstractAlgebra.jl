```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Direct Sums

AbstractAlgebra allows the construction of the external direct sum of any
nonempty vector of finitely presented modules.

As well as implementing the entire Module interface, AbstractAlgebra
direct sums also provide the following interface.

Note that external direct sums are considered equal iff they are the same
object.

## Constructors

```@docs
DirectSum(::Vector{<:AbstractAlgebra.FPModule{T}}) where T <: RingElement
DirectSum(::AbstractAlgebra.FPModule{T}...) where T <: RingElement
```

**Examples**

```jldoctest
julia> F = FreeModule(ZZ, 5)
Free module of rank 5 over Integers

julia> m1 = F(BigInt[4, 7, 8, 2, 6])
(4, 7, 8, 2, 6)

julia> m2 = F(BigInt[9, 7, -2, 2, -4])
(9, 7, -2, 2, -4)

julia> S1, f1 = sub(F, [m1, m2])
(Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: Free module of rank 5 over Integers)

julia> m1 = F(BigInt[3, 1, 7, 7, -7])
(3, 1, 7, 7, -7)

julia> m2 = F(BigInt[-8, 6, 10, -1, 1])
(-8, 6, 10, -1, 1)

julia> S2, f2 = sub(F, [m1, m2])
(Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: Free module of rank 5 over Integers)

julia> m1 = F(BigInt[2, 4, 2, -3, -10])
(2, 4, 2, -3, -10)

julia> m2 = F(BigInt[5, 7, -6, 9, -5])
(5, 7, -6, 9, -5)

julia> S3, f3 = sub(F, [m1, m2])
(Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: Free module of rank 5 over Integers)

julia> D, f = DirectSum(S1, S2, S3)
(DirectSumModule over Integers, AbstractAlgebra.Generic.ModuleHomomorphism{BigInt}[Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: DirectSumModule over Integers, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: DirectSumModule over Integers, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: DirectSumModule over Integers], AbstractAlgebra.Generic.ModuleHomomorphism{BigInt}[Module homomorphism with
Domain: DirectSumModule over Integers
Codomain: Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: DirectSumModule over Integers
Codomain: Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: DirectSumModule over Integers
Codomain: Submodule over Integers with 2 generators and no relations
])
```

## Functionality for direct sums

In addition to the Module interface, AbstractAlgebra direct sums implement the
following functionality.

### Basic manipulation

```@docs
summands(::AbstractAlgebra.Generic.DirectSumModule{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> F = FreeModule(ZZ, 5)
Free module of rank 5 over Integers

julia> m1 = F(BigInt[4, 7, 8, 2, 6])
(4, 7, 8, 2, 6)

julia> m2 = F(BigInt[9, 7, -2, 2, -4])
(9, 7, -2, 2, -4)

julia> S1, f1 = sub(F, [m1, m2])
(Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: Free module of rank 5 over Integers)

julia> m1 = F(BigInt[3, 1, 7, 7, -7])
(3, 1, 7, 7, -7)

julia> m2 = F(BigInt[-8, 6, 10, -1, 1])
(-8, 6, 10, -1, 1)

julia> S2, f2 = sub(F, [m1, m2])
(Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: Free module of rank 5 over Integers)

julia> m1 = F(BigInt[2, 4, 2, -3, -10])
(2, 4, 2, -3, -10)

julia> m2 = F(BigInt[5, 7, -6, 9, -5])
(5, 7, -6, 9, -5)

julia> S3, f3 = sub(F, [m1, m2])
(Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: Free module of rank 5 over Integers)

julia> D, f = DirectSum(S1, S2, S3)
(DirectSumModule over Integers, AbstractAlgebra.Generic.ModuleHomomorphism{BigInt}[Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: DirectSumModule over Integers, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: DirectSumModule over Integers, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: DirectSumModule over Integers], AbstractAlgebra.Generic.ModuleHomomorphism{BigInt}[Module homomorphism with
Domain: DirectSumModule over Integers
Codomain: Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: DirectSumModule over Integers
Codomain: Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: DirectSumModule over Integers
Codomain: Submodule over Integers with 2 generators and no relations
])

julia> summands(D)
3-element Array{AbstractAlgebra.Generic.Submodule{BigInt},1}:
 Submodule over Integers with 2 generators and no relations

 Submodule over Integers with 2 generators and no relations

 Submodule over Integers with 2 generators and no relations
```

